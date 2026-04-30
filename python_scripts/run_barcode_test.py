import argparse
import difflib
import hashlib
import os
import subprocess
from pathlib import Path
from typing import Sequence


FASTQ_SUFFIXES = {".fastq", ".fq"}


def run(cmd: Sequence[str]) -> None:
    """Run a command and stream output; raise if non-zero."""
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def collect_fastq_files(output_dir: Path) -> dict[Path, Path]:
    """
    Return {relative_path: absolute_path} for all FASTQ files under output_dir.
    """
    result: dict[Path, Path] = {}

    for p in output_dir.rglob("*"):
        if not p.is_file():
            continue

        if p.suffix.lower() in FASTQ_SUFFIXES:
            rel = p.relative_to(output_dir)
            result[rel] = p

    return result


def print_text_diff_preview(a: Path, b: Path, max_lines: int = 200) -> None:
    """
    Print a small unified diff preview for debugging.
    Avoid reading huge FASTQ files completely.
    """
    try:
        with a.open("r", encoding="utf-8", errors="replace") as fa:
            a_lines = []
            for _, line in zip(range(max_lines), fa):
                a_lines.append(line)

        with b.open("r", encoding="utf-8", errors="replace") as fb:
            b_lines = []
            for _, line in zip(range(max_lines), fb):
                b_lines.append(line)

        diff = difflib.unified_diff(
            a_lines,
            b_lines,
            fromfile=str(a),
            tofile=str(b),
            lineterm="",
        )

        print("---- diff preview ----")
        for line in diff:
            print(line.rstrip("\n"))
        print("---- end diff preview ----")

    except Exception as e:
        print(f"Failed to print diff preview: {e}")


def assert_fastq_dirs_equal(
    current_dir: Path,
    stable_dir: Path,
) -> None:
    """
    Compare FASTQ files under two output directories.

    Checks:
    1. The relative FASTQ file list is identical.
    2. Every corresponding FASTQ file has identical sha256.
    """
    current_files = collect_fastq_files(current_dir)
    stable_files = collect_fastq_files(stable_dir)

    current_set = set(current_files.keys())
    stable_set = set(stable_files.keys())

    missing_in_current = sorted(stable_set - current_set)
    extra_in_current = sorted(current_set - stable_set)

    if missing_in_current or extra_in_current:
        print("FASTQ file list mismatch.")

        if missing_in_current:
            print("\nFiles missing in current output:")
            for p in missing_in_current:
                print(f"  - {p}")

        if extra_in_current:
            print("\nExtra files in current output:")
            for p in extra_in_current:
                print(f"  - {p}")

        raise AssertionError("FASTQ file list differs between current and stable outputs")

    if not current_files:
        raise AssertionError(f"No FASTQ files found under {current_dir}")

    mismatch_count = 0

    for rel in sorted(current_set):
        current_file = current_files[rel]
        stable_file = stable_files[rel]

        current_hash = sha256_file(current_file)
        stable_hash = sha256_file(stable_file)

        if current_hash != stable_hash:
            mismatch_count += 1

            print(f"\nFASTQ content mismatch: {rel}")
            print(f"  current: {current_file}")
            print(f"  stable : {stable_file}")
            print(f"  current sha256: {current_hash}")
            print(f"  stable  sha256: {stable_hash}")

            print_text_diff_preview(current_file, stable_file)

    if mismatch_count > 0:
        raise AssertionError(f"{mismatch_count} FASTQ file(s) differ")

    print(f"✅ FASTQ no-diff passed: {len(current_files)} file(s) are identical")


def docker_barcode_demux(
    docker_image: str,
    infile: Path,
    output_dir: Path,
    pattern: Path,
    max_distance: int,
    pipeline_version: int,
    output_format: str,
) -> None:
    """
    Build and run a docker command for barcode_demux.

    Assumption:
    Host paths live under /data and are visible in the container via:
        -v /data:/data
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "docker", "run", "--rm", "-t",
        "-v", "/data:/data",
        docker_image,
        "barcode_demux",
        "-i", str(infile),
        "-o", str(output_dir),
        "--barcode", str(pattern),
        "--log-folder", str(output_dir),
        "--max-distance", str(max_distance),
        "--pipeline-version", str(pipeline_version),
        "--output-format", output_format,
    ]

    run(cmd)


def run_pipeline1_fastq_no_diff_for_file(
    infile: Path,
    test_file_dir: Path,
    current_image: str,
    stable_image: str,
    pattern: Path,
    max_distance: int,
) -> None:
    """
    Run current image and stable image on one BAM,
    then compare generated FASTQ files.
    """
    sample_name = infile.name
    if sample_name.endswith(".bam"):
        sample_name = sample_name[:-4]

    base_output_dir = test_file_dir / "output" / "no_diff_pipeline1_fastq" / sample_name

    current_output_dir = base_output_dir / "current"
    stable_output_dir = base_output_dir / "stable"

    print(f"\n➡️  Running current image for pipeline 1 FASTQ: {infile}")
    docker_barcode_demux(
        docker_image=current_image,
        infile=infile,
        output_dir=current_output_dir,
        pattern=pattern,
        max_distance=max_distance,
        pipeline_version=1,
        output_format="fastq",
    )

    print(f"\n➡️  Running stable image for pipeline 1 FASTQ: {infile}")
    docker_barcode_demux(
        docker_image=stable_image,
        infile=infile,
        output_dir=stable_output_dir,
        pattern=pattern,
        max_distance=max_distance,
        pipeline_version=1,
        output_format="fastq",
    )

    print(f"\n🔍 Comparing FASTQ outputs for {infile.name}")
    assert_fastq_dirs_equal(
        current_dir=current_output_dir,
        stable_dir=stable_output_dir,
    )


def mode_barcode_test(
    test_file_dir: Path,
    docker_image: str,
    pipeline_version: int = 3,
    output_format: str = "bam",
    max_distance: int = 1,
) -> None:
    """
    Basic smoke test: only run current image.
    """
    test_files = sorted(test_file_dir.glob("*smc_all_reads.bam"))
    pattern = test_file_dir / "ABarcode.fasta"

    if not test_files:
        raise FileNotFoundError(f"No *smc_all_reads.bam found under {test_file_dir}")

    if not pattern.exists():
        raise FileNotFoundError(f"Barcode pattern file not found: {pattern}")

    output_dir = test_file_dir / "output" / "smoke"

    for file in test_files:
        print(f"\n➡️  Running barcode_demux smoke test on {file}")

        sample_output_dir = output_dir / file.stem

        docker_barcode_demux(
            docker_image=docker_image,
            infile=file,
            output_dir=sample_output_dir,
            pattern=pattern,
            max_distance=max_distance,
            pipeline_version=pipeline_version,
            output_format=output_format,
        )


def mode_pipeline1_fastq_no_diff_test(
    test_file_dir: Path,
    docker_image: str,
    stable_image: str,
    max_distance: int = 1,
) -> None:
    """
    Regression no-diff test:
    current image vs stable image,
    pipeline 1,
    output format FASTQ.
    """
    test_files = sorted(test_file_dir.glob("*smc_all_reads.bam"))
    pattern = test_file_dir / "ABarcode.fasta"

    if not test_files:
        raise FileNotFoundError(f"No *smc_all_reads.bam found under {test_file_dir}")

    if not pattern.exists():
        raise FileNotFoundError(f"Barcode pattern file not found: {pattern}")

    for file in test_files:
        run_pipeline1_fastq_no_diff_for_file(
            infile=file,
            test_file_dir=test_file_dir,
            current_image=docker_image,
            stable_image=stable_image,
            pattern=pattern,
            max_distance=max_distance,
        )


def smoke_and_regression_test(
    test_file_dir: Path,
    docker_image: str,
    stable_image: str | None,
) -> None:
    """
    Recommended behavior:
    1. Run smoke test with current image.
    2. If stable image is provided, also run pipeline 1 FASTQ no-diff test.
    """
    mode_barcode_test(
        test_file_dir=test_file_dir,
        docker_image=docker_image,
        pipeline_version=3,
        output_format="bam",
        max_distance=1,
    )

    if stable_image is not None:
        mode_pipeline1_fastq_no_diff_test(
            test_file_dir=test_file_dir,
            docker_image=docker_image,
            stable_image=stable_image,
            max_distance=1,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run SMC smoke/regression tests")

    parser.add_argument(
        "--test-file-dir",
        required=True,
        type=Path,
        dest="test_file_dir",
    )

    parser.add_argument(
        "--docker-image",
        required=True,
        dest="docker_image",
        help="Current image to test",
    )

    parser.add_argument(
        "--stable-image",
        default=None,
        dest="stable_image",
        help="Previous packaged stable image for no-diff regression test",
    )

    parser.add_argument(
        "--test-type",
        required=True,
        type=str,
        dest="test_type",
        choices=["smoke", "regression", "no-diff"],
    )

    args = parser.parse_args()

    test_file_dir: Path = args.test_file_dir
    docker_image: str = args.docker_image
    stable_image: str | None = args.stable_image
    test_type: str = args.test_type.strip().lower()

    try:
        if test_type == "smoke":
            mode_barcode_test(
                test_file_dir=test_file_dir,
                docker_image=docker_image,
                pipeline_version=3,
                output_format="bam",
                max_distance=1,
            )

        elif test_type == "regression":
            if stable_image is None:
                raise ValueError("--stable-image is required for regression test")

            smoke_and_regression_test(
                test_file_dir=test_file_dir,
                docker_image=docker_image,
                stable_image=stable_image,
            )

        elif test_type == "no-diff":
            if stable_image is None:
                raise ValueError("--stable-image is required for no-diff test")

            mode_pipeline1_fastq_no_diff_test(
                test_file_dir=test_file_dir,
                docker_image=docker_image,
                stable_image=stable_image,
                max_distance=1,
            )

        else:
            raise ValueError(f"Unknown test type: {test_type}")

    except Exception as e:
        print(f"\n❌ Error: {e}")

        try:
            run(["docker", "rmi", docker_image])
        except Exception:
            pass

        sys_exit_code = 1
        os._exit(sys_exit_code)


if __name__ == "__main__":
    main()