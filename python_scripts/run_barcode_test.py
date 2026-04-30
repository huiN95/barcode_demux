import argparse
import difflib
import functools
import hashlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence


def run(cmd: Sequence[str]) -> None:
    """Run a command and stream output; raise if non-zero."""
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_capture(cmd: Sequence[str]) -> str:
    """
    Run command and capture stdout + stderr.
    Do not raise on non-zero, because some --help implementations exit non-zero.
    """
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return (p.stdout or "") + "\n" + (p.stderr or "")


def is_fastq_file(path: Path) -> bool:
    name = path.name.lower()
    return (
        name.endswith(".fastq")
        or name.endswith(".fq")
        or name.endswith(".fastq.gz")
        or name.endswith(".fq.gz")
    )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)

    return h.hexdigest()


def collect_fastq_files(output_dir: Path) -> Dict[Path, Path]:
    """
    Return:
        {relative_path: absolute_path}

    Only FASTQ files are collected.
    Logs, metrics, json, bam, bai, etc. are ignored.
    """
    result: Dict[Path, Path] = {}

    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")

    for p in output_dir.rglob("*"):
        if not p.is_file():
            continue

        if is_fastq_file(p):
            rel = p.relative_to(output_dir)
            result[rel] = p

    return result


def print_text_diff_preview(a: Path, b: Path, max_lines: int = 200) -> None:
    """
    Print a small unified diff preview for debugging.

    For .gz files, only hash mismatch is reported.
    """
    if a.name.lower().endswith(".gz") or b.name.lower().endswith(".gz"):
        print("Skip text diff preview for gzipped FASTQ file.")
        return

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


def assert_fastq_dirs_equal(current_dir: Path, stable_dir: Path) -> None:
    """
    Compare FASTQ files under two output directories.

    Checks:
      1. Relative FASTQ file list is identical.
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
        raise AssertionError(f"No FASTQ files found under current output: {current_dir}")

    mismatch_count = 0

    for rel in sorted(current_set):
        current_file = current_files[rel]
        stable_file = stable_files[rel]

        current_hash = sha256_file(current_file)
        stable_hash = sha256_file(stable_file)

        if current_hash != stable_hash:
            mismatch_count += 1

            print(f"\nFASTQ content mismatch: {rel}")
            print(f"  current file = {current_file}")
            print(f"  stable file  = {stable_file}")
            print(f"  current sha256 = {current_hash}")
            print(f"  stable  sha256 = {stable_hash}")

            print_text_diff_preview(current_file, stable_file)

    if mismatch_count > 0:
        raise AssertionError(f"{mismatch_count} FASTQ file(s) differ")

    print(f"✅ FASTQ no-diff passed: {len(current_files)} file(s) are identical")


def choose_flag(help_text: str, kebab: str, snake: str) -> str:
    """
    Choose CLI flag style according to docker image help output.

    Old stable image may use:
        --log_folder
        --max_distance
        --pipeline_version
        --output_format

    New image may use:
        --log-folder
        --max-distance
        --pipeline-version
        --output-format
    """
    if snake in help_text and kebab not in help_text:
        return snake

    if kebab in help_text:
        return kebab

    if snake in help_text:
        return snake

    # Default to new style if not found in help.
    return kebab


@functools.lru_cache(maxsize=16)
def detect_cli_flags(docker_image: str) -> Dict[str, str]:
    """
    Detect which CLI flag style the image supports.

    This allows current image and stable image to use different argument styles.
    """
    help_cmd = [
        "docker",
        "run",
        "--rm",
        docker_image,
        "barcode_demux",
        "--help",
    ]

    help_text = run_capture(help_cmd)

    flags = {
        "log_folder": choose_flag(help_text, "--log-folder", "--log_folder"),
        "max_distance": choose_flag(help_text, "--max-distance", "--max_distance"),
        "pipeline_version": choose_flag(help_text, "--pipeline-version", "--pipeline_version"),
        "output_format": choose_flag(help_text, "--output-format", "--output_format"),
    }

    print(f"Detected CLI flags for image {docker_image}:")
    for k, v in flags.items():
        print(f"  {k} = {v}")

    return flags


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
    Run barcode_demux in docker.

    Important:
      - /data is mounted read-only for test input.
      - output_dir is mounted separately and is writable.
      - output_dir should usually be under CI_PROJECT_DIR, not under test data dir.
      - CLI flag style is detected per image for stable/new compatibility.
    """
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    flags = detect_cli_flags(docker_image)

    cmd = [
        "docker",
        "run",
        "--rm",

        # Input data. Read-only is safer for CI datasets.
        "-v",
        "/data:/data:ro",

        # Writable output directory.
        "-v",
        f"{str(output_dir)}:{str(output_dir)}",

        # Avoid root-owned files in GitLab CI workspace.
        "--user",
        f"{os.getuid()}:{os.getgid()}",

        docker_image,
        "barcode_demux",

        "-i",
        str(infile),

        "-o",
        str(output_dir),

        "--barcode",
        str(pattern),

        flags["log_folder"],
        str(output_dir),

        flags["max_distance"],
        str(max_distance),

        flags["pipeline_version"],
        str(pipeline_version),

        flags["output_format"],
        output_format,
    ]

    run(cmd)


def mode_barcode_test(
    test_file_dir: Path,
    output_root: Path,
    docker_image: str,
    pipeline_version: int = 3,
    output_format: str = "bam",
    max_distance: int = 1,
) -> None:
    """
    Basic smoke test.

    Only runs the current image.
    """
    test_files = sorted(test_file_dir.glob("*smc_all_reads.bam"))
    pattern = test_file_dir / "ABarcode.fasta"

    if not test_files:
        raise FileNotFoundError(f"No *smc_all_reads.bam found under {test_file_dir}")

    if not pattern.exists():
        raise FileNotFoundError(f"Barcode pattern file not found: {pattern}")

    output_dir = output_root / "smoke"

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


def run_pipeline1_fastq_no_diff_for_file(
    infile: Path,
    output_root: Path,
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

    base_output_dir = output_root / "no_diff_pipeline1_fastq" / sample_name

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


def mode_pipeline1_fastq_no_diff_test(
    test_file_dir: Path,
    output_root: Path,
    docker_image: str,
    stable_image: str,
    max_distance: int = 1,
) -> None:
    """
    No-diff regression test.

    Compares:
      current image pipeline 1 FASTQ output
      vs
      stable image pipeline 1 FASTQ output
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
            output_root=output_root,
            current_image=docker_image,
            stable_image=stable_image,
            pattern=pattern,
            max_distance=max_distance,
        )


def smoke_and_regression_test(
    test_file_dir: Path,
    output_root: Path,
    docker_image: str,
    stable_image: Optional[str],
) -> None:
    """
    Recommended regression behavior:
      1. Run smoke test with current image.
      2. If stable image is provided, run pipeline 1 FASTQ no-diff test.
    """
    mode_barcode_test(
        test_file_dir=test_file_dir,
        output_root=output_root,
        docker_image=docker_image,
        pipeline_version=3,
        output_format="bam",
        max_distance=1,
    )

    if stable_image is not None:
        mode_pipeline1_fastq_no_diff_test(
            test_file_dir=test_file_dir,
            output_root=output_root,
            docker_image=docker_image,
            stable_image=stable_image,
            max_distance=1,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run barcode_demux smoke/regression tests")

    parser.add_argument(
        "--test-file-dir",
        required=True,
        type=Path,
        dest="test_file_dir",
        help="Directory containing test BAM files and ABarcode.fasta",
    )

    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        dest="output_root",
        help="Writable root directory for test outputs",
    )

    parser.add_argument(
        "--docker-image",
        required=True,
        dest="docker_image",
        help="Current docker image to test",
    )

    parser.add_argument(
        "--stable-image",
        default=None,
        dest="stable_image",
        help="Previous stable docker image for no-diff regression test",
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
    output_root: Path = args.output_root
    docker_image: str = args.docker_image
    stable_image: Optional[str] = args.stable_image
    test_type: str = args.test_type.strip().lower()

    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    try:
        if test_type == "smoke":
            mode_barcode_test(
                test_file_dir=test_file_dir,
                output_root=output_root,
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
                output_root=output_root,
                docker_image=docker_image,
                stable_image=stable_image,
            )

        elif test_type == "no-diff":
            if stable_image is None:
                raise ValueError("--stable-image is required for no-diff test")

            mode_pipeline1_fastq_no_diff_test(
                test_file_dir=test_file_dir,
                output_root=output_root,
                docker_image=docker_image,
                stable_image=stable_image,
                max_distance=1,
            )

        else:
            raise ValueError(f"Unknown test type: {test_type}")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()