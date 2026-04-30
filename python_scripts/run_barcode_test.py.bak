import argparse
import glob
import json
import subprocess
import sys
import os
from pathlib import Path
from typing import Sequence


def run(cmd: Sequence[str]) -> None:
    """Run a command and stream output; raise if non-zero."""
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True)


def docker_barcode_demux(
    docker_image: str,
    infile: Path,
    out_prefix: Path,
    pattern: Path,
    max_distance: int,
    # *,
    # keep_adapter: bool = False,
    # barcode_file: Path | None = None,
    # common_adapter_file: Path | None = None,
    # save_adapter: bool = True,
) -> None:
    """
    Build and run a docker command for adapter_demux.
    Assumes host paths live under /data and are visible in the container via -v /data:/data.
    """
    # All args must be str
    cmd = [
        "docker", "run", "--rm", "-t",
        "-v", "/data:/data",
        # "--user", "$(id -u):$(id -g)",
        docker_image,
        "barcode_demux",
        "-i", str(infile),
        "-o", str(out_prefix),
        "--barcode", str(pattern),
        "--log-folder", str(out_prefix),
        "--max-distance", str(max_distance),
        "--pipeline-version", "3",
    ]

    run(cmd)


def mode_barcode_test(test_file_dir: Path, docker_image: str) -> None:
    test_files = sorted(Path(test_file_dir).glob("*smc_all_reads.bam"))
    full_length_pattern = Path(test_file_dir) / "ABarcode.fasta"

    for file in test_files:
        print(f"➡️  Running keep_adapter toggles on {file}")
        # prefix = Path(str(file).rsplit(".", maxsplit=1)[0])
        # out_prefix = Path(f"{prefix}.keep_adapter.bam")
        output_dir = test_file_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        docker_barcode_demux(
            docker_image,
            infile=file,
            out_prefix=output_dir,
            pattern=full_length_pattern,
            max_distance=1
        )



def smoke_and_regression_test(test_file_dir: Path, docker_image: str) -> None:
    # single_adapter_dir = Path(test_file_dir) / "single"
    # multiple_adapter_dir = Path(test_file_dir) / "multiple"

    mode_barcode_test(test_file_dir, docker_image)
    # mode_multiple_test(multiple_adapter_dir, docker_image)


def main():
    parser = argparse.ArgumentParser(description="Run SMC smoke/regression tests")
    parser.add_argument("--test-file-dir", required=True, type=Path, dest="test_file_dir")
    parser.add_argument("--docker-image", required=True, dest="docker_image")
    parser.add_argument("--test-type", required=True, type=str, dest="test_type")
    parser.add_argument("--stable-image", default=None, dest="stable_image")
    args = parser.parse_args()

    test_file_dir: Path = args.test_file_dir
    docker_image: str = args.docker_image
    test_type: str = args.test_type.strip().lower()

    try:
        if test_type in {"smoke", "regression"}:
            smoke_and_regression_test(test_file_dir, docker_image)
        else:
            print(f"Unknown test type: {test_type}")
            return
    except Exception as e:
        print(f"Error: {e}")
        # Try to remove the temp image; ignore errors so we don't mask the real cause.
        try:
            run(["docker", "rmi", docker_image])
        except Exception as _:
            pass
        os._exit(1)


if __name__ == "__main__":
    main()