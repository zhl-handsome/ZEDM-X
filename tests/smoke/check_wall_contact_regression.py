import argparse
import pathlib
import re
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    root = pathlib.Path(args.root).resolve()
    config = "config/wall_regression.txt"
    completed = subprocess.run(
        [str(exe), "--config", config],
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        raise SystemExit(f"{config} failed with exit code {completed.returncode}")

    match = re.search(r"step=0 contacts=(\d+)", completed.stdout)
    if not match:
        print(completed.stdout)
        raise SystemExit("missing step=0 contacts output")

    contacts = int(match.group(1))
    if contacts != 1:
        print(completed.stdout)
        raise SystemExit(f"expected a single aggregated wall contact, got {contacts}")

    print("wall contact regression passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
