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
    completed = subprocess.run(
        [str(exe), "--config", "config/wall_multiloop_regression.txt"],
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        raise SystemExit(f"wall_multiloop_regression failed with exit code {completed.returncode}")

    if "section_loops=2" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit("expected wall debug output to report section_loops=2")

    match = re.search(r"step=0 contacts=(\d+)", completed.stdout)
    if not match:
        print(completed.stdout)
        raise SystemExit("missing step=0 contacts output")
    if int(match.group(1)) != 1:
        print(completed.stdout)
        raise SystemExit("expected multiple wall loops to be simplified to one force contact")

    print("wall multiloop regression passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
