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

    # v8 unified per-vertex contact (5c9595e) replaced the per-loop section
    # telemetry with one wall_contact debug line: multiple contained vertices
    # (n_cross >= 2) collapsing into a single force contact.
    match = re.search(r"wall_contact: particle=\d+ wall=\d+ n_cross=(\d+)", completed.stdout)
    if not match:
        print(completed.stdout)
        raise SystemExit("expected wall contact debug output with n_cross count")
    if int(match.group(1)) < 2:
        print(completed.stdout)
        raise SystemExit("expected multiple contained vertices in the wall contact")

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
