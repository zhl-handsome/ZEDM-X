import argparse
import pathlib
import subprocess
import sys


def run_case(exe: pathlib.Path, root: pathlib.Path, config: str) -> None:
    cmd = [str(exe), "--config", config]
    completed = subprocess.run(
        cmd,
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        raise SystemExit(f"{config} failed with exit code {completed.returncode}")
    if "done" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit(f"{config} did not print done")
    if "step=" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit(f"{config} did not print step progress")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    root = pathlib.Path(args.root).resolve()
    cases = [
        "config/smoke_no_contact.txt",
        "config/smoke_wall_contact.txt",
    ]
    for case in cases:
        run_case(exe, root, case)
    print("smoke tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
