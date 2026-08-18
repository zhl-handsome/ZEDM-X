import argparse
import pathlib
import re
import shutil
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    root = pathlib.Path(args.root).resolve()

    if shutil.which("nvidia-smi") is None:
        print("gpu_smoke: SKIP (no NVIDIA driver/device on this machine)")
        return 0

    # Reuse the CPU multiloop regression scenario (steps=1, gravity=0, static
    # contact at step 0) with an isolated output dir so GPU and CPU test
    # artifacts never collide.
    src = (root / "config" / "wall_multiloop_regression.txt").read_text()
    src = re.sub(r"output_dir\s*=.*",
                 "output_dir = build/test_output/gpu_smoke", src)
    cfg = root / "build" / "test_output" / "gpu_smoke.txt"
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text(src)

    completed = subprocess.run(
        [str(exe), "--config", str(cfg)],
        cwd=str(root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=120,
    )
    if completed.returncode != 0 or "done" not in completed.stdout:
        print(completed.stdout)
        raise SystemExit(f"gpu_smoke failed with exit code {completed.returncode}")

    match = re.search(r"step=0 contacts=(\d+)", completed.stdout)
    if not match or int(match.group(1)) < 1:
        print(completed.stdout)
        raise SystemExit("gpu_smoke expected contacts >= 1 at step 0")

    print(f"gpu smoke passed (contacts={match.group(1)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
