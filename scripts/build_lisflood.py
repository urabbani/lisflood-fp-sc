#!/usr/bin/env python3
"""
Build LISFLOOD-FP from integrated source.

Usage:
    python scripts/build_lisflood.py              # CPU-only build
    python scripts/build_lisflood.py --cuda        # GPU build (requires CUDA toolkit)
    python scripts/build_lisflood.py --force       # Rebuild even if binary exists
    python scripts/build_lisflood.py --clean       # Clean build directory first
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = PROJECT_ROOT / "src" / "lisflood-fp"
BUILD_DIR = SOURCE_DIR / "build"


def find_binary() -> Path | None:
    """Find compiled LISFLOOD-FP binary."""
    candidates = [
        BUILD_DIR / "lisflood",
        BUILD_DIR / "lisflood.exe",
        BUILD_DIR / "Release" / "lisflood.exe",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def check_dependencies() -> list[str]:
    """Check for required build tools. Returns list of missing."""
    missing = []

    # Check cmake
    try:
        result = subprocess.run(["cmake", "--version"], capture_output=True, text=True)
        if result.returncode != 0:
            missing.append("cmake")
    except FileNotFoundError:
        missing.append("cmake")

    # Check C++ compiler
    system = platform.system()
    if system == "Linux" or system == "Darwin":
        for compiler in ["g++", "c++"]:
            try:
                subprocess.run([compiler, "--version"], capture_output=True, text=True)
                break
            except FileNotFoundError:
                continue
        else:
            missing.append("g++ or compatible C++ compiler")
    elif system == "Windows":
        # MSVC is typically found by cmake
        pass

    return missing


def build(clean: bool = False, cuda: bool = False, jobs: int | None = None) -> Path:
    """
    Build LISFLOOD-FP from source.

    Args:
        clean: Remove build directory first
        cuda: Enable CUDA GPU acceleration
        jobs: Number of parallel jobs (default: CPU count)

    Returns:
        Path to compiled binary

    Raises:
        RuntimeError: If build fails
    """
    if not SOURCE_DIR.exists():
        raise FileNotFoundError(
            f"LISFLOOD-FP source not found at {SOURCE_DIR}\n"
            f"Source should be in src/lisflood-fp/"
        )

    # Check dependencies
    missing = check_dependencies()
    if missing:
        raise RuntimeError(
            f"Missing build dependencies: {', '.join(missing)}\n"
            f"Install them before building."
        )

    # Clean if requested
    if clean and BUILD_DIR.exists():
        print(f"Cleaning build directory: {BUILD_DIR}")
        shutil.rmtree(BUILD_DIR)

    # Create build directory
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # Determine config file
    config_dir = VENDOR_DIR / "config"
    system = platform.system()

    if system == "Linux" or system == "Darwin":
        config_file = config_dir / "local"
    else:
        config_file = config_dir / "local"

    # Check if config exists, fall back to creating one
    if not config_file.exists():
        config_file = config_dir / "local"
        config_file.parent.mkdir(parents=True, exist_ok=True)
        config_file.write_text("# Auto-generated local config\n")
        print(f"Created default config: {config_file}")

    # CMake configure
    print("\n--- Configuring LISFLOOD-FP ---")
    cmake_args = [
        "cmake", str(VENDOR_DIR),
        f"-D_CONFIG=../config/local",
    ]

    if cuda:
        cmake_args.append("-D_CUDA=1")
        print("   CUDA GPU acceleration: ENABLED")
    else:
        print("   CUDA GPU acceleration: disabled (CPU-only build)")

    print(f"   Build directory: {BUILD_DIR}")
    print(f"   Platform: {system}")

    result = subprocess.run(
        cmake_args,
        cwd=BUILD_DIR,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"\nCMake configure failed:\n{result.stderr}")
        raise RuntimeError(f"CMake configure failed. See output above.")

    print("   CMake configure: OK")

    # Build
    print("\n--- Building LISFLOOD-FP ---")

    if jobs is None:
        try:
            jobs = os.cpu_count() or 4
        except Exception:
            jobs = 4

    make_result = subprocess.run(
        ["make", f"-j{jobs}"],
        cwd=BUILD_DIR,
        capture_output=True,
        text=True,
    )

    if make_result.returncode != 0:
        print(f"\nBuild failed:\n{make_result.stderr}")
        raise RuntimeError(f"Build failed. See output above.")

    # Find the binary
    binary = find_binary()
    if binary is None:
        raise RuntimeError(
            f"Build completed but binary not found in {BUILD_DIR}\n"
            f"Check the build output above."
        )

    # Make executable (Unix)
    if platform.system() != "Windows":
        binary.chmod(binary.stat().st_mode | 0o755)

    print(f"\n--- Build successful ---")
    print(f"   Binary: {binary}")
    print(f"   Size: {binary.stat().st_size / 1024 / 1024:.1f} MB")

    return binary


def main():
    parser = argparse.ArgumentParser(
        description="Build LISFLOOD-FP from vendored source",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--cuda", action="store_true",
        help="Enable CUDA GPU acceleration (requires CUDA toolkit)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Rebuild even if binary already exists",
    )
    parser.add_argument(
        "--clean", action="store_true",
        help="Clean build directory before building",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=None,
        help="Number of parallel build jobs (default: CPU count)",
    )

    args = parser.parse_args()

    # Check if binary already exists
    if not args.force and not args.clean:
        existing = find_binary()
        if existing:
            print(f"LISFLOOD-FP binary already exists: {existing}")
            print("Use --force to rebuild or --clean for a fresh build.")
            return 0

    try:
        binary = build(clean=args.clean, cuda=args.cuda, jobs=args.jobs)
        print(f"\nReady to use. Set in config:")
        print(f'  executable: "{binary}"')
        print(f'  # or use:')
        print(f'  executable: "auto"')
        return 0
    except (FileNotFoundError, RuntimeError) as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
