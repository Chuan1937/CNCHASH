"""
Command-line interface for CNCHASH.
"""

import argparse
import os
import sys

from . import __version__, driver, io


def _load_velocity_model(filename):
    depth, velocity = io.read_velocity_model(filename)
    return depth, velocity


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="CNCHASH - Earthquake focal mechanism inversion (Python version of HASH v1.2)"
    )

    parser.add_argument("input_file", help="HASH input file (like example.inp)")

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "--velocity",
        action="append",
        metavar="FILE",
        help=(
            "1D velocity model file (depth velocity, two columns). "
            "Repeat for multiple models (rotated per Monte Carlo trial). "
            "Default: models listed at the end of the input file."
        ),
    )
    parser.add_argument("--version", action="version", version=f"CNCHASH {__version__}")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)

    velocity_models = None
    if args.velocity:
        velocity_models = []
        for filename in args.velocity:
            if not os.path.exists(filename):
                print(f"Error: Velocity model file '{filename}' not found")
                sys.exit(1)
            velocity_models.append(_load_velocity_model(filename))

    # Run HASH
    print(f"Running HASH on input file: {args.input_file}")

    try:
        results = driver.run_hash_from_file(args.input_file, velocity_models=velocity_models)

        print(f"\nProcessed {len(results)} events")

        # Print summary
        n_success = sum(1 for r in results if r.get("success", False))
        n_failed = len(results) - n_success

        print(f"  Successful: {n_success}")
        print(f"  Failed: {n_failed}")

        if args.verbose:
            print("\nEvent results:")
            for i, result in enumerate(results):
                if result.get("success"):
                    print(
                        f"  Event {i + 1}: "
                        f"strike={result.get('strike_avg', 0):.1f}, "
                        f"dip={result.get('dip_avg', 0):.1f}, "
                        f"rake={result.get('rake_avg', 0):.1f}, "
                        f"quality={result.get('quality', '?')}"
                    )
                else:
                    print(f"  Event {i + 1}: Failed ({result.get('quality', 'F')})")

    except Exception as e:
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
