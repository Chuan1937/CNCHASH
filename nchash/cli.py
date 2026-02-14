"""
Command-line interface for NCHASH.
"""

import argparse
import os
import sys

from . import driver


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="NCHASH - Earthquake focal mechanism inversion (Python version of HASH v1.2)"
    )

    parser.add_argument("input_file", help="HASH input file (like example.inp)")

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    parser.add_argument("--version", action="version", version="NCHASH 1.0.0")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)

    # Run HASH
    print(f"Running HASH on input file: {args.input_file}")

    try:
        results = driver.run_hash_from_file(args.input_file)

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
