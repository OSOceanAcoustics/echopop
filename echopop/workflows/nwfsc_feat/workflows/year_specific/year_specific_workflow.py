# CLI hook
import argparse
import runpy
import sys
import os

# Set up argument parser and parsing
parser = argparse.ArgumentParser()
parser.add_argument("--year", required=True, help="Which workflow script to run (e.g. hake_1995)")
parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
args, unknown = parser.parse_known_args()

# Pass additional CLI args to the target script
sys.argv = [args.year + ".py"] + unknown
if args.verbose:
    sys.argv.append("--verbose")

# Get the script path
script_path = os.path.join(os.path.dirname(__file__), f"{args.year}.py")
runpy.run_path(script_path, run_name="__main__")