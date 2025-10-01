import argparse

def get_verbose():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args, unknown = parser.parse_known_args()
    return args.verbose