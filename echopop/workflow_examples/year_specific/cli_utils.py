import argparse

def get_verbose():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args, unknown = parser.parse_known_args()
    return args.verbose

def get_compare():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compare", action="store_false", help="Compare generated reports to EchoPro"
    )
    args, unknown = parser.parse_known_args()
    return args.compare