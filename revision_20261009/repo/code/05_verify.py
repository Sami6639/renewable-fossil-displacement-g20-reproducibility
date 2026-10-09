import json

from pipeline import verify_outputs


if __name__ == "__main__":
    print(json.dumps(verify_outputs(), indent=2))
