import glob
from typing import List

def find_matching_files(path_pattern: str) -> List[str]:
    """Find all files matching the given wildcard path pattern."""
    return glob.glob(path_pattern)

def main():

    path_pattern = "trial/data/reads/*"
    matching_files = glob.glob(path_pattern)

    matching_files_list = (" ".join(matching_files))

    print(matching_files_list.split())

if __name__ == "__main__":
    main()
