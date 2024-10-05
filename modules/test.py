import os.path


def get_basename_and_ext(path) -> tuple:
    filename = os.path.splitext(os.path.basename(path))
    return filename

input_path1 = "/home/vs/github/varcall/python/varcall/"
input_path2 = "/home/vs/github/varcall/python/varcall/trial/results/bam/aligned.index.bam"
input_path_tuple = get_basename_and_ext(input_path2)
last = len(input_path_tuple)-1
output_path = f"{input_path_tuple[0]}.sorted{input_path_tuple[last]}"
print(output_path)
