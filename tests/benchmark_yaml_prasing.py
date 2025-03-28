#!/usr/bin/env python3
import timeit
import yaml  # Standard PyYAML
import ruamel.yaml as ruyaml
import ryaml

file = "../src/varcall/config/config.yaml"
def benchmark_parser_pyyaml():
    yaml.safe_load(file)

def benchmark_parser_ruyaml():
    ruyaml.YAML().load(file)

def benchmark_parser_ryaml():
    ryaml.loads(file)

if __name__ == "__main__":

    pyyaml = timeit.timeit(stmt=benchmark_parser_pyyaml, number=10000)
    pyyaml = pyyaml * 1000 / 10000

    ru_yaml = timeit.timeit(stmt=benchmark_parser_ruyaml, number=10000)
    ru_yaml = ru_yaml * 1000 / 10000

    r_yaml = timeit.timeit(stmt=benchmark_parser_ryaml, number=10000)
    r_yaml = r_yaml * 1000 / 10000

    print(f"Time taken by pyyaml: {pyyaml} ms")
    print(f"Time taken by ryaml: {ru_yaml} ms")
    print(f"Time taken by ryaml: {r_yaml} ms")
