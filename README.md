# varcall

A genome variant calling tool (wrapper) written in python using textual lib

[![PyPI - Version](https://img.shields.io/pypi/v/varcall.svg)](https://pypi.org/project/varcall)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/varcall.svg)](https://pypi.org/project/varcall)

-----

## Table of Contents

- [Installation](#installation)
- [License](#license)

## Installation

> [!NOTE]
> **🚧 Under early development 🚧**

```console
pip install varcall
```
### Contributions are welcomed

> [!NOTE]
> It is recommended to run `varcall` in a conda environment

To install conda refer [miniconda](https://docs.anaconda.com/miniconda/)

`varcall` requires variant calling tools to run, you can download them by:
```console
conda create -n varcall -y
conda install -c bioconda fastqc multiqc bwa samtools bcftools
```

![Home tab image](./doc/images/Hometab-image.jpg?raw=true "Home tab")



![Samtools tab image](./doc/images/Samtoolstab-image.jpg?raw=true "Samtools tab")



![Bcftools tab image](./doc/images/Bcftoolstab-image.jpg?raw=true "Bcftools tab")



![GATK tab image](./doc/images/gatktab-image.jpg?raw=true "GATK tab")



![Pipeline tab image](./doc/images/Pipelinetab-image.jpg?raw=true "Pipeline")



![Help tab image](./doc/images/Helptab-image.jpg?raw=true "Help tab")

## License

`varcall` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
