# VCF2GSGP
The Python implementation of the GSGP algorithm.

Given a VCF file, it calculates the GSGPs for all the samples & genes & signatures.

Contents:
- vcf2gsgp.py: The main process of GSGP. Python>=3.8 on Linux or MacOS system is needed. Detailed options can be obtained by running the program with the -h parameter.
- generate_gtf.py: The Python script used to generate the gene annotation file needed for the above program.
- model: [The COSMIC SBS mutational signatures data files](https://cancer.sanger.ac.uk/signatures/downloads/).
- env.yaml: The conda enviroment dependences.

## install
conda env create -f env.yaml
