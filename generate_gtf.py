# coding: utf-8
"""
@author: Xiangwen Ji
@email: jxw01@pku.edu.cn
@file: generate_gtf.py
"""
import os
import re
import argparse
import pandas as pd
from tqdm import tqdm
from tempfile import NamedTemporaryFile
from pybedtools import BedTool

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(
    description="Generate GTF files used by vcf2gsgp.py.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-G",
    "--gtf",
    type=str,
    required=True,
    help="Input GTF file name, gzip format is also supported.",
    dest="input_gtf",
)
parser.add_argument(
    "-O",
    "--output",
    type=str,
    default=os.path.join(CURRENT_DIR, "output"),
    help="The prefix of output files. Two files named '{prefix}.exon.gtf' and '{prefix}.gene.gtf' will be generated.",
    dest="output_prefix",
)
parser.add_argument(
    "--keep-pseudo",
    action="store_true",
    default=False,
    help="Keep pseudogenes.",
    dest="keep_pseudo",
)
parser.add_argument(
    "--add-chr",
    action="store_true",
    default=False,
    help="Add 'chr' prefix before chromosomes.",
    dest="add_chr",
)
parser.add_argument(
    "--tmp",
    "--tmp-dir",
    type=str,
    default=CURRENT_DIR,
    help="The directory of temporary files.",
    dest="tmp_dir",
)
args = parser.parse_args()
args.input_gtf = os.path.abspath(args.input_gtf)
assert os.path.isfile(args.input_gtf), f"--gtf {args.input_gtf} does not exist."
args.output_prefix = os.path.abspath(args.output_prefix)
if args.output_prefix.endswith(("/", "\\")):
    if not os.path.isdir(args.output_prefix):
        os.makedirs(args.output_prefix)
else:
    if not os.path.isdir(output_dir := os.path.dirname(args.output_prefix)):
        os.makedirs(output_dir)
        del output_dir
make_tmp_dir = False
args.tmp_dir = os.path.abspath(args.tmp_dir)
if not os.path.isdir(args.tmp_dir):
    os.makedirs(args.tmp_dir)
    make_tmp_dir = True


def read_gtf(filepath):
    df = pd.read_csv(
        filepath,
        delimiter="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        dtype={0: str},
    )
    df.columns = ["chromosome", "type", "start", "end", "strand", "attributes"]
    attr_pattern = re.compile(r'^(.+) "(.+?)"$')
    df_attr = pd.DataFrame(
        [
            {
                re_result.group(1): re_result.group(2)
                for i in attr.split(";")
                if (re_result := attr_pattern.search(i.strip()))
            }
            for attr in df["attributes"]
        ]
    )
    df[df_attr.columns] = df_attr
    del df["attributes"]
    return df


def generate_gtf():
    print("Read GTF file...")
    df_gtf = read_gtf(args.input_gtf)
    df_gtf = df_gtf[df_gtf["type"].isin(["gene", "exon"])]
    if not args.keep_pseudo:
        df_gtf = df_gtf[~df_gtf["gene_biotype"].str.contains("pseudogene")]
    df_gtf = df_gtf[df_gtf["gene_id"] != ""]
    df_gtf = df_gtf[df_gtf["strand"].isin(["+", "-"])]
    if args.add_chr:
        df_gtf["chromosome"] = "chr" + df_gtf["chromosome"].astype(str)
    exon_tmp = NamedTemporaryFile(
        mode="w+", dir=args.tmp_dir, encoding="utf-8", newline="\n", suffix=".exon.gtf"
    )
    gene_tmp = NamedTemporaryFile(
        mode="w+", dir=args.tmp_dir, encoding="utf-8", newline="\n", suffix=".gene.gtf"
    )
    print("Generate temporary GTF file...")
    f = {
        "gene": gene_tmp,
        "exon": exon_tmp,
    }
    for _, row in tqdm(df_gtf.iterrows(), total=len(df_gtf)):
        line = "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t.\n".format(
            row["chromosome"],
            row["gene_id"],
            row["type"],
            row["start"],
            row["end"],
            row["strand"],
        )
        f[row["type"]].write(line)
    f["gene"].seek(0)
    f["exon"].seek(0)
    print("Sort GTF file...")
    exon_gtf = f"{args.output_prefix}.exon.gtf"
    gene_gtf = f"{args.output_prefix}.gene.gtf"
    with open(exon_gtf, "w", encoding="utf-8", newline="\n") as f_exon, open(
        gene_gtf, "w", encoding="utf-8", newline="\n"
    ) as f_gene:
        f_exon.write(str(BedTool(f["exon"]).sort()))
        f_gene.write(str(BedTool(f["gene"]).sort()))
    f["gene"].close()
    f["exon"].close()
    print("Output files have been writen to {} and {}.".format(exon_gtf, gene_gtf))


if __name__ == "__main__":
    generate_gtf()
    if make_tmp_dir:
        os.removedirs(args.tmp_dir)
