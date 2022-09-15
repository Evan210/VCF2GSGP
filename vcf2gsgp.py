# coding: utf-8
"""
@author: Xiangwen Ji
@email: jxw01@pku.edu.cn
@file: vcf2gsgp.py
@date: 2021/12/04
@desc: Calculate gene somatic genome patterns (GSGPs) using VCF file
"""
import os
import pysam
import argparse
import traceback
import logging.config
import pandas as pd
from math import ceil
from datetime import datetime
from time import sleep, time
from gzip import GzipFile
from multiprocessing import Pool, Process, Manager
from pybedtools import BedTool
from collections import OrderedDict
from tempfile import NamedTemporaryFile, mkdtemp
from io import StringIO
from sklearn.decomposition import NMF
from shutil import rmtree

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_FILE = os.path.join(
    CURRENT_DIR,
    "vcf2gsgp_{}.log".format(datetime.now().strftime("%Y%m%d%H%M%S")),
)

parser = argparse.ArgumentParser(
    description="Calculate gene somatic genome patterns (GSGPs) using VCF file.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-I",
    "--input",
    type=str,
    required=True,
    help="Input VCF file name, gzip format is also supported.",
    dest="input",
)
parser.add_argument(
    "-R",
    "--ref",
    "--reference",
    type=str,
    required=True,
    help="The path of reference FASTA.",
    dest="reference",
)
parser.add_argument(
    "-G",
    "--gtf-prefix",
    type=str,
    required=True,
    help="The prefix of GTF file generated by 'generate_gtf.py'.",
    dest="gtf_prefix",
)
parser.add_argument(
    "-O",
    "--output",
    "--output-dir",
    type=str,
    required=True,
    help="The directory of output files.",
    dest="output_dir",
)
parser.add_argument(
    "-g",
    "--genome",
    type=str,
    choices=["GRCh37", "GRCh38", "mm9", "mm10", "rn6"],
    default="GRCh38",
    help="The genome of the COSMIC SBS signature.",
    dest="genome",
)
parser.add_argument(
    "-t", "--threads", type=int, default=4, help="Number of processes.", dest="threads"
)
parser.add_argument(
    "-d",
    "--deep",
    type=int,
    default=None,
    help="The minimum value of read depth. No filtering by default.",
    dest="deep",
)
parser.add_argument(
    "-f",
    "--vaf",
    type=float,
    default=None,
    nargs="+",
    help="The interval of variant allele frequency (VAF), every two numbers is an interval. "
    + "Only the mutations with VAF between one of the given intervals are retained. "
    + "Given only one value, it means the minimum value of VAF. No filtering by default.",
    dest="vaf",
)
parser.add_argument(
    "--up",
    "--upstream",
    type=int,
    default=121925,
    help="The maximum base pairs that allows the variant to be located upstream of the gene.",
    dest="upstream",
)
parser.add_argument(
    "--down",
    "--downstream",
    type=int,
    default=121388,
    help="The maximum base pairs that allows the variant to be located downstream of the gene.",
    dest="downstream",
)
parser.add_argument(
    "-w",
    "--weight",
    type=float,
    default=None,
    nargs="+",
    help="The weights of variant types when a variant can be annotated to multi-genes. "
    + "Default (None) or given only one value, different variant types will have the same weights. "
    + "Given two values, the first one is the weight of exon or intron variants, the second one is the weight of upstream/downstream variants. "
    + "Given three values, 1st: exon, 2nd: intron, 3rd: upstream/downstream. "
    + "Given four values or more, 1st: exon, 2nd: intron, 3rd: upstream, 4th: downstream, the rest: ignored.",
    dest="weight",
)
parser.add_argument(
    "--max-iter",
    type=int,
    default=1000,
    help="Maximum number of iterations, used in NMF.",
    dest="max_iter",
)
parser.add_argument(
    "--random-seed",
    type=int,
    default=None,
    help="Random seed of NMF.",
    dest="random_seed",
)
parser.add_argument(
    "--add-chr",
    action="store_true",
    default=False,
    help="Add 'chr' before chromosome. "
    + "This option is used to handle the inconsistency between hg19 and GRCh37 (hg38 and GRCh38).",
    dest="add_chr",
)
parser.add_argument(
    "--save-x",
    action="store_true",
    default=False,
    help="Save the matrix X used for NMF.",
    dest="save_x",
)
parser.add_argument(
    "--log",
    type=str,
    default=LOG_FILE,
    help="The path of log file.",
    dest="log",
)
parser.add_argument(
    "--silent",
    action="store_true",
    default=False,
    help="Do not output logs to the console. This will not cancel the log file set by --log.",
    dest="silent",
)
parser.add_argument(
    "--tmp",
    "--tmp-dir",
    type=str,
    default=None,
    help="The directory of temporary files.",
    dest="tmp_dir",
)
args = parser.parse_args()

args.input = os.path.abspath(args.input)
assert os.path.isfile(args.input), f"--input: file ({args.input}) does not exist."
args.reference = os.path.abspath(args.reference)
assert os.path.isfile(
    args.reference
), f"--reference: file ({args.input}) does not exist."
args.gtf_prefix = os.path.abspath(args.gtf_prefix)
GENE_GTF = args.gtf_prefix + ".gene.gtf"
EXON_GTF = args.gtf_prefix + ".exon.gtf"
assert os.path.isfile(
    GENE_GTF
), f"--gtf-prefix: gene GTF file ({GENE_GTF}) does not exist."
assert os.path.isfile(
    EXON_GTF
), f"--gtf-prefix: exon GTF file ({EXON_GTF}) does not exist."
GENE_GTF = BedTool(GENE_GTF)
EXON_GTF = BedTool(EXON_GTF)
args.output_dir = os.path.abspath(args.output_dir)
if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)
SBS_FILE = os.path.join(CURRENT_DIR, f"model/COSMIC_v3.2_SBS_{args.genome}.txt")
assert os.path.isfile(
    SBS_FILE
), f"--genome: COSMIC SBS signatures file ({SBS_FILE}) does not exist."
assert args.threads > 0, "--threads should be greater than 0."
assert args.deep >= 0, "--deep should not less than 0."
if args.vaf:
    if len(args.vaf) == 1:
        args.vaf = [[args.vaf[0], 1]]
    else:
        args.vaf = [list(i) for i in zip(args.vaf[0::2], args.vaf[1::2])]
assert args.upstream >= 0, "--upstream should not be less than 0."
assert args.downstream >= 0, "--downstream should not be less than 0."
if isinstance(args.weight, list):
    if len(args.weight) <= 1:
        args.weight = None
    elif len(args.weight) == 2:
        args.weight = {
            "exon": args.weight[0],
            "intron": args.weight[0],
            "upstream": args.weight[1],
            "downstream": args.weight[1],
        }
    elif len(args.weight) == 3:
        args.weight = {
            "exon": args.weight[0],
            "intron": args.weight[1],
            "upstream": args.weight[2],
            "downstream": args.weight[2],
        }
    else:
        args.weight = {
            "exon": args.weight[0],
            "intron": args.weight[1],
            "upstream": args.weight[2],
            "downstream": args.weight[3],
        }
if args.weight is None:
    args.weight = {"exon": 1, "intron": 1, "upstream": 1, "downstream": 1}
assert args.max_iter > 0, "--max-iter should be greater than 0."
args.log = os.path.abspath(args.log)
if not os.path.isdir(log_dir := os.path.dirname(args.log)):
    os.makedirs(log_dir)
logging.config.dictConfig(
    {
        "version": 1,
        "disable_existing_loggers": True,
        "formatters": {
            "verbose": {
                "format": "%(asctime)s [%(levelname)s] %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
        "handlers": {
            "console": {
                "level": "DEBUG",
                "class": "logging.StreamHandler",
                "formatter": "verbose",
            },
            "file": {
                "level": "DEBUG",
                "class": "logging.FileHandler",
                "delay": True,
                "filename": args.log,
                "formatter": "verbose",
            },
        },
        "loggers": {
            "": {
                "handlers": ["file"] + ([] if args.silent else ["console"]),
                "level": "DEBUG",
            },
        },
    }
)
logger = logging.getLogger()
make_tmp_dir = False
if args.tmp_dir:
    args.tmp_dir = os.path.abspath(args.tmp_dir)
    if not os.path.isdir(args.tmp_dir):
        os.makedirs(args.tmp_dir)
        make_tmp_dir = True
else:
    args.tmp_dir = mkdtemp(dir=CURRENT_DIR)
    make_tmp_dir = True
BASE_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
DF_PATTERN = pd.DataFrame(
    [
        [f"{a}[{c}>{d}]{b}"]
        for a in ["A", "G", "C", "T"]
        for b in ["A", "G", "C", "T"]
        for c in ["C", "T"]
        for d in ["A", "G", "C", "T"]
        if c != d
    ],
    columns=["MUT"],
)


class Progress:
    def __init__(self, queue, total, interval=0.01):
        self.queue = queue
        self.total = total
        self.interval = interval
        self.prog = 0
        self.next_prog = ceil(self.total * interval)
        logger.info("Progress: 0.0% (0 / {})".format(self.total))
        self.start_time = time()
        self.process = Process(target=self.update_progress)
        self.process.start()

    def update(self, n=1):
        self.queue.put(n)

    def update_progress(self):
        while True:
            if self.queue.empty():
                sleep(1)
            else:
                if (update := self.queue.get()) > 0:
                    self.prog += update
                    if self.prog >= self.next_prog:
                        elapsed_time = time() - self.start_time
                        left_time = (
                            elapsed_time * (self.total - self.prog) / self.prog
                        )
                        logger.info(
                            "Progress: {:.1%} ({} / {}) [{} < {}]".format(
                                self.prog / self.total,
                                self.prog,
                                self.total,
                                self.sec2time(elapsed_time),
                                self.sec2time(left_time),
                            )
                        )
                        self.next_prog = ceil(
                            self.total
                            * (int(self.prog / self.total / self.interval) + 1)
                            * self.interval
                        )
                else:
                    break

    @staticmethod
    def sec2time(seconds):
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        return "{:0>2d}:{:0>2d}:{:0>2d}".format(int(h), int(m), int(s))


def open_vcf(vcf_path):
    gz = True
    if vcf_path.lower().endswith(".gz"):
        try:
            with GzipFile(vcf_path, "r") as f:
                header = False
                for line in f:
                    if not (header or line.startswith(b"#CHROM")):
                        continue
                    header = True
                    yield line.decode("utf-8").strip()
        except:
            gz = False
    else:
        gz = False
    if not gz:
        with open(vcf_path, "r") as f:
            header = False
            for line in f:
                if not (header or line.startswith("#CHROM")):
                    continue
                header = True
                yield line.strip()


def get_sequence(pysam_obj, chromosome, position):
    try:
        chromosome = str(chromosome)
        position = int(position)
        return pysam_obj.fetch(chromosome, position - 2, position + 1).upper()
    except:
        chromosome = (
            chromosome[3:] if chromosome.startswith("chr") else f"chr{chromosome}"
        )
        return pysam_obj.fetch(chromosome, position - 2, position + 1).upper()


def get_mut_pattern(
    pysam_obj, chromosome: str, position: int, ref: str, alt: str
):
    try:
        ref_pattern = get_sequence(pysam_obj, chromosome, position)
        if ref in "AG":
            return "{}[{}>{}]{}".format(
                BASE_COMPLEMENT[ref_pattern[2]],
                BASE_COMPLEMENT[ref],
                BASE_COMPLEMENT[alt],
                BASE_COMPLEMENT[ref_pattern[0]],
            )
        else:
            return "{}[{}>{}]{}".format(ref_pattern[0], ref, alt, ref_pattern[2])
    except:
        return ""


def _filter(row, header):
    format = {fmt: ind for ind, fmt in enumerate(row[header["FORMAT"]].split(":"))}
    filter_rst = []
    for sample in row[header["FORMAT"] + 1 :]:
        sample = sample.split(":")
        if len(sample) != len(format):
            filter_rst.append(False)
            continue
        if "1" not in sample[format["GT"]]:
            filter_rst.append(False)
            continue
        if args.deep is None and args.vaf is None:
            filter_rst.append(True)
            continue
        try:
            deep = int(sample[format["DP"]])
            if args.deep and deep < args.deep:
                filter_rst.append(False)
                continue
            if args.vaf:
                split_rst = sample[format["AD"]].split(",")
                if len(split_rst) == 1:
                    vaf = int(split_rst[0]) / deep
                elif len(split_rst) == 2:
                    vaf = int(split_rst[1]) / deep
                else:
                    filter_rst.append(False)
                    continue
                if not any((freq[0] <= vaf <= freq[1] for freq in args.vaf)):
                    filter_rst.append(False)
                    continue
        except:
            filter_rst.append(False)
            continue
        filter_rst.append(True)
    return filter_rst


def _annotate(
    chromosome, position
):
    tmp_vcf = NamedTemporaryFile(
        mode="w+",
        dir=args.tmp_dir,
        encoding="utf-8",
        newline="\n",
        prefix=f"tmp.{os.getpid()}.",
        suffix=".vcf",
    )
    tmp_vcf.write("##fileformat=VCFv4.1\n")
    tmp_vcf.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n")
    tmp_vcf.write(
        "{}\t{}\t.\t.\t.\t.\t.\t.\t.\n".format(
            "chr" + chromosome
            if args.add_chr and not chromosome.startswith("chr")
            else chromosome,
            position,
        )
    )
    tmp_vcf.seek(0)
    snps = BedTool(tmp_vcf)
    intersect_rst = EXON_GTF.intersect(snps, wa=True)
    hit_genes = {i[1] for i in intersect_rst}
    tmp_vcf.seek(0)
    window_rst = str(
        GENE_GTF.window(snps, l=args.upstream, r=args.downstream, sw=True)
    ).strip()
    tmp_vcf.close()
    if not window_rst:
        return ["none_gene"], ["exon"]
    df_window = pd.read_csv(
        StringIO(window_rst),
        delimiter="\t",
        header=None,
        usecols=[1, 3, 4, 6, 10],
        dtype={1: str, 3: int, 4: int, 10: int},
    )
    df_window.columns = ["gene_id", "start", "end", "strand", "position"]
    df_window["gene_id"] = df_window["gene_id"].astype(str)
    df_window["position"] = df_window["position"].astype(int)
    df_window["start"] = df_window["start"].astype(int)
    df_window["end"] = df_window["end"].astype(int)
    df_window.drop_duplicates(subset=["gene_id"], inplace=True)
    gene_ids, mut_types = [], []
    for _, row in df_window.iterrows():
        if row["gene_id"] in hit_genes:
            mut_type = "exon"
        else:
            if row["position"] < row["start"]:
                if row["strand"] == "+":
                    mut_type = "upstream"
                else:
                    mut_type = "downstream"
            elif row["position"] > row["end"]:
                if row["strand"] == "+":
                    mut_type = "downstream"
                else:
                    mut_type = "upstream"
            else:
                mut_type = "intron"
        gene_ids.append(row["gene_id"])
        mut_types.append(mut_type)
    return gene_ids, mut_types


def _filter_and_annotate(
    row, header, queue=None
):
    try:
        filter_rst = _filter(row=row, header=header)
        if not any(filter_rst):
            return []
        anno_rst = list(_annotate(row[header["#CHROM"]], row[header["POS"]]))
        return anno_rst + filter_rst
    except:
        for line in traceback.format_exc().splitlines():
            logger.error(line)
        return []
    finally:
        if queue is not None:
            queue.put(1)


def filter_and_annotate(vcf_path):
    bases = {"A", "G", "C", "T"}
    pysam_obj = pysam.FastaFile(args.reference)
    header = {}
    total_vars, total_tasks = 0, 0
    results = {"MUT": [], "task": []}
    p = Pool(args.threads)
    queue = Manager().Queue()
    logger.info("Reading VCF file...")
    for line in open_vcf(vcf_path):
        if line.startswith("#CHROM") and not header:
            header = OrderedDict(
                ((col, ind) for ind, col in enumerate(line.split("\t")))
            )
        else:
            total_vars += 1
            if len(row := line.split("\t")) == len(header):
                if not (row[header["REF"]] in bases and row[header["ALT"]] in bases):
                    continue
                mut_pattern = get_mut_pattern(
                    pysam_obj,
                    row[header["#CHROM"]],
                    row[header["POS"]],
                    row[header["REF"]],
                    row[header["ALT"]],
                )
                if not mut_pattern:
                    continue
                task = p.apply_async(_filter_and_annotate, (row, header, queue))
                results["MUT"].append(mut_pattern)
                results["task"].append(task)
                total_tasks += 1
    logger.info(
        "Detect {} variants, {} SNVs are selected for further filtering and annotation.".format(
            total_vars, total_tasks
        )
    )
    progress = Progress(queue, total_tasks)
    p.close()
    p.join()
    queue.put(-1)
    progress.process.join()
    df = pd.DataFrame(
        [
            [results["MUT"][ind]] + rst
            for ind in range(len(results["MUT"]))
            if (rst := results["task"][ind].get())
        ],
        columns=["MUT", "GENE", "TYPE"] + list(header.keys())[header["FORMAT"] + 1 :],
    )
    return df


def _get_gsgp(df_mut, sample_id, queue=None):
    try:
        df_mut = pd.DataFrame(
            [
                row[:2].tolist() + list(z)
                for row in df_mut.values
                for z in zip(*row[2:])
            ],
            columns=df_mut.columns,
        )
        df_mut["weight"] = df_mut["TYPE"].map(args.weight)
        df_mut = pd.merge(
            df_mut,
            df_mut.groupby(["index", "MUT"])["weight"].sum(),
            on=["index", "MUT"],
            suffixes=("", "_sum"),
        )
        df_mut["weight"] /= df_mut["weight_sum"]
        df_mut = (
            df_mut.groupby(["MUT", "GENE"])["weight"]
            .sum()
            .reset_index(drop=False)
            .pivot_table(index="MUT", columns="GENE", values="weight")
        )
        df_mut["sample_level"] = df_mut.sum(axis=1)
        df_mut = pd.merge(
            DF_PATTERN,
            df_mut,
            on="MUT",
            how="left",
        ).fillna(0)
        df_mut.set_index("MUT", inplace=True)
        df_mut = (
            df_mut * 2 / df_mut.sum().sum()
        )
        df_mut = df_mut.T
        df_mut = df_mut[df_mut.sum(axis=1) != 0]
        H = pd.read_csv(SBS_FILE, delimiter="\t", low_memory=False, index_col=0).T
        df_mut = df_mut[list(H.columns)]
        if args.save_x:
            df_mut.to_csv(
                os.path.join(args.output_dir, f"{sample_id}.X.txt.gz"),
                sep="\t",
                index_label="GENE",
                compression="gzip",
                line_terminator="\n",
            )
        model = NMF(
            n_components=len(H),
            random_state=args.random_seed,
            max_iter=args.max_iter,
            init="custom",
            tol=1e-16,
        )
        model.n_components_ = len(H)
        model.components_ = H.values
        W = model.transform(df_mut.values)
        df_signature = pd.DataFrame(
            W,
            index=list(df_mut.index),
            columns=list(H.index),
        )
        df_signature.to_csv(
            signature_file := os.path.join(
                args.output_dir, f"{sample_id}.signatures.txt.gz"
            ),
            sep="\t",
            index_label="GENE",
            compression="gzip",
            line_terminator="\n",
        )
        return signature_file
    except:
        logger.error(f"Sample {sample_id} failed in calculating GSGPs.")
        for line in traceback.format_exc().splitlines():
            logger.error(line)
        return ""
    finally:
        if queue is not None:
            queue.put(1)


def get_gsgp(df_filter):
    p = Pool(args.threads)
    queue = Manager().Queue()
    n_tasks = 0
    for sample in list(df_filter.columns)[3:]:
        df_mut = df_filter[["MUT", "GENE", "TYPE"]][df_filter[sample]].reset_index(
            drop=False
        )
        if len(df_mut) > 0:
            logger.info(f"{len(df_mut)} SNP(s) in sample {sample} passed filter.")
        else:
            logger.warning(
                f"All SNPs in sample {sample} failed to pass filter, skip calculating GSGPs."
            )
            continue
        p.apply_async(_get_gsgp, (df_mut, sample, queue))
        n_tasks += 1
    logger.info(f"{n_tasks} samples were used for calculation of GSGPs.")
    progress = Progress(queue, n_tasks)
    p.close()
    p.join()
    queue.put(-1)
    progress.process.join()
    logger.info("Results were saved in {}.".format(args.output_dir))


if __name__ == "__main__":
    logger.info("Parameters:")
    for k, v in args.__dict__.items():
        logger.info(f"{k}: {v}")
    get_gsgp(filter_and_annotate(args.input))
    if make_tmp_dir:
        rmtree(args.tmp_dir)
    logger.info("Done!")