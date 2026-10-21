#!/usr/bin/env python3
# sim_reads.py — 从参考基因组随机截取 reads，生成「练习用」FASTQ
# 用途：在你还没下载真实 RNA-seq 数据前，先用它熟悉 fastqc / trim / 比对 全流程。
# 注意：这是「模拟」数据（reads 来自你已有的 sarscov2.fa），仅用于练手，
#       不是真实测序结果；但流程形状与真实数据完全一致。
#
# 用法：
#   python3 sim_reads.py <参考.fa> <输出.fq> [reads数] [读长] [错误率]
# 例：
#   python3 sim_reads.py sarscov2.fa reads.fq 500 75 0.02

import sys
import random

BASES = "ACGT"


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def qual_encode(phreds):
    # Phred+33 编码（Sanger / Illumina 1.8+）
    out = []
    for q in phreds:
        q = max(0, min(93, q))
        out.append(chr(q + 33))
    return "".join(out)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    ref = sys.argv[1]
    out = sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    L = int(sys.argv[4]) if len(sys.argv) > 4 else 75
    err = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02

    seq = read_fasta(ref)
    if len(seq) < L:
        print("参考序列太短，无法截取读长 %d" % L)
        sys.exit(1)

    random.seed(42)
    with open(out, "w") as o:
        for i in range(n):
            start = random.randint(0, len(seq) - L)
            sub = seq[start:start + L]
            # 以 err 概率引入错误碱基，模拟测序错误
            err_sub = "".join(
                random.choice(BASES) if random.random() < err else b
                for b in sub
            )
            # 质量值：30~40 之间，模拟较高质量测序
            phred = [random.randint(30, 40) for _ in range(L)]
            o.write("@%s_read%d length=%d\n" % (out, i, L))
            o.write(err_sub + "\n")
            o.write("+\n")
            o.write(qual_encode(phred) + "\n")

    print("已生成 %d 条 reads（读长 %d，错误率 %.2f）-> %s" % (n, L, err, out))


if __name__ == "__main__":
    main()
