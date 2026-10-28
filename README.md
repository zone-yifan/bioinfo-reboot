bioinfo-reboot
生物信息学「重新出发」——一条从零跑通的 RNA-seq 最小分析流水线练习仓库。
从模拟测序数据出发，完成 质控 → 修剪 → 比对 → 计数 全流程，输出基因表达计数矩阵。
一句话
用自包含的模拟 reads 跑通 fastqc → trim_galore → minimap2 → samtools → featureCounts → multiqc，
理解真实 RNA-seq 分析每一步的输入、输出与原理，为后续差异表达（DESeq2）打基础。
背景与目标
- 为什么叫 reboot：中断后重新系统学习生物信息学，按「可复现、可追溯、先跑通最小闭环」的思路推进。
- 为什么用模拟数据：真实 RNA-seq 原始文件动辄几百 MB～几 GB，新手阶段先下载容易卡在网络/磁盘。
这里用脚本 sim_reads.py 从参考基因组截取生成练习 reads，零下载、自包含、流程形状与真实数据完全一致。
- 一个彩蛋：reads 本就截自 sarscov2.fa，比对回去映射率接近 100%，正好用来验证管线没有接错。
分析流程
sim_reads.py (生成模拟 reads)
        │
        ▼
fastqc_raw/           —— 原始 reads 质量评估
        │
        ▼
trim_galore (trimmed/)—— 修剪接头与低质量末端
        │
        ▼
fastqc_trim/          —— 修剪后质量复查
        │
        ▼
minimap2 比对 (aln.sam → aln.bam → aln.bam.bai)
        │
        ▼
featureCounts (counts.txt) —— 基因表达计数矩阵
        │
        ▼
multiqc_report/       —— 全部质控报告汇总
目录结构
bioinfo-reboot/
├── sim_reads.py                  # 模拟测序 reads 生成脚本
├── sarscov2.fa / sarscov2.gtf    # 参考基因组（NC_045512.2）与最小注释
├── reads.fq                      # 生成的练习 reads（500×75bp）
├── fastqc_raw/                   # 原始 reads FastQC 报告
├── trimmed/                      # trim_galore 修剪后 reads
├── fastqc_trim/                  # 修剪后 FastQC 报告
├── aln.sam / aln.bam / .bai      # minimap2 比对与排序、索引结果
├── counts.txt / .summary         # featureCounts 计数矩阵
├── multiqc_report/               # MultiQC 汇总报告
├── multiqc_config.yml            # MultiQC 配置
└── 第3-4周生信核心实操清单.md     # 完整操作步骤与报告解读
依赖工具
- FastQC —— 测序质量评估
- Trim Galore —— 接头/低质量修剪
- minimap2 —— 长读/短读比对（本流程用 -ax sr 短读模式）
- samtools —— SAM/BAM 处理与 flagstat 统计
- featureCounts（Subread）—— 基因水平定量
- MultiQC —— 多报告聚合
- seqkit —— 序列统计
建议在 Conda 环境（如 bioinfo）中通过 Bioconda 安装：
mamba install -y -c bioconda -c conda-forge fastqc multiqc trim-galore minimap2 samtools subread seqkit
快速复现
# 1. 生成模拟 reads（500 条、读长 75、错误率 2%）
python3 sim_reads.py sarscov2.fa reads.fq 500 75 0.02

# 2. 原始 reads 质控
fastqc reads.fq --outdir fastqc_raw

# 3. 修剪
trim_galore reads.fq --output_dir trimmed

# 4. 修剪后质控
fastqc trimmed/reads_trimmed.fq --outdir fastqc_trim

# 5. 建索引并比对
minimap2 -d sarscov2.mmi sarscov2.fa
minimap2 -ax sr sarscov2.mmi trimmed/reads_trimmed.fq.gz > aln.sam

# 6. 排序 + 索引 + 统计（映射率应接近 100%）
samtools sort -o aln.bam aln.sam
samtools index aln.bam
samtools flagstat aln.bam

# 7. 基因计数
featureCounts -a sarscov2.gtf -o counts.txt -t exon -g gene_id aln.bam

# 8. 汇总报告
multiqc fastqc_raw fastqc_trim -o multiqc_report
结果示例
counts.txt 输出（featureCounts v2.1.1），得到 5 个 SARS-CoV-2 基因的 reads 计数：
Geneid
长度 (bp)
计数 (reads)
ORF1ab
21290
345
S
3822
71
N
1260
24
M
669
7
E
228
5
计数高低与基因长度大致正相关（ORF1ab 最长、覆盖最多），符合随机均匀取样的预期，验证了比对与定量步骤正确。
质控结论（Week 3 核心）
- 模拟数据质量分恒定在 Q34–35，Per Base Sequence Quality 全程绿色 PASS；
- trim_galore 成功把部分 reads 从 75 bp 剪短到 71–75 bp，说明修剪生效；
- GC 约 37%，与 SARS-CoV-2 基因组 GC（37.97%）吻合，验证模拟是按真实参考截取；
- 个别琥珀告警（长度分布、过表达序列）均源于「小样本模拟随机数据」，非真实质控失败。
进度与下一步
[x] Week 1–2：WSL2 环境搭建、Git/GitHub、seqkit 实战
[x] Week 3：质控与修剪闭环（fastqc → trim → fastqc → multiqc）
[x] Week 4：比对与计数（minimap2 → samtools → featureCounts）
[ ] Week 5+：差异表达分析（DESeq2，R 语言），需引入真实数据与生物学重复
关于作者
生物分子 + 生物信息学背景，正在系统重建生物信息学技能栈。
本仓库是「最小闭环」复现练习，详细步骤与报告解读见
第3-4周生信核心实操清单.md。
