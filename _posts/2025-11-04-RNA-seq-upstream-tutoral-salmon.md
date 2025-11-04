---
title: RNA-seq analysis tutorial (salmon)
date: 2025-11-4 11:45
categories: [NGS Tutorial,RNA-seq Tutorial]
tags: [rna-seq,rna]
---

RNA-seq（RNA测序）是研究转录组的重要高通量技术，其标准分析流程通常包括：RNA提取、文库构建、高通量测序、数据质控、序列比对（或准比对）、基因/转录本表达定量、差异表达分析及功能富集分析等步骤。传统分析依赖于“比对→计数”的策略，即先将测序reads比对到参考基因组或转录组，再通过计数工具（如featureCounts）统计每个基因的reads数，进而进行表达量计算与下游分析。这一流程虽然成熟，但比对步骤耗时较长，且对计算资源要求较高。

相比之下，Salmon 作为一款基于“准映射”（quasi-mapping）的非比对型转录本定量工具，在基因表达定量环节展现出显著优势。它不依赖传统的序列比对，而是通过构建转录组索引，直接将测序reads快速映射到转录本集合，结合概率模型进行转录本丰度的精确估计。这一策略使Salmon具备“闪电般的速度”与低资源消耗，在处理大规模数据时显著提升分析效率，尤其适合多样本、高通量的研究场景。同时，Salmon采用选择性对齐和并行化变分推断技术，在保证超高准确性的同时，有效避免了错误归因，即使对低表达转录本也能提供可靠的定量结果。此外，其支持“诱饵”序列添加和多种索引构建策略，进一步增强了分析的灵活性与鲁棒性。因此，Salmon不仅简化了RNA-seq分析流程，更以快速、准确、低门槛的特点，成为现代转录组研究中不可或缺的利器。


# 1. 数据准备与质控（Quality Control）
- 获取原始测序数据（FASTQ 格式），通常为双端测序（_R1 和 _R2）。
- 去接头（Adapter Trimming）：测序过程中可能因插入片段较短导致读段读到接头序列，若不剔除会影响后续比对或定量准确性。
- 质量修剪：去除低质量碱基（如 Phred 质量值 < 20）可减少错误匹配，提高定量可靠性。
- 长度过滤：保留足够长度的 reads（如 ≥50 bp），确保其具备足够的比对信息。

命令示例（双端测序）：
```
fastp -l 25 -w 4 \
--detect_adapter_for_pe \
-i input_R1.fastq \
-I input_R2.fastq \
-o output_R1.fastq \
-O output_R2.fastq
```
参数解释：
- -l 表示保留最终序列长度大于25的reads
-  -w 表示线程
-  --detect\_adapter\_for_pe 双端测序参数
-  -i, -I 分别表示read1，read2两个输入文件
-  -o, -O 分别表示read1，read2对应的输出文件

# 2. 构建转录组索引（Build Transcriptome Index）
- 下载参考基因组对应的转录本序列（FASTA 格式，如从 Ensembl 或 NCBI 获取）。
- 索引构建是定量的前提：Salmon 不直接比对 reads 到基因组，而是基于转录本集合建立哈希表索引，实现快速查找。
- “准映射”机制：通过 k-mer 匹配将 reads 映射到可能的转录本，无需全局比对，极大提升速度。
- 轻量高效：相比 HISAT2/STAR 建索引更快速，占用内存更少。

命令示例：
```
salmon index -p 4 -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index -k 31 --type quasi
```
参数解释：
- -t：输入转录本序列文件（包含所有已知转录本的cDNA序列）。
- -i：输出索引目录名。
- --type quasi：使用“准映射”模式，是 Salmon 的核心加速机制。
- -k：k-mer 长度，默认 31，适用于大多数测序数据。

<details markdown="1">
<summary>为什么你的定量结果中没有基因ID?</summary>
如果定量时不提供转录本ID-基因ID对照表，那么你的定量结果只有转录本水平的定量结果。  
转录谱文件可以到[Ensembl](https://useast.ensembl.org/info/data/ftp/index.html)下载, 转录谱参考文件内容如下：
  
```
>ENST00000632684.1 cdna chromosome:GRCh38:7:142786213:142786224:1 gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T cell receptor beta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12158]
GGGACAGGGGGC
```

第一行是header信息，以`>`开头，第一列 `ENST00000632684.1` 为转录本ID，salmon默认对转录本进行定量，但是一个基因有多个转录本，因此需要把同一基因的各个转录本定量结果合并在一起，所以需要提取一个转录本-基因对应的信息文件。`gene_symbol:TRBD1`中的`TRBD1 `即为基因名字，把第一列和基因列提取出来就行了。

提取命令如下：
```
zgrep ">" Homo_sapiens.GRCh38.cdna.all.fa.gz | sed 's/>//g' | sed 's/cdna.*gene_symbol://g' | sed 's/description.*//g' > gene_map.txt
```
</details>

# 3. 使用 Salmon 进行基因/转录本表达定量
- 将清洗后的 clean reads 映射到转录组索引，并估计每个转录本的表达丰度。
- 非比对式定量：Salmon 不生成 SAM/BAM 文件，而是直接输出每个转录本的表达量（TPM、NumReads）。
- 概率分配机制：对于多映射 reads（比对到多个转录本），Salmon 使用 EM 算法进行概率性分配，避免误计。
- 输出标准化结果：直接提供 TPM（Transcripts Per Million），适合跨样本表达量比较。
- 支持可变剪接建模：能区分来自同一基因的不同转录本，实现转录本水平定量。
  
命令示例（双端测序）：
```
salmon quant -p 4 \
--validateMappings \
-i salmon_index \
-l A -g gene_map.txt \
-1 output_R1.fastq \
-2 output_R2.fastq \
-o Salmon_out/output_name
```
参数解释：
- -p	表示线程
-  --validateMappings	启用选择性比对到转录本。能够提高比对的灵敏度和特异性从而提高了定量的准确度。
-  -l	library 类型，默认A，自动检测
-  -g	转录本-基因之间的对应关系，如果不指定，不会输出对基因表达的定量， 而是输出对转录本表达的定量
-  -1,-2 质控后的read1，read2两个fastq文件
-  -o 输出文件夹，Salmon_out是所有样本输出的目录，output_name是本次running的样本定量结果输出目录

输出目录结构如下：

```
Salmon_out
├── CC_Cas9-1
├── CC_Cas9-2
├── CC_Cas9-3
├── CC_Cas9_Cre-1
├── CC_Cas9_Cre-2
├── CC_Cas9_Cre-3
├── CC_Cas9_Flpo-1
├── CC_Cas9_Flpo-2
└── CC_Cas9_Flpo-3
```
每个目录下文件结构如下：

```
Salmon_out/CC_Cas9-1
├── aux_info
├── cmd_info.json
├── libParams
├── lib_format_counts.json
├── logs
├── quant.genes.sf
└── quant.sf
```
其中quant.sf为转录本定量结果：

```
Name    Length  EffectiveLength TPM     NumReads
ENST00000474626.5       639     522.973 75.863174       32.032
ENST00000476222.5       685     690.231 8.012039        4.465
ENST00000496813.5       642     654.652 2.670470        1.411
ENST00000490191.5       996     844.000 0.000000        0.000
ENST00000472218.1       1109    957.000 0.000000        0.000
ENST00000479802.1       413     261.136 281.314078      59.311
```
如果指定了参数-g，则会输出基因定量结果，存储在文件quant.genes.sf中：

```
Name    Length  EffectiveLength TPM     NumReads
AC008626.1      255     775.152 202.283 126.596
VN1R107P        209     80      0       0
Z97832.1        328     179     0       0
OR7E13P 910     554.755 2.23267 1
AC006511.1      315     167     0       0
MTCO3P12        547     355.804 8.42352 2.42
MTCO2P12        682     531     0       0
MTND2P28        1044    730.596 29.969  17.678
```
# 4. 合并样本定量结果，构建表达矩阵
#### 用shell脚本提取基因名称和基因定量两列数据，以table分割
```
#!/bin/bash
# 脚本名称：extract_gene_cols.sh

for i in `ls $1`
do
    sample_name=$i
    awk '{print$1,$5}' $1/$i//quant.genes.sf | grep -v 'Name' | sed 's/ /\t/' > $1/$i/${i}.count
done
```
使用方法: `bash extract_gene_cols.sh Salmon_out`, 参数Salmon\_out即为所有样本的输出文件夹。

#### 用python脚本合并创建表达矩阵
```
#!/usr/local/bin/python
# -*- coding:utf-8 -*-
# 脚本名：generate_count_matrix.py
import sys
import pandas as pd
from glob import iglob

sample_path = iglob(sys.argv[1] + "/*/*.count")
dfs = []

for df in sample_path:
    fn = df.split('/')[-1]
    sample_name = fn[:-6]
    dfs.append(pd.read_table(df, header=None, index_col=0,sep='\t', names=["geneid", sample_name]))

matrix = pd.concat(dfs, axis=1)
matrix.to_csv(sys.argv[2])
```
使用方法: `python generate_count_matrix.py Salmon_out Count_matrix.csv`, 第一个参数Salmon\_out即为所有样本的输出文件夹, 第二个参数Count\_matrix.csv为输出表达矩阵文件名字。

最终生成的Count_matrix.csv表达矩阵为行为基因，列为样本的表达矩阵：

```
geneid,CC_Cas9_Flpo-1,CC_Cas9-1,CC_Cas9_Flpo-2,CC_Cas9_Cre-1,CC_Cas9-2,CC_Cas9_Cre-3,CC_Cas9_Flpo-3,CC_Cas9-3,CC_Cas9_Cre-2
AC217785.3,2.639,0.0,1.9909999999999999,0.0,0.0,0.0,0.0,0.0,0.0
FRG1KP,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.915,0.0
AC007322.4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
RPS6P20,511.55699999999996,975.587,1573.53,1309.35,689.59,1340.45,308.182,752.914,217.21400000000003
OR13I1P,0.986,6.224,20.303,5.135,4.225,10.315,0.0,3.5069999999999997,2.539
```
> # 一个完整的snakemake流程脚本

```
from glob import iglob
import pandas as pd

r1 = iglob('SeqData/RNA/*/*_R1.fq.gz')

sampleInfo = pd.DataFrame()
sampleInfo['r1'] = [i for i in r1]
sampleInfo['r1'] = [i for i in r1]
sampleInfo['name'] = sampleInfo.r1.str.extract('SeqData/RNA/.*/(.*)_R1.fq.gz', expand=True)

rule all:
        input:
                "Ref/genemap.txt",
                "Ref/salmon_index",
                "qc_metric.txt",
                "RunData/RNA/Count_matrix.csv",
                "RunData/RNA/TPM_matrix.csv",
                "RunData/RNA/FPKM_matrix.csv"

rule buildGM:
        input:
                "Ref/GCF_000001635.27_GRCm39_rna.fna.gz"
        output:
                "Ref/genemap.txt"
        shell:
                '''
                        zgrep ">" {input} | sed 's/>//g' | sed 's/ .*//g' > Ref/transID.txt
                        zgrep ">" {input} | sed 's/.*(//g' | sed 's/).*//g' > Ref/geneID.txt
                        paste -d '\t' Ref/transID.txt Ref/geneID.txt > {output}
                        rm Ref/transID.txt Ref/geneID.txt
                '''
rule buildIndex:
        input:
                "Ref/GCF_000001635.27_GRCm39_rna.fna.gz"
        output:
                "Ref/salmon_index"
        threads: 4
        shell:
                '''
                        docker run --rm \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/salmon salmon index -p {threads} \
                        -t $PWD/Ref/GCF_000001635.27_GRCm39_rna.fna.gz \
                        -i $PWD/Ref/salmon_index -k 31
                '''
rule fastp:
        input:
                r1 = "SeqData/RNA/{name}/{name}_R1.fq.gz",
                r2 = "SeqData/RNA/{name}/{name}_R2.fq.gz"
        output:
                r1 = "RunData/RNA/Fastp/{name}_R1.fq.gz",
                r2 = "RunData/RNA/Fastp/{name}_R2.fq.gz"
        threads: 4
        shell:
                '''
                        docker run --rm --user root \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/fastp fastp -w {threads} \
                        --detect_adapter_for_pe \
                        --compression=2 \
                        -i $PWD/SeqData/RNA/{wildcards.name}/{wildcards.name}_R1.fq.gz \
                        -I $PWD/SeqData/RNA/{wildcards.name}/{wildcards.name}_R2.fq.gz \
                        -o $PWD/RunData/RNA/Fastp/{wildcards.name}_R1.fq.gz \
                        -O $PWD/RunData/RNA/Fastp/{wildcards.name}_R2.fq.gz \
                        -h $PWD/RunData/RNA/Fastp/{wildcards.name}.html \
                        -j $PWD/RunData/RNA/Fastp/{wildcards.name}.json
                        -h $PWD/RunData/RNA/Fastp/{wildcards.name}.html \
                        -j $PWD/RunData/RNA/Fastp/{wildcards.name}.json
                '''
rule salmon:
        input:
                r1 = "RunData/RNA/Fastp/{name}_R1.fq.gz",
                r2 = "RunData/RNA/Fastp/{name}_R2.fq.gz",
                genemap = "Ref/genemap.txt",
                index = "Ref/salmon_index"
        output:
                "RunData/RNA/Salmon/{name}/quant.genes.sf"
        threads: 4
        shell:
                '''
                        docker run --rm --user root \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/salmon salmon quant -p {threads} \
                        --validateMappings \
                        -i $PWD/Ref/salmon_index \
                        -l A -g $PWD/Ref/genemap.txt \
                        -1 $PWD/RunData/RNA/Fastp/{wildcards.name}_R1.fq.gz \
                        -2 $PWD/RunData/RNA/Fastp/{wildcards.name}_R2.fq.gz \
root@insight:/mnt/Data2/home/Goushixue/project-ISSAAC-seq-embryo-test-data# vim Snakefile.rna
from glob import iglob
import pandas as pd

r1 = iglob('SeqData/RNA/*/*_R1.fq.gz')

sampleInfo = pd.DataFrame()
sampleInfo['r1'] = [i for i in r1]
sampleInfo['name'] = sampleInfo.r1.str.extract('SeqData/RNA/.*/(.*)_R1.fq.gz', expand=True)

rule all:
        input:
                "Ref/genemap.txt",
                "Ref/salmon_index",
                "qc_metric.txt",
                "RunData/RNA/Count_matrix.csv",
                "RunData/RNA/TPM_matrix.csv",
                "RunData/RNA/FPKM_matrix.csv"

rule buildGM:
        input:
                "Ref/GCF_000001635.27_GRCm39_rna.fna.gz"
        output:
        output:
                "Ref/genemap.txt"
        shell:
                '''
                        zgrep ">" {input} | sed 's/>//g' | sed 's/ .*//g' > Ref/transID.txt
                        zgrep ">" {input} | sed 's/.*(//g' | sed 's/).*//g' > Ref/geneID.txt
                        paste -d '\t' Ref/transID.txt Ref/geneID.txt > {output}
                        rm Ref/transID.txt Ref/geneID.txt
                '''
rule buildIndex:
        input:
                "Ref/GCF_000001635.27_GRCm39_rna.fna.gz"
        output:
                "Ref/salmon_index"
        threads: 4
        shell:
                '''
                        docker run --rm \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/salmon salmon index -p {threads} \
                        -t $PWD/Ref/GCF_000001635.27_GRCm39_rna.fna.gz \
                        -i $PWD/Ref/salmon_index -k 31
                '''
rule fastp:
        input:
                r1 = "SeqData/RNA/{name}/{name}_R1.fq.gz",
                r2 = "SeqData/RNA/{name}/{name}_R2.fq.gz"
        output:
                r1 = "RunData/RNA/Fastp/{name}_R1.fq.gz",
                r2 = "RunData/RNA/Fastp/{name}_R2.fq.gz"
        threads: 4
        shell:
                '''
                        docker run --rm --user root \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/fastp fastp -w {threads} \
                        --detect_adapter_for_pe \
                        --compression=2 \
                        -i $PWD/SeqData/RNA/{wildcards.name}/{wildcards.name}_R1.fq.gz \
                        -I $PWD/SeqData/RNA/{wildcards.name}/{wildcards.name}_R2.fq.gz \
                        -o $PWD/RunData/RNA/Fastp/{wildcards.name}_R1.fq.gz \
                        -O $PWD/RunData/RNA/Fastp/{wildcards.name}_R2.fq.gz \
                        -h $PWD/RunData/RNA/Fastp/{wildcards.name}.html \
                        -j $PWD/RunData/RNA/Fastp/{wildcards.name}.json
                '''
rule salmon:
        input:
                r1 = "RunData/RNA/Fastp/{name}_R1.fq.gz",
                r2 = "RunData/RNA/Fastp/{name}_R2.fq.gz",
                genemap = "Ref/genemap.txt",
                index = "Ref/salmon_index"
        output:
                "RunData/RNA/Salmon/{name}/quant.genes.sf"
        threads: 4
        shell:
                '''
                        docker run --rm --user root \
                        -v /mnt/Data2:/mnt/Data2 \
                        biodockers/salmon salmon quant -p {threads} \
                        --validateMappings \
                        -i $PWD/Ref/salmon_index \
                        -l A -g $PWD/Ref/genemap.txt \
                        -1 $PWD/RunData/RNA/Fastp/{wildcards.name}_R1.fq.gz \
                        -2 $PWD/RunData/RNA/Fastp/{wildcards.name}_R2.fq.gz \
                        -o $PWD/RunData/RNA/Salmon/{wildcards.name}
                '''
rule qcmetric:
        input:
                expand("RunData/RNA/Salmon/{name}/quant.genes.sf", name = sampleInfo['name'])
        output:
                "qc_metric.txt"
        shell:
                '''
                        echo -e "Name\tMappingRate\tMappedRead\tTotalReads" > qc_metric.txt

                        for i in RunData/RNA/Salmon/*/logs/salmon_quant.log; \
                        do name=$(echo $i | sed 's/\/logs.*//g' | sed 's/.*Salmon\///g'); \
                        mprate=$(grep "Mapping rate" $i | sed 's/.*= //g' | sed 's/\%//g'); \
                        mpreads=$(grep "Counted" $i | sed 's/.*Counted //g' | sed 's/ total.*//g'); \
                        ttreads=$(grep "Observed" $i | sed 's/Observed //g' | sed 's/ total.*//g'); \
                        echo -e "${{name}}\t${{mprate}}\t${{mpreads}}\t${{ttreads}}" >> qc_metric.txt; done
                '''
rule getMatrix:
        input:
                expand("RunData/RNA/Salmon/{name}/quant.genes.sf", name = sampleInfo['name'])
        output:
                "RunData/RNA/TPM_matrix.csv",
                "RunData/RNA/FPKM_matrix.csv",
                "RunData/RNA/Count_matrix.csv"
        run:
                sample_path = iglob("RunData/RNA/Salmon/*/quant.genes.sf")
                dfs_mtx = dfs_tpm = dfs_fpkm = []
                for sp in sample_path:
                        sample_name = sp.split('/')[-2]
                        temp_df = pd.read_table(sp, header = 0, sep="\t", index_col = 0)
                        temp_df['FPKM'] = (temp_df.NumReads / temp_df.NumReads.sum()) / temp_df.EffectiveLength * 10e9
                        dfs_mtx.append(pd.DataFrame(temp_df.loc[:,"NumReads"]).rename(columns={"NumReads": sample_name}))
                        dfs_tpm.append(pd.DataFrame(temp_df.loc[:,"TPM"]).rename(columns={"TPM": sample_name}))
                        dfs_fpkm.append(pd.DataFrame(temp_df.loc[:,"FPKM"]).rename(columns={"FPKM": sample_name}))
                pd.concat(dfs_mtx, axis=1).to_csv("RunData/RNA/Count_matrix.csv")
                pd.concat(dfs_tpm, axis=1).to_csv("RunData/RNA/TPM_matrix.csv")
                pd.concat(dfs_fpkm, axis=1).to_csv("RunData/RNA/FPKM_matrix.csv")
```
