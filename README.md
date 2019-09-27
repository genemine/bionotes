
# 1. SNP calling pipeline

使用 *bwa* + *samtools* + *bcftools* 的组合做 SNP calling 的工作。

输入数据是 fastq 类型的测序数据，输出是 SNP 文件。

Directory is [./SNPcalling](https://github.com/genemine/bionotes/tree/master/SNPcalling)

Documents is [SNPcalling.md](https://github.com/genemine/bionotes/blob/master/SNPcalling/SNPcalling.md)



# 2. RNA-seq based gene/isoform expression calculation

目的:利用STAR、eXpress计算单细胞gene和isoform的表达量

整体流程：fastq + STAR → （gene expression）
                      → （bam） + eXpress → （isoform expression）

Directory is [./RNA-seq expression calculate](https://github.com/genemine/bionotes/tree/master/RNA-seq%20expression%20calculate)

Documents is [calculation.md](https://github.com/genemine/bionotes/blob/master/RNA-seq%20expression%20calculate/calculation.md)


# 3. Single cell RNA-seq Clustering
包含单细胞聚类相关算法工具（R语言）

目的: 集合基于单细胞RNAseq数据的多种聚类算法

整体流程：gene expression matrix (row/genes col/cells) → 多种算法 → 聚类结果评价 (Evaluation)

Directory is [./Single-cell-RNA-seq-Clustering](https://github.com/genemine/bionotes/tree/master/Single-cell-RNA-seq-Clustering)

Documents is [Readme.md](https://github.com/genemine/bionotes/blob/master/Single-cell-RNA-seq-Clustering/README.md)


# 4. Download human brain-specific RNA-seq data

目的：从NCBI下载人脑组织RNA-seq数据

整体流程：下载可能的人脑组织的RNA-seq元数据 → 按一定规则过滤掉一些experiment和dataset → 下载每个experiment所用sample的网页 → 提取组织信息并过滤 → 下载RNA-seq data
Directory is 
Documents is 
