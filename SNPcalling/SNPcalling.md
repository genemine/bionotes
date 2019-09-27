

# SNP calling 的流程

整体流程：BWA-MEM + Samtools + bcftools

为什么选择该流程的参考综述  [http://www.nature.com/articles/srep17875](http://www.nature.com/articles/srep17875)

**安装：**

```linux
sudo apt install bwa
sudo apt install samtools
# 也可以源码安装
```

## 1. BWA-MEM

使用BWA整个比对过程主要分为两步，第一步建索引，第二步使用BWA MEM进行比对

### 1.1 建立索引

```linux
bwa index -a bwtsw hg38.fa
```

对参考序列建立索引，参考序列下载地址[http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)

### 1.2 比对

双端测序的就会有read1.fastq和read2.fastq

单端测序的就是read.fastq

```linux
bwa mem -t 20 ref.fa read1.fastq read2.fastq > aln-pe.sam
```

## 2. Samtools + bcftools SNP calling



```linux
# 2.1 将sam转换成bam格式
samtools view -bS aln-pe.sam > aln-pe.bam
# 2.2 排序 （-@是线程数）
samtools sort -@ 20 aln-pe.bam  > aln-pe.sort.bam
# 2.3 索引
samtools index aln-pe.sort.bam
# 2.4 使用bcftools工具call SNP
bcftools mpileup -Ou -f ref.fa aln-pe.sort.bam | bcftools call -Ou -mv | bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > var.flt.vcf
#如果多样本就直接将所有的bam文件 都放到 aln-pe.sort.bam 的位置
```



说明文档：

 [http://samtools.github.io/bcftools/bcftools.html](http://samtools.github.io/bcftools/bcftools.html)

[https://samtools.github.io/bcftools/howtos/variant-calling.html](https://samtools.github.io/bcftools/howtos/variant-calling.html)

