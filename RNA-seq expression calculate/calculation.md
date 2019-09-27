

# RNA-seq gene/isoform expression calculation 的流程

整体流程：fastq + STAR → （gene expression）
                      → （bam） + eXpress → （isoform expression）

**安装：**

```linux
1.SRAtoolkit安装地址：https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.4/sratoolkit.2.9.4-ubuntu64.tar.gz
2.STAR安装：https://github.com/alexdobin/STAR
3.eXpress安装:https://pachterlab.github.io/eXpress/overview.html 
```


## 1. fastq文件
以Deng数据集为例，可以在 获得Deng数据集所有测序数据的下载地址，然后通过SRAtoolkit中的fastq-dump工具下载得到fastq文件

```linux
#若为双端测序数据
fastq-dump --gzip --split-files -O fastq  + 数据地址
#若为单端测序数据
fastq-dump --gzip -O fastq + 数据地址
```

## 2. STAR+fastq
利用STAR和fastq可以生成gene expression文件和bam文件，此步骤需要用到Transcriptom annotation (gtf)文件，gtf文件的获得可以联系yingy.l@qq.com

### 2.1 建立人/鼠索引

以鼠为例

```linux
#--runThreadN 为线程数，--genomeDir 为生成的索引文件的地址, --genomeFastaFiles 为fa文件存放的地址，--sjdbGTFfile 为gtf文件存放的地址
STAR \
--runThreadN 16 \   
--runMode genomeGenerate \
--genomeDir ./index_mus \
--genomeFastaFiles ./genome_mus/Mus_musculus.GRCm38.75.dna.SORTED.fa \
--sjdbGTFfile ./genome_mus/Mus_musculus.GRCm38.75.gtf \ 
--sjdbOverhang 50
```

### 2.2 比对

双端测序的就会有read1.fastq和read2.fastq

单端测序的就是read.fastq

```linux
#prefix 为STAR生成文件的前缀名，--sjdbFileChrStartEnd 和 --genomeDir 为相关索引文件的地址
STAR --runThreadN 16 \
            --quantMode TranscriptomeSAM GeneCounts \
            --outFileNamePrefix prefix \
            --sjdbFileChrStartEnd ./index_mus/sjdbList.out.tab \
            --outSAMtype BAM \
            SortedByCoordinate \
            --genomeDir ./index_mus \
            --readFilesCommand gunzip -c --readFilesIn read1.fastq read2.fastq

```

## 3. eXpress+bam
利用eXpress和bam可以生成isoform expression文件，此步骤需要用到人/鼠的Isoform Sequence文件，文件的获得可以联系yingy.l@qq.com

```linux
#./Deng_star_output 为文件输出地址，./genome_mus/mouse_tcx.fa 为Isoform Sequence文件地址，tcxbam 为STAR生成的bam文件

express -o ./Deng_star_output ./genome_mus/mouse_tcx.fa tcxbam
```

说明文档：

[https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

[https://pachterlab.github.io/eXpress/overview.html](https://pachterlab.github.io/eXpress/overview.html)

