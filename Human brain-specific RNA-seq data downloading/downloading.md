

#  从NCBI下载人脑组织特异性数据的流程

整体流程：下载可能的人脑组织的RNA-seq元数据 → 按一定规则过滤掉一些experiment和dataset → 下载每个experiment所用sample的网页 → 提取组织信息并过滤 → 下载RNA-seq data



## 1.从NCBI下载可能的人脑组织的RNA-seq元数据

    (1)In NCBI, select SRA, use key words “(rna-seq illumina) AND "Homo sapiens"[orgn:__txid9606] AND (brain) ” to search
    
    (2)save search results to SraRunInfo
    


## 2. 过滤SraRunInfo
  
    (1)过滤掉reads（spots）的数目少于2,000,000的experiment
    
    (2)过滤掉不同experiment的数目少于10的dataset(SRAStudy)
    
代码：SraRunInfo_filter.py



## 3. 下载SraRunInfo中剩下的SRA对应的网页

代码：download_tissue_information_exp.pl



## 4. 从3中下载的网页中提取出sampleID，并下载sample对应的网页

代码：Extract_SampleID_and_download.py



## 5. 从4中下载的网页提取组织信息

代码：extract_tissue_information.py



## 6. 人工筛选tissue为brain相关的SRA并下载
