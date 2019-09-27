import re
import os
import pandas as pd

#extract //tissue  //<th>sample composition</th><td>      //<th>brain region</th><td>  //<tr><th>source name</th><td>
# //<th>cell type</th><td>  //<dt>Description</dt><dd class="wide"><p class="first">  //<th>organoid culture</th><td>   //<th>starting cell</th><td>
def extract_tissue_information(path):
    f = open(path)
    web_page = f.read()
    f.close()
    #extract tissue
    tissue = re.search('<th>tissue</th><td>[^<]*',web_page)
    sample_composition = re.search('<th>sample composition</th><td>[^<]*',web_page)
    brain_region = re.search('<th>brain region</th><td>[^<]*',web_page)
    source_name = re.search('<tr><th>source name</th><td>[^<]*',web_page)
    cell_type = re.search('<th>cell type</th><td>[^<]*',web_page)
    organoid_culture = re.search('<th>organoid culture</th><td>[^<]*',web_page)
    starting_cell = re.search('<th>starting cell</th><td>[^<]*',web_page)
    description = re.search('<dt>Description</dt><dd class="wide"><p class="first">[^<]*',web_page)
    if tissue != None:
        tissue = tissue.group()[19:]
    if sample_composition != None:
        sample_composition = sample_composition.group()[31:]
    if brain_region != None:
        brain_region = brain_region.group()[25:]
    if source_name != None:
        source_name = source_name.group()[28:]
    if cell_type != None:
        cell_type = cell_type.group()[22:]
    if organoid_culture != None:
        organoid_culture = organoid_culture.group()[29:]
    if starting_cell != None:
        starting_cell = starting_cell.group()[26:]
    if description != None:
        description = description.group()[54:]
    information = [tissue,sample_composition,brain_region,source_name,cell_type,organoid_culture,starting_cell,description]
    return information

#folder path
folder_path = '/home/wangys/isoform-network/Crawling/Sample'
lisr_dir = os.listdir(folder_path)

#extract tissue information and storage in a dict
tissue_dataframe = pd.DataFrame(columns=['tissue','sample_composition','brain_region','source_name','cell_type','organoid_culture','starting_cell','description'])
for SampleID in lisr_dir:
    path = '/home/wangys/isoform-network/Crawling/Sample/' + SampleID
    information = extract_tissue_information(path)
    tissue_dataframe = tissue_dataframe.append(pd.DataFrame({'tissue':information[0],'sample_composition':information[1],'brain_region':information[2],'source_name':information[3],'cell_type':information[4],'organoid_culture':information[5],'starting_cell':information[6],'description':information[7]},index=[SampleID]))


tissue_dataframe.to_csv('/home/wangys/isoform-network/Crawling/Sample_tissue.csv')