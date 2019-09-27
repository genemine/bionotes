import re
import os
import pandas as pd

#extract SampleID from html file and return
def extract(path):
    #read html file
    f = open(path)
    html = f.read()
    f.close()
    #extract SampleID[abc]
    Sample_ID = re.search(r'title="Link to BioSample">[^<]+',html)
    'title = "Link to BioSample" > SAMN03658256'
    Sample_ID = Sample_ID.group()[26:]
    return Sample_ID

#make a list of specified folder
path = '/home/wangys/isoform-network/Crawling/my_tissue_html/'
list_dir = os.listdir(path)

#reserve the relationship between ExperimentID and SampleID in a dict
dict = {}

for filename in list_dir:
    SampleID = extract(path+filename)
    i = re.search('[^.]*', filename)
    ExperimentID = i.group()
    dict.update({ExperimentID:[SampleID]})

#write the dict to SampleList.txt
data = pd.DataFrame(dict)
data.to_csv('/home/wangys/isoform-network/Crawling/Sample/SampleList.txt',index=False)



###############################
###############################
###down Sample web page

#path of storage
save_path = '/home/wangys/isoform-network/Crawling/Sample'

#path of download
down_path = 'https://www.ncbi.nlm.nih.gov/biosample/'
#download
for SampleID in dict.values():
    command = 'wget -P'+' '+save_path+' '+down_path+SampleID[0]
    os.system(command)