
rm(list=ls())
library(vioplot)
library(sva)


# folder setting and file list
setwd("~/Desktop/combat_usage")
load("testdata.Rdata")
ir=as.matrix(ir) # rows are genes, columns are samples
prin(dim(ir))


# COMBAT: here lab is treated as the batch variable.

batchVariable = meta$lab
mod1 = model.matrix(~1, data=meta)
ir_corrected = ComBat(dat=ir, batch=batchVariable, mod=mod1)






