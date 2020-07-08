
rm(list=ls())
library(vioplot)
library(sva)


# folder setting and file list
setwd("~/Desktop/combat_usage")
load("testdata.Rdata")



# COMBAT: lab+pop. Correct ir to ir_2 (genes in row, samples in column)
batch1 = meta$lab
mod1 = model.matrix(~1+pop, data=meta)
ir_1 = ComBat(dat=ir, batch=batch1, mod=mod1)



batch2 = meta$pop
mod2 = model.matrix(~1, data=meta)
ir_2 = ComBat(dat=ir_1, batch=batch2, mod=mod2)
ir_2 = round(ir_2,4)  # corrected data to use




