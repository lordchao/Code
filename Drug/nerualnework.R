setwd("F:/Project/NCI60_SuperLearner-master")
install.packages("caret")
install.packages("h2o")

library(caret)

VariantTable <- read.csv(file.path("Data", "VariantTableByGene.csv"), row.names = 1)

## Load drug/agent GI50 data
DrugDat <- read.csv(file.path("Data", "AOD_IOA_GI50.csv"))

# find overlap in cell lines
setdiff(rownames(VariantTable), unique(DrugDat$CELL))
setdiff(unique(DrugDat$CELL), rownames(VariantTable))

# fix names in exome data
rownames(VariantTable) <- sub("A549_ATCC", "A549/ATCC", rownames(VariantTable))
rownames(VariantTable) <- sub("COLO-205", "COLO 205", rownames(VariantTable))
rownames(VariantTable) <- sub("DU145", "DU-145", rownames(VariantTable))
rownames(VariantTable) <- sub("HCC2998", "HCC-2998", rownames(VariantTable))
rownames(VariantTable) <- sub("HL-60", "HL-60(TB)", rownames(VariantTable))
rownames(VariantTable) <- sub("Hs_578T", "HS 578T", rownames(VariantTable))
rownames(VariantTable) <- sub("IGR-OV1", "IGROV1", rownames(VariantTable))
rownames(VariantTable) <- sub("K562", "K-562", rownames(VariantTable))
rownames(VariantTable) <- sub("LOX_IMVI", "LOX IMVI", rownames(VariantTable))
rownames(VariantTable) <- sub("MDA-MB-231", "MDA-MB-231/ATCC", rownames(VariantTable))
rownames(VariantTable) <- sub("NCI-ADR-RES", "NCI/ADR-RES", rownames(VariantTable))
rownames(VariantTable) <- sub("RXF-393", "RXF 393", rownames(VariantTable))

# restrict drug data to overlap with exome data
DrugDat <- DrugDat[which(DrugDat$CELL %in% rownames(VariantTable)), ]

GI50Wide <- reshape(DrugDat[, c("NSC", "CELL", "NLOGGI50")], direction = "wide", timevar = "NSC", idvar = "CELL")
colnames(GI50Wide) <- sub("NLOGGI50.", "NSC", colnames(GI50Wide))

# use cell line name as rowname
rownames(GI50Wide) <- GI50Wide[, 1]
# remove cell line name, and line up cell lines with Variant Table
GI50Wide <- GI50Wide[rownames(VariantTable), -1]
all.equal(rownames(VariantTable), rownames(GI50Wide))

## filter VariantTable
CountVar <- colSums(VariantTable)
head(CountVar,5)
VariantTable_sub <- as.data.frame(VariantTable[, CountVar > 5])
dim(VariantTable_sub)

## clean variable names in variant table
colnames(VariantTable_sub) <- make.names(colnames(VariantTable_sub)) 

#分割数据集
NSC752<-cbind(GI50Wide[,1],VariantTable_sub)
names(NSC752)[1]<-"GI50"
NSC752<-NSC752[!is.na(NSC752$GI50),]

inTrain<-createDataPartition(NSC752$GI50,times = 1,p=0.8,list = FALSE)
train_NSC752<- NSC752[inTrain,]
test_NSC752<- NSC752[-inTrain,]

#NERUAL NETWORK
library(nnet)
nn<-nnet(train_NSC752$GI50~.,train_NSC752,size = 10,decay = 0,01,maxit = 1000,trace=F)
# SVM
library(e1071)
sv<-svm(train_NSC752$GI50~.,data=train_NSC752,gamma = 0.001402525,cost = 0.1 )
tune.out<-tune(svm,train_NSC752$GI50~.,train_NSC752$GI50,data = train_NSC752,ranges = list(epsilon = seq(0, 1, 0.1), cost = 2^(2:9)))
bestmod<-tune.out$best.model

s.pred<-predict(sv,test_NSC752)

RMSE(s.pred,NSC752$GI50)




