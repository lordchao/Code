setwd("F:/Project/Data/????????????")
set.seed(2016)
library(h2o)
h2o.init()
y.trian<- read.table("training.txt", header = T, row.names = 2)
y.test<- read.table("test.txt", header = T, row.names = 2)
trainingset<-read.csv("training.csv", header = T, row.names = 1)
testset<-testset<-read.csv("test.csv", header = T, row.names = 1)
#combine daraframe
dataset<-rbind(testset[1:23,],trainingset)
dataset<- dataset[order(row.names(dataset)),]
y<-rbind(y.test[1:23,],y.trian)
y<- y[order(row.names(y)),]
#fix name
rownames(dataset) <- sub("./", " ", rownames(dataset))
rownames(dataset) <- sub(".hin", " ", rownames(dataset))
#bind
dataset<-cbind(y$Exp,dataset)

#spilt dataset
train.data<-dataset[1:86,]
test.data<-dataset[87:109,]

#normalization
maxs<-apply(dataset,2,max)
mins<-apply(dataset,2,min)
norm.data<-as.data.frame(scale(dataset,center = mins,scale = maxs-mins))
#delete column

norm.data<-norm.data[, colSums(is.na(norm.data)) != nrow(norm.data)]
norm.data<-norm.data[,colSums(norm.data==0)<4]
norm.data<-norm.data[,sapply(norm.data,is.numeric)]

#PCA
y.norm<-norm.data[,1]
norm.data<-norm.data[,-1]
pca<-prcomp(norm.data)
norm.pca<-as.data.frame(pca$x[,1:10])   #PCA??????????????????
norm.pca<-cbind(y.norm,norm.pca)

train.pca<-norm.pca[1:86,]
test.pca<-norm.pca[87:109,]
#??????
train.data.hex<-as.h2o(train.data)
test.data.hex<-as.h2o(test.data)
test.pca.hex<-as.h2o(test.pca)
train.pca.hex<-as.h2o(train.pca)
#??????
drug.dl.train <- h2o.deeplearning(x = 2:1876, y = 1,
                                  training_frame = train.data.hex,
                                  validation_frame=train.data.hex,
                                  hidden = (c(200,200,200)))
drug.dl.test <- h2o.deeplearning(x = 2:1876, y = 1, training_frame = test.data.hex)

drug.pca.train <- h2o.deeplearning(x = 2:11, y = 1,
                                   training_frame = test.pca.hex,
                                   hidden = (c(200,200,200)))
drug.pca.test <- h2o.deeplearning(x = 2:11, y = 1, training_frame = test.pca.hex,nfolds = 5)

#gbm
gbm <- h2o.gbm(x = 2:1876, y = 1, training_frame = train.data.hex, validation_frame=test.data.hex, ntrees=500, learn_rate=0.01, score_each_iteration = F)
gbm <- h2o.gbm(x = 2:11, y = "y", training_frame = train.pca.hex, validation_frame = valid, ntrees=500, learn_rate=0.01, score_each_iteration = TRUE)
#random forest
rf<-h2o.randomForest(2:11,y=1,training_frame = train.pca.hex,ntrees = 1000)              
#svm
svm<-svm(y.norm~.,data = train.pca)

#??????
predictions.svm<-predict(svm,test.pca)
predictions.train <- h2o.predict(drug.dl.train, test.data.hex)
predictions.test <- h2o.predict(drug.dl.test, test.data.hex)
predictions.test.pca <- h2o.predict(drug.pca.test, test.datap.hex)
predictions.train.pca <- h2o.predict(drug.pca.train, test.datap.hex)

as.data.frame(predictions.train) 
performance = h2o.performance(model = drug.dl.train)  
print(performance) 

#auc
h2o.auc(drug.dl.train, train = FALSE, valid = FALSE, xval = FALSE)
