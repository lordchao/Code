library("RMySQL")


setwd("/home/chao/Documents/Project/Data/Compounds Data/")
set.seed(2016)
chembl = dbConnect(MySQL(), user='root', password='liuchao', dbname='chembl_22', host='127.0.0.1')
dbListTables(chembl)
dbListFields(chembl,'compound_properties')
#fetch compound data
compound_properties<-dbSendQuery(chembl,"select * from compound_properties")
compound_properties<-fetch(compound_properties,n=Inf)
head(compound_properties)
#fetch label data
DRUG_MECHANISM<-dbSendQuery(chembl,"select * from drug_mechanism")
DRUG_MECHANISM<-fetch(DRUG_MECHANISM,n=Inf)
DRUG_MECHANISM<-DRUG_MECHANISM[order(DRUG_MECHANISM[,3]),]
DRUG_MECHANISM<-DRUG_MECHANISM[!duplicated(DRUG_MECHANISM[,3]),]
head(DRUG_MECHANISM)
compound_name<-DRUG_MECHANISM$molregno
compound_name
#select rigt compound corresponding to label
compound_properties<-compound_properties[compound_name,]
head(compound_properties)
#creat dataset
dataset<-cbind(DRUG_MECHANISM$direct_interaction,compound_properties)
colnames(dataset)[1]<-"label"
#transform into factor
dataset[,1]<-as.factor(dataset[,2]) 
head(dataset)
#delete column
dataset<-dataset[,-grep("full_molformula",colnames(dataset))]
dataset<-dataset[,-grep("ro3_pass",colnames(dataset))]


#assign value to the char class in properties
as.numeric(dataset[,14])
for(x in 2:2892){
  if(dataset[x,14]=='NEUTRAL') dataset[x,14]<-1
}

#normaliation
maxs<-apply(dataset,2,max)
mins<-apply(dataset,2,min)
dataset<-as.data.frame(scale(dataset,center = mins,scale = maxs-mins))

#creat negative sample
negative<-dataset[sample(nrow(dataset)),]
dataset<-rbind(negative,dataset)
dim(dataset)

#transform into h2o file
train.hex<-as.h2o(train)
test.hex<-as.h2o(test)
validation.hex<-as.h2o(validation)

#build model

library('h2o')
h2o.init()
#three_layer_model
m1 <- h2o.deeplearning(x = predictors, y = response,
                             model_id="dl_three_layers_200_units_classification",
                             training_frame = train.hex,
                             validation_frame=validation.hex,
                             distribution = "bernoulli",
                             hidden = (c(200,200,200)))

#four layers
m2 <- h2o.deeplearning(
  training_frame=train.hex, 
  validation_frame=validation.hex,
  model_id="dl_four_layers_200_units_classification",
  x=predictors,
  y=response,
  hidden=c(200,200,200,200),                  # small network, runs faster
  epochs=1000000,                      # hopefully converges earlier...
  score_validation_samples=10000,      # sample the validation dataset (faster)
  stopping_rounds=2,
  stopping_metric="misclassification", 
  stopping_tolerance=0.01
)
summary(m2)
plot(m2)
#regularized model
m3 <- h2o.deeplearning(
  model_id="dl_model_reglarized_classification", 
  training_frame=train.hex, 
  validation_frame=validation.hex, 
  x=predictors, 
  y=response, 
  overwrite_with_best_model=F,    ## Return the final model after 10 epochs, even if not the best
  hidden=c(200,200,200,200),          ## more hidden layers -> more complex interactions
  epochs=10,                      ## to keep it short enough
  score_validation_samples=10000, ## downsample validation set for faster scoring
  score_duty_cycle=0.025,         ## don't score more than 2.5% of the wall time
  adaptive_rate=F,                ## manually tuned learning rate
  rate=0.01, 
  rate_annealing=2e-6,            
  momentum_start=0.2,             ## manually tuned momentum
  momentum_stable=0.4, 
  momentum_ramp=1e7, 
  l1=1e-5,                        ## add some L1/L2 regularization
  l2=1e-5,
  max_w2=10                       ## helps stability for Rectifier
) 
summary(m3)
#  
#    
####Checkpointing
#Let's continue training the manually tuned model from before, for 2 more epochs. Note that since many important parameters such as `epochs, l1, l2, max_w2, score_interval, train_samples_per_iteration, input_dropout_ratio, hidden_dropout_ratios, score_duty_cycle, classification_stop, regression_stop, variable_importances, force_load_balance` can be modified between checkpoint restarts, it is best to specify as many parameters as possible explicitly.
#
max_epochs <- 12 ## Add two more epochs
m_cont <- h2o.deeplearning(
  model_id="dl_model_tuned_continued_classification", 
  checkpoint="dl_model_tuned", 
  training_frame=train.hex, 
  validation_frame=validation.hex, 
  x=predictors, 
  y=response, 
  hidden=c(200,200,200),          ## more hidden layers -> more complex interactions
  epochs=max_epochs,              ## hopefully long enough to converge (otherwise restart again)
  stopping_metric="logloss",      ## logloss is directly optimized by Deep Learning
  stopping_tolerance=1e-2,        ## stop when validation logloss does not improve by >=1% for 2 scoring events
  stopping_rounds=2,
  score_validation_samples=10000, ## downsample validation set for faster scoring
  score_duty_cycle=0.025,         ## don't score more than 2.5% of the wall time
  adaptive_rate=F,                ## manually tuned learning rate
  rate=0.01, 
  rate_annealing=2e-6,            
  momentum_start=0.2,             ## manually tuned momentum
  momentum_stable=0.4, 
  momentum_ramp=1e7, 
  l1=1e-5,                        ## add some L1/L2 regularization
  l2=1e-5,
  max_w2=10                       ## helps stability for Rectifier
) 
summary(m_cont)
plot(m_cont)
####Cross-Validation   
dlmodel <- h2o.deeplearning(
  x=predictors,
  y=response, 
  training_frame=train.hex,
  hidden=c(10,10),
  epochs=1,
  nfolds=5,
  fold_assignment="Modulo" # can be "AUTO", "Modulo", "Random" or "Stratified"
)
dlmodel
##try other model
#svm
model.svm<-svm(subset(train,select = 3:1877),subset(train,select = 1),
               na.action = na.omit,
               nfolds= 86,
               momentum_start=0.5,
               cost = 0.5,
               coef0 = 1,
               scale = TRUE,
               epsilon = 0.01)
#build basic binomial glm model
m4 <- h2o.glm(x = predictors,
                    y = response,
                    training_frame = train.hex,
                    validation_frame = validation.hex,
                    model_id = "glm",
                    family = "gaussian",
                    nfold=10,
                    lambda_search = TRUE) 
summary(m4)
h2o.performance(model = m4,
                newdata = test.hex)
h2o.auc(m4,train = TRUE)
# Random Forest
m5 <- h2o.randomForest(x = predictors,
                            y = response,
                            training_frame = train.hex,
                            model_id = "rf",
                            validation_frame = validation.hex,  #only used if stopping_rounds > 0
                            ntrees = 100,
                            seed = 1,
                            nfolds = 10)
summary(m5)
