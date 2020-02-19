
#++++++++++++++++++++++++++++++++EXPLORATORY ANALYSIS++++++++++++++++++++++++++++++


library(tidyverse)
library(ggcorrplot)
library(lattice)
library(tidyr)
library(caret)
library(e1071)
library(ggplot2)
library(rpart.plot)

breastCancer<- read.csv("breast-cancer-wisconsin.csv")

colnames(breastCancer) <- c("ID", "Clump_thickness", "Uniformity_of_Cell_Size", "Uniformity_of_Cell_Shape", 
                            "Marginal_Adhesion", "Single_Epithelial_Cell_Size", 
                            "Bare_Nuclei", "Bland_Chrobreast_cancer_wisconsinmatin", "Normal_Nucleoli", "Mitosis", "Class")


breastCancer$Class <- ifelse(breastCancer$Class == "2", "benign",
                             ifelse(breastCancer$Class == "4", "malignant", NA))

breastCancer[breastCancer == "?"] <- NA

str(breastCancer)
# how many NAs are in the data
length(which(is.na(breastCancer)))

# impute missing data


breastCancer[,2:10] <- apply(breastCancer[, 2:10], 2, function(x) as.numeric(as.character(x)))
library(mice)
dataset_impute <- mice(breastCancer[,2:10],  print = FALSE)
breastCancer <- cbind(breastCancer[, 11, drop = FALSE], mice::complete(dataset_impute, 1))

breastCancer$ID=NULL
str(breastCancer)
breastCancer$Class <- as.factor(breastCancer$Class)

# how many benign and malignant cases are there?
summary(breastCancer$Class)

str(breastCancer)

bc_data_gather <- breastCancer %>%
  gather(measure, value, Clump_thickness:Mitosis)

ggplot(bc_data_gather, aes(x = value, fill = Class, color = Class)) +
  geom_density(alpha = 0.3, size = 1) +
  geom_rug() +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap( ~ measure, scales = "free_y", ncol = 3)
#==================== Correlation================================
library(dplyr)
library(GGally)
bc<-breastCancer

bc$Class<-as.numeric(bc$Class)

ggcorr(bc,  method = c("pairwise", "pearson"),nbreaks =7, palette = "RdBu", geom = "tile",
       midpoint = 0, legend.position = "left", label_size = 2, label_alpha = "0", layout.exp = 0)



#====================================Decision Trees========================================
set.seed(123,sample.kind = "Rejection")
ind <- sample(2, nrow(breastCancer),replace = TRUE, prob = c(0.7,0.3))
train <- breastCancer[ind==1,]
head(train)
test <- breastCancer[ind==2,]
head(test)

library(rpart)

#####################by information gain ##################
latlontree <- rpart(Class ~.,train,parms = list(prior = c(.65,.35), split = "information"),maxdepth=5)
set.seed(123)
latlontree
print(latlontree)
plot(latlontree)
text(latlontree, digits = 2)
predicted=predict(latlontree,test, "class") 
summary(predicted)

conf<-table(predicted, test$Class)
conf
missclassificationerror<-1-(sum(diag(conf))/sum(conf))
missclassificationerror
#Accuracy caculation

confusionMatrix(conf)
rpart.plot(latlontree, extra=2,
           main="Breast Cancer Min Error classification tree using Information Gain")



######################entropy####################################
latlontree <- rpart(Class ~.,train,parms = list(prior = c(.65,.35), split = "entropy"),maxdepth=5)
set.seed(123)
latlontree
print(latlontree)
plot(latlontree)
text(latlontree, digits = 2)
predicted=predict(latlontree,test, "class") 
summary(predicted)

conf<-table(predicted, test$Class)
conf
missclassificationerror<-1-(sum(diag(conf))/sum(conf))
missclassificationerror
#Accuracy caculation

confusionMatrix(conf)




#################### Gini ##################
latlontree <- rpart(Class ~.,train,parms = list(prior = c(.65,.35), split = "gini"),maxdepth=5)
set.seed(123)
latlontree
print(latlontree)
plot(latlontree)
text(latlontree, digits = 2)
predicted=predict(latlontree,test, "class") 
summary(predicted)

conf1<-table(predicted, test$Class)
conf1
missclassificationerror<-1-(sum(diag(conf1))/sum(conf1))
missclassificationerror
#Accuracy caculation

confusionMatrix(conf1)

rpart.plot(latlontree, extra=2,
           main="Breast Cancer Min Error classification tree using Gini")







#======================================Random Forest===========================================
library(dplyr)

library(corrplot)

dim(breastCancer)

#check for unique values in the dataset
apply(breastCancer,2,function(x)length(unique(x)))


#Split dataset into 70% train and 30% test data
set.seed(123)
library(caTools)
ind<-sample.split(Y=breastCancer$Class,SplitRatio = 0.70)
traindf<-breastCancer[ind,]
str(traindf)
dim(traindf)
testdf<-breastCancer[!ind,]
dim(testdf)


#Find best mtry
library(randomForest)
bestmtry <- tuneRF(traindf , traindf$Class , ntreeTry = 200 , stepFactor = 1.2 , improve = 0.01, trace = T , plot = T)
bestmtry

# Fit the model
model<-randomForest(Class~. , data = traindf, mtry=3, ntree=20)
model
print(model)


# Find importance of the model and predictors 
importance(model)
varImpPlot(model)

#Prediction on Test dataset
predictdf<-predict(model,testdf,type ='class')
predictdf


#Table metric
t<-table(predictions=predictdf , actual=testdf$Class)
t

confusionMatrix(t)



