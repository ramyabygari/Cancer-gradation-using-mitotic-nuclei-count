trainw1 <- read.csv("C:/Users/sabee/Documents/R/trainset1.csv", stringsAsFactors = FALSE)
str(trainw1)
test <- read.csv("C:/Users/sabee/Documents/R/testset1.csv", stringsAsFactors = FALSE)
str(test)
trainw1$selection <- factor(trainw1$selection)

str(trainw1$selection)
table(trainw1$selection)



library(randomForest)
rf <- randomForest(selection ~ ., data = trainw1,ntree=1000,mtry=sqrt(20))
library(caret)
ctrl <- trainControl(method = "LGOCV",p=0.75)
#number = 25, repeats = 10
# auto-tune a random forest
grid_rf <- expand.grid(.mtry = c(2, 4,8,16,18))
m_rf <- train(selection ~ ., data = trainw1, method = "rf",metric = "Kappa", trControl = ctrl,tuneGrid = grid_rf)
#m_rf
predtn<-predict(m_rf,test[-21])
#predtn <-as.data.frame( predtn, drop=false);
#test$selection<-as.data.frame( test$selection, drop=false);
confusionMatrix(predtn, test$selection,positive = '1')
library(gmodels)
CrossTable(x = test$selection, y = predtn,prop.chisq = FALSE);
library(caret)

predtn <-as.factor( predtn);
test$selection<-as.factor( test$selection);

rec<-sensitivity(predtn, test$selection, positive = '1')
spec<-specificity(predtn, test$selection, negative = '0')
prec<-posPredValue(predtn, test$selection, positive = '1')
f <- (2 * prec * rec) / (prec + rec)
library(vcd)
Kappa(table(test$selection, predtn))
