trainnew4 <- read.csv("C:/Users/train4setxl.csv", stringsAsFactors = FALSE)
str(trainnew4)
testnew2 <- read.csv("C:/Users/sabee/Documents/R/datanew/secondtestdataxl.csv", stringsAsFactors = FALSE)
str(testnew2)
trainnew4$selection <- factor(trainnew4$selection)
str(trainnew4$selection)
table(trainnew4$selection)
m <- train(selection ~ ., data = trainnew4, method = "C5.0")

# summary of tuning results


# apply the best C5.0 candidate model to make predictions
p <- predict(m, testnew2[-21])
table(p, testnew2$selection)
ctrl <- trainControl(method = "cv", number = 10,
                     selectionFunction = "oneSE")

# use expand.grid() to create grid of tuning parameters
grid <- expand.grid(.model = "tree",
                    .trials = c(1, 5, 10, 15, 20, 25, 30, 35),
                    .winnow = "FALSE")
m_tre<- train(selection ~ ., data = trainnew4, method = "C5.0",metric = "Kappa", trControl = ctrl,tuneGrid = grid)
p <- predict(m_tre, testnew2[-21])
confusionMatrix(p, testnew2$selection,positive = '1')
predtn <-as.factor( p);
testnew2$selection<-as.factor( testnew2$selection);
rec<-sensitivity(predtn, testnew2$selection, positive = '1')
spec<-specificity(predtn, testnew2$selection, negative = '0')
prec<-posPredValue(predtn, testnew2$selection, positive = '1')
f <- (2 * prec * rec) / (prec + rec)
#####################################################################3
library(ipred)

mybag <- bagging(selection ~ ., data = trainnew4, nbagg = 25)
credit_pred <- predict(mybag, testnew2[-21])
table(credit_pred, testnew2$selection)

# estimate performance of ipred bagged trees

library(caret)

ctrl <- trainControl(method = "cv", number = 10)
m2<-train(selection ~ ., data = trainnew4, method = "treebag",
      trControl = ctrl)
credit_pred <- predict(m2, testnew2[-21])
table(credit_pred, testnew2$selection)

confusionMatrix(credit_pred, testnew2$selection,positive = '1')
predtn <-as.factor( credit_pred);
testnew2$selection<-as.factor( testnew2$selection);
rec<-sensitivity(predtn, testnew2$selection, positive = '1')
spec<-specificity(predtn, testnew2$selection, negative = '0')
prec<-posPredValue(predtn, testnew2$selection, positive = '1')
f <- (2 * prec * rec) / (prec + rec)
## Using AdaBoost.M1
library(adabag)
m_adaboost <- boosting(selection ~ ., data = trainnew4)
p_adaboost <- predict(m_adaboost, testnew2[-21])
p_adaboost$confusion
table(p_adaboost, testnew2$selection)
p_adaboost<-as.factor(p_adaboost)
confusionMatrix(p_adaboost, testnew2$selection,positive = '1')
adaboost_cv <- boosting.cv(selection ~ ., data = trainnew4)
p_adaboost <- predict(adaboost_cv, testnew2[-21])
adaboost_cv$confusion

# calculate kappa
library(vcd)
Kappa(adaboost_cv$confusion)
