trainw1 <- read.csv("C:/Users/sabee/Documents/R/trainset1.csv", stringsAsFactors = FALSE)
str(trainw1)
test <- read.csv("C:/Users/sabee/Documents/R/testset1.csv", stringsAsFactors = FALSE)
str(test)
trainw1$selection <- factor(trainw1$selection)

str(trainw1$selection)
table(trainw1$selection)
library(kernlab)
svmmdl <- ksvm(selection  ~ ., data = trainw1,kernel = "vanilladot")
predictions <- predict(svmmdl, test[-21])

head(predictions)

table(predictions, test$selection)



# look only at agreement vs. non-agreement
# construct a vector of TRUE/FALSE indicating correct/incorrect predictions
agreement <- predictions == test$selection
table(agreement)
prop.table(table(agreement))

## Step 5: Improving model performance ----
#set.seed(12345)
svmmdl2 <- ksvm(selection  ~ ., data = trainw1,kernel = "rbfdot")
predictions <- predict(svmmdl2, test[-21])


agreement <- predictions == test$selection
table(agreement)
prop.table(table(agreement))



svmmdl3 <- ksvm(selection  ~ ., data = trainw1,kernel = "polydot")
predictions <- predict(svmmdl3, test[-21])


agreement <- predictions == test$selection
table(agreement)
prop.table(table(agreement))

svmmdl4 <- ksvm(selection  ~ ., data = trainw1,kernel = "tanhdot")
predictions <- predict(svmmdl4, test[-21])


agreement <- predictions == test$selection
table(agreement)
prop.table(table(agreement))

####################################
str(svmBag)

svmBag$fit

bagctrl <- bagControl(fit = svmBag$fit,
                      predict = svmBag$pred,
                      aggregate = svmBag$aggregate)

# fit the bagged svm model
set.seed(300)
svmbag <- train(selection ~ ., data = trainw1, "bag",trControl = ctrl, bagControl = bagctrl)

svmbag
