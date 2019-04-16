# Author: Camilla Fløjgaard
# Script for analysis for paper: 
# Title:  Predicting provenance of forensic soil samples: 
#         Linking soil to ecological habitats by metabarcoding and supervised classification
# Submitted to Proc B 

Data <- read.table('SoilTackerDataforPublication.txt', header = TRUE)

library(MASS)
library(mgcv)

# Create a result file: 
Prediction <- data.frame(Data[,1])

#### Predicting environmental gradients ####
primers <- c("fung","insec","eukar","plant")

# Fetility
# Model selection: 
Result.list.N <- list()
Summary.list.N <- list()
for (i in seq_along(primers)){
  pred.vector <- paste(primers[i],c(1:4),sep="")
  Model <- stepAIC(lm(as.formula(paste("EIV_N~", paste(pred.vector, collapse = "+"))), data=Data), direction="backward")
  best.vars <- as.vector(attr(Model$coefficients,"names")[2:length(Model$coefficients)])
  Result.list.N[[i]] <- lm(as.formula(paste("EIV_N~", paste(best.vars, collapse = "+"))), data=Data)
  Summary.list.N[[i]] <- summary(lm(as.formula(paste("EIV_N~", paste(best.vars, collapse = "+"))), data=Data))
}

r.squared.vector <- c(Summary.list.N[[1]]$r.squared,Summary.list.N[[2]]$r.squared,Summary.list.N[[3]]$r.squared,Summary.list.N[[4]]$r.squared)
primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.N[[which(r.squared.vector[] == max(r.squared.vector))]]
# Check residuals of the best model:
par(mfrow=c(1,1))
plot(Result.list.N[[4]]$fitted,Result.list.N[[4]]$residuals)

# plant1 is non-significant (see also AIC), remove plant1 from model and check again:
test <- lm(EIV_N~plant2+plant3+plant4, data=Data)
Summary.list.N.test <- summary(test)
Summary.list.N.test$r.squared
plot(test$fitted,test$residuals)

# So... 
Result.list.N[[4]] <- lm(EIV_N~plant2+plant3+plant4, data=Data)
Summary.list.N[[4]] <- summary(test)
r.squared.vector[4] <- Summary.list.N.test$r.squared

# Predict values using leave-one-out and best model
Prediction$N <- NA
Prediction$NSE <- NA
pred.vector <- as.vector(names(Result.list.N[[4]]$coefficients)[2:length(Result.list.N[[4]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_N~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$N[i] <- Pred[1]
  Prediction$NSE[i] <- Pred[1]-Pred[2]
}

# Moisture
# Model selection
Result.list.M <- list()
Summary.list.M <- list()
for (i in seq_along(primers)){
  pred.vector <- paste(primers[i],c(1:4),sep="")
  Model <- stepAIC(lm(as.formula(paste("EIV_M~", paste(pred.vector, collapse = "+"))), data=Data), direction="backward")
  best.vars <- as.vector(attr(Model$coefficients,"names")[2:length(Model$coefficients)])
  Result.list.M[[i]] <- lm(as.formula(paste("EIV_M~", paste(best.vars, collapse = "+"))), data=Data)
  Summary.list.M[[i]] <- summary(lm(as.formula(paste("EIV_M~", paste(best.vars, collapse = "+"))), data=Data))
}

r.squared.vector <- c(Summary.list.M[[1]]$r.squared,Summary.list.M[[2]]$r.squared,Summary.list.M[[3]]$r.squared,Summary.list.M[[4]]$r.squared)
primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.M[[which(r.squared.vector[] == max(r.squared.vector))]]
# Check residuals of the best model:
plot(Result.list.M[[1]]$fitted,Result.list.M[[1]]$residuals)

# eukar4 is non-significant (see also AIC), remove eukar4 from model and check again:
test <- lm(EIV_M~eukar2+eukar3, data=Data)
Summary.list.M.test <- summary(test)
Summary.list.M.test$r.squared
plot(test$fitted,test$residuals)

# So... 
Result.list.M[[3]] <- lm(EIV_M~eukar2+eukar3, data=Data)
Summary.list.M[[3]] <- summary(test)
r.squared.vector[3] <- Summary.list.N.test$r.squared

# Predict values using leave-one-out and best model
Prediction$M <- NA
Prediction$MSE <- NA
pred.vector <- as.vector(names(Result.list.M[[3]]$coefficients)[2:length(Result.list.M[[3]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_M~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$M[i] <- Pred[1]
  Prediction$MSE[i] <- Pred[1]-Pred[2]
}

# pH
# model selection: 
Result.list.R <- list()
Summary.list.R <- list()
for (i in seq_along(primers)){
  pred.vector <- paste(primers[i],c(1:4),sep="")
  Model <- stepAIC(lm(as.formula(paste("EIV_Ph~", paste(pred.vector, collapse = "+"))), data=Data), direction="backward")
  best.vars <- as.vector(attr(Model$coefficients,"names")[2:length(Model$coefficients)])
  Result.list.R[[i]] <- lm(as.formula(paste("EIV_Ph~", paste(best.vars, collapse = "+"))), data=Data)
  Summary.list.R[[i]] <- summary(lm(as.formula(paste("EIV_Ph~", paste(best.vars, collapse = "+"))), data=Data))
}

r.squared.vector <- c(Summary.list.R[[1]]$r.squared,Summary.list.R[[2]]$r.squared,Summary.list.R[[3]]$r.squared,Summary.list.R[[4]]$r.squared)
primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.R[[which(r.squared.vector[] == max(r.squared.vector))]]

# Check residuals of the best model:
plot(Result.list.R[[1]]$fitted,Result.list.R[[1]]$residuals)

# Predict values using leave-one-out and best model: 
Prediction$R <- NA
Prediction$RSE <- NA
pred.vector <- as.vector(names(Result.list.R[[1]]$coefficients)[2:length(Result.list.R[[1]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_Ph~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$R[i] <- Pred[1]
  Prediction$RSE[i] <- Pred[1]-Pred[2]
}

# Light
# Model selection
Result.list.L <- list()
Summary.list.L <- list()
for (i in seq_along(primers)){
  pred.vector <- paste(primers[i],c(1:4),sep="")
  Model <- stepAIC(lm(as.formula(paste("EIV_L~", paste(pred.vector, collapse = "+"))), data=Data), direction="backward")
  best.vars <- as.vector(attr(Model$coefficients,"names")[2:length(Model$coefficients)])
  Result.list.L[[i]] <- lm(as.formula(paste("EIV_L~", paste(best.vars, collapse = "+"))), data=Data)
  Summary.list.L[[i]] <- summary(lm(as.formula(paste("EIV_L~", paste(best.vars, collapse = "+"))), data=Data))
}

r.squared.vector <- c(Summary.list.L[[1]]$r.squared,Summary.list.L[[2]]$r.squared,Summary.list.L[[3]]$r.squared,Summary.list.L[[4]]$r.squared)
primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.L[[which(r.squared.vector[] == max(r.squared.vector))]]

# Check residuals of the best model:
plot(Result.list.L[[1]]$fitted,Result.list.L[[1]]$residuals)
# Some non-linear pattern...

#Check if GAM improves the model: 
gam1 <- gam(EIV_L~s(fung1, k=3)+s(fung2, k=3)+s(fung3, k=3)+s(fung4, k=3), data=Data)
par(mfrow=c(2,2))
plot(gam1)
summary(gam1)
lm1 <- lm(EIV_L~fung1+fung2+fung3+fung4, data=Data)
AIC(gam1)
AIC(lm1)
plot(gam1$fitted.values ,gam1$residuals)

# Predict values using leave-one-out and best model: 
Prediction$L <- NA
Prediction$LSE <- NA
pred.vector <- as.vector(names(Result.list.L[[1]]$coefficients)[2:length(Result.list.L[[1]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.gam(gam(EIV_L~s(fung1, k=3)+s(fung2, k=3)+s(fung3, k=3)+s(fung4, k=3), data=Data[-i,]),
                             Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$L[i] <- Pred[1]
  Prediction$LSE[i] <- Pred[2] ## SE
}

head(Prediction)
write.table(Prediction, file="Data/PredictionEllenbergValuesGAM2012218.txt")


#### Predicting habitat characteristics ####

# function for step-selection of best model: 
# Author: Jesper Moeslund
stepCV.qda <- function(resp.var, pred.vector, data){
  goal.reached <- FALSE #used to flag when no more variables should be removed in the while loop below
  cur.preds <- NULL
  preds.not.added.currently <- pred.vector
  add.these.preds.to.get.best.model <- NULL #Variables that should be kicked out will be added to this and this vector will be returned
  max.no.of.preds <- min(table(data[,resp.var])) #There can be no more predictors than the "minimum number of observations minus 1" in the smallest class in the response variable (qda limitation)
  print(paste0("Max no of preds to test is: ",max.no.of.preds-1 , " because the no of observation the class having the lowest no of observations is: ",max.no.of.preds))
  while(!goal.reached){
    if(!is.null(cur.preds)){
      cur.form <- as.formula(paste(resp.var, "~", paste(cur.preds,collapse = "+")))
      print(paste("Now fitting these predictors as base model: ", paste(cur.preds,collapse = "+")))
      cur.base.model <- qda(cur.form, data = data, CV=TRUE)
      cr.table <- table(cur.base.model$class, data[,resp.var])
      cur.preds.best.cv <- sum(diag(prop.table(cr.table))) #percent correctly classified
    }
    else {cur.preds.best.cv <- -99999999}
    pred.to.add <- ""
    for (pred in preds.not.added.currently){ #Comparing cv performance of models with one predictor thrown out at a time
      print(paste("Now testing if ", pred, " should be added"))
      #new.preds <- cur.preds[!cur.preds %in% pred]
      min.var.in.groups <- min(aggregate(data[,pred], by=list(Category=data[,resp.var]), FUN=sd)[,2]) #Calculating the minimum sd value within each group. If zero no variation in that group and the corresponding predictor cannot be tested meaningfully (and isn't a candidate at all!)
      #print(min.var.in.groups)
      if(min.var.in.groups > 0){
        if(!is.null(cur.preds)){
          new.form <- as.formula(paste(resp.var, "~", paste(cur.preds,collapse = "+"), "+",pred))
        }
        else {
          new.form <- as.formula(paste(resp.var, "~", pred))
        }
        print(toString(new.form))
        new.model <- qda(new.form, data = data, CV = TRUE)
        new.cr.table <- table(new.model$class, data[,resp.var])
        print(new.cr.table)
        new.model.cv.score <- sum(diag(prop.table(new.cr.table))) #percent correctly classified
        print(new.model.cv.score)
        print(cur.preds.best.cv)
        if (new.model.cv.score > cur.preds.best.cv){ #If performance better, mark predictor as candidate for being thrown out. In each step in the while loop only the worst one is kicked out during this for loop
          cur.preds.best.cv <- new.model.cv.score
          pred.to.add <- pred
          print(paste("Candidate for being added: ", pred))
        }
      }
      else{print(paste(pred," does not vary in one or more groups and therefore can't be tested!"))}
    }
    if((pred.to.add != "")){ #add the best predictor to the vector being added and prepare new set of predictors for the new base model
      add.these.preds.to.get.best.model <- append(add.these.preds.to.get.best.model,pred.to.add)
      cur.preds <- append(cur.preds, pred.to.add)
      preds.not.added.currently <- preds.not.added.currently[preds.not.added.currently != pred.to.add]
      if(length(cur.preds) == max.no.of.preds-1){ #Don't run again if number of allowable predictors is exceeded see comment for max.no.of.preds above
        goal.reached <- TRUE
        print("Goal reached, exiting...")
      }
    }
    else{ #if no predictors were candidates for being thrown out (i.e. if model performance was best with all predictors), end the loop
      goal.reached <- TRUE
      print("Goal reached, exiting...")
    }
  }
  return(add.these.preds.to.get.best.model)
}

# High forest
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("HighForest", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("HighForest ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"HighForest"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("HighForest", pred.vector, Data)
    best.qda <- qda(as.formula(paste("HighForest ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"HighForest"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting using leave-one-out
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("HighForest", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$HighForest[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

for (j in 1:4){
  pred.vector <- paste(a$best.vars[j])
  for (i in 1:130){
    column.NA <- which(colnames(Data)==paste(unlist(strsplit(pred.vector, split=" "))))
    if(anyNA(Data[column.NA])) {
      if(is.na(Data[i,column.NA])){
        Prediction[i,2+2*j] <- NA
        Prediction[i,3+2*j] <- NA
      } else{
        best.qda <- qda(as.matrix(Data[-c(i,which(is.na(Data[,column.NA]))),paste(unlist(strsplit(pred.vector, split=" ")))]), total$Hoejskov[-c(i,which(is.na(total[,column.NA])))])
        Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
        Prediction[i,2+2*j] <- Pred$class
        Prediction[i,3+2*j] <- Pred$posterior[2]
      }
    } else {
      best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$HighForest[-i])
      Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
      Prediction[i,2+2*j] <- Pred$class
      Prediction[i,3+2*j] <- Pred$posterior[2]
    }
  }
}

# Forest
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Forest", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Forest ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Forest"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Forest", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Forest ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Forest"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting using leave-one-out
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Forest", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Forest[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

#Agriculture
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Agriculture", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Agriculture ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Agriculture"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Agriculture", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Agriculture ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Agriculture"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Agriculture", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Agriculture[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Coniferous
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Coniferous", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Coniferous ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Coniferous"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Coniferous", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Coniferous ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Coniferous"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Coniferous", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Coniferous[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Beech
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Beech", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Beech ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Beech"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Beech", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Beech ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Beech"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Beech", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Beech[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Oak
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Oak", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Oak ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Oak"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Oak", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Oak ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Oak"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Oak", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Oak[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Willow
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Willow", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Willow ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Willow"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Willow", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Willow ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Willow"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Willow", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Willow[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Heathland
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Heathland", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Heathland ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Heathland"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Heathland", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Heathland ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Heathland"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Heathland", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Heathland[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Dwarfshrubs
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Dwarfshrubs", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Dwarfshrubs ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Dwarfshrubs"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Dwarfshrubs", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Dwarfshrubs ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Dwarfshrubs"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Dwarfshrubs", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Dwarfshrubs[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Atlantic
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Atlantic", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Atlantic ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Atlantic"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Atlantic", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Atlantic ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Atlantic"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Atlantic", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Atlantic[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}

# Jutland
# Model selection: 
Result.list.S <- list()
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  column.NA <- which(colnames(Data)==paste(primers[i],1,sep=""))
  if(anyNA(Data[column.NA])) {
    best.vars <- stepCV.qda("Jutland", pred.vector, Data[-which(is.na(Data[column.NA])),])
    best.qda <- qda(as.formula(paste("Jutland ~ ", paste(best.vars, collapse = "+"))), data = Data[-which(is.na(Data[column.NA])),], CV = TRUE)
    best.cr <- table(best.qda$class, Data[-which(is.na(Data[column.NA])),"Jutland"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  } else {
    best.vars <- stepCV.qda("Jutland", pred.vector, Data)
    best.qda <- qda(as.formula(paste("Jutland ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
    best.cr <- table(best.qda$class, Data[,"Jutland"])
    Result.list.S[[i]] <- sum(diag(prop.table(best.cr)))
  }
}

# Predicting...
Prob.vector <- c(Result.list.S[[1]],Result.list.S[[2]],Result.list.S[[3]],Result.list.S[[4]])
primers[which(Prob.vector[] == max(Prob.vector))]
pred.vector <- paste(primers[which(Prob.vector[] == max(Prob.vector))],c(1:4),sep="")
for (i in 1:130){
  best.vars <- stepCV.qda("Jutland", pred.vector, Data)
  best.qda <- qda(Data[-i,c(paste(best.vars))], Data$Jutland[-i])
  Pred <- predict(best.qda, newdata=Data[i,c(paste(best.vars))])
  Prediction$SProb[i] <- round(Pred$posterior[2], digits=3)
}
