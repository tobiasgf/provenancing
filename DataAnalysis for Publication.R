# Analysis SoilTracker:
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Camilla Fløjgaard
# Date: 23-04-2019
library(readxl)
library(here)
library(MASS)
library(mgcv)


#### Data: environmental data, NMS-ordination-axes and habitat characteristics
Data <- read.table(here::here("data","SoilTrackerData.txt"), header = TRUE) # Data from file "ordination_results.txt" were added to this file manually first

# Create a result file: 
Prediction <- data.frame(Data[,1])

#### Predicting environmental gradients ####
primers <- c("fung","insect","eukar","plant")

########### Fetility ###########
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
x <- primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.N[[which(r.squared.vector[] == max(r.squared.vector))]]
# Check residuals of the best model:
par(mfrow=c(1,1))
plot(Result.list.N[[which(primers==x)]]$fitted,
     Result.list.N[[which(primers==x)]]$residuals)

# Predict values using leave-one-out and best model
Prediction$N <- NA
Prediction$NSE <- NA
pred.vector <- as.vector(names(Result.list.N[[which(primers==x)]]$coefficients)[2:length(Result.list.N[[which(primers==x)]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_N~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$N[i] <- Pred[1]
  Prediction$NSE[i] <- Pred[1]-Pred[2]
}

########### Moisture ##############
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
x<- primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.M[[which(r.squared.vector[] == max(r.squared.vector))]]
# Check residuals of the best model:
plot(Result.list.M[[which(primers==x)]]$fitted,Result.list.M[[which(primers==x)]]$residuals)

# Predict values using leave-one-out and best model
Prediction$M <- NA
Prediction$MSE <- NA
pred.vector <- as.vector(names(Result.list.M[[which(primers==x)]]$coefficients)[2:length(Result.list.M[[which(primers==x)]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_M~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$M[i] <- Pred[1]
  Prediction$MSE[i] <- Pred[1]-Pred[2]
}

########### pH ###########
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
x <- primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.R[[which(r.squared.vector[] == max(r.squared.vector))]]

# Check residuals of the best model:
plot(Result.list.R[[1]]$fitted,Result.list.R[[1]]$residuals)

# Predict values using leave-one-out and best model: 
Prediction$R <- NA
Prediction$RSE <- NA
pred.vector <- as.vector(names(Result.list.R[[which(primers==x)]]$coefficients)[2:length(Result.list.R[[which(primers==x)]]$coefficients)])
for (i in 1:130){
  Pred <- unlist(predict.lm(lm(as.formula(paste("EIV_Ph~", paste(pred.vector, collapse = "+"))), data=Data[-i,]),
                            Data[i,],se.fit=TRUE, interval="predict"))
  Prediction$R[i] <- Pred[1]
  Prediction$RSE[i] <- Pred[1]-Pred[2]
}

########### Light ###########
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
x <-primers[which(r.squared.vector[] == max(r.squared.vector))]
Summary.list.L[[which(r.squared.vector[] == max(r.squared.vector))]]

# Check residuals of the best model:
plot(Result.list.L[[which(primers==x)]]$fitted,Result.list.L[[which(primers==x)]]$residuals)
# Some non-linear pattern...

#Check if GAM improves the model: 
gam1 <- gam(EIV_L~s(fung1, k=3)+s(fung2, k=3)+s(fung3, k=3), data=Data)
par(mfrow=c(2,2))
plot(gam1)
summary(gam1)
lm1 <- lm(EIV_L~fung1+fung2+fung3, data=Data)
AIC(gam1)
AIC(lm1)
plot(gam1$fitted.values ,gam1$residuals)
summary(gam1)$r.sq

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
write.table(Prediction, file=here::here("data","PredictionEllenbergValuesGAM07032019.txt"))

########### Predicting habitat characteristics ###########

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

correctly.predicted <- as.vector(c(1,2,3,4))
true.positives <- as.vector(c(1,2,3,4))
false.neg.sites <- as.vector(c(1,2,3,4))
best.vars <- as.vector(c(1,2,3,4))
a <- data.frame(cbind(correctly.predicted,best.vars,true.positives,false.neg.sites))

########### High forest ###########
PredictionBi <- Data[,c("Site_nr","site_name","HighForest")]

# Model selection:Find the best set of predictors for each primer set 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("HighForest", pred.vector, Data)
  best.qda <- qda(as.formula(paste("HighForest ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"HighForest"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$HighForest==1)],collapse=";")
}
a
a$habchar <- "HighForest" #fung maximizes both corretly.predicted and true positives.
ResultBi <- a
# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$HighForest[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$HFclass[i] <- Pred$class
  PredictionBi$HFpost[i] <- Pred$posterior[2]
}

########### Forest ###########
PredictionBi$Forest <- Data$Forest
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Forest", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Forest ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Forest"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Forest==1)],collapse=";")
}
a #fung maximizes both corretly.predicted and true positives. 
a$habchar <- "Forest"
ResultBi <- rbind(ResultBi,a)
# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Forest[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Foclass[i] <- Pred$class
  PredictionBi$Fopost[i] <- Pred$posterior[2]
}

########### Agriculture ###########
# Model selection: 
PredictionBi$Agriculture <- Data$Agriculture
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Agriculture", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Agriculture ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Agriculture"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Agriculture==1)],collapse=";")
}
a #fung maximizes both corretly.predicted and true positives. 
a$habchar <- "Agriculture"
ResultBi <- rbind(ResultBi,a)
# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[4])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Agriculture[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Agclass[i] <- Pred$class
  PredictionBi$Agpost[i] <- Pred$posterior[2]
}

########### Coniferous ###########
PredictionBi$Coniferous <- Data$Coniferous
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Coniferous", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Coniferous ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Coniferous"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Coniferous==1)],collapse=";")
}
a #No model is good for coniferous - try sequnces....  
a$habchar <- "Coniferous"
ResultBi <- rbind(ResultBi,a)

### Can sequences improve model fit? 
#get sequence data from plants:
Annotated <- read.table(here::here("data","plant_table_final.txt"), header=TRUE, sep="\t")
freq.sums <- colSums(Annotated[1:130])
OTUplanter.frequencies <- Annotated[,1:130]
for (i in 1:130){
  OTUplanter.frequencies[,i] <- Annotated[,i]/freq.sums[i]
}
OTUplanter.frequencies2 <- cbind(Annotated[,137:139],OTUplanter.frequencies)
Data$ConifSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$family %in% c("Pinaceae","Cupressaceae")),4:133])

pred.vector <- c("fung1","ConifSeq")
best.qda <- qda(as.formula(paste("Coniferous ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Coniferous"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Coniferous==1)],collapse=";")

a #does not increase true positives to over 50%, don't use this for prediction.  
a$habchar <- "ConiferousSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)]

########### Beech ###########
# Model selection: 
PredictionBi$Beech <- Data$Beech
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Beech", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Beech ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Beech"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Beech==1)],collapse=";")
}

a #fung maximizes corretly predicted. 
a$habchar <- "Beech"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]

# is model improved by adding fagus sequences?
Data$BeechSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Fagus"),4:133])

pred.vector <- c("fung2","fung4","BeechSeq")
best.qda <- qda(as.formula(paste("Beech ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Beech"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Beech==1)],collapse=";")

a  
a$habchar <- "BeechSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)]
#Does not improve from model without seq  

# Predict (first, run the model selection without seq)
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Coniferous[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Beclass[i] <- Pred$class
  PredictionBi$Bepost[i] <- Pred$posterior[2]
}


########### Oak ###########
PredictionBi$Oak <- Data$Oak
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Oak", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Oak ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Oak"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Oak==1)],collapse=";")
}

a #fung maximizes both corretly.predicted and true pos. 
a$habchar <- "Oak"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]
# is model improved by adding quercus sequences?
Data$OakSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Quercus"),4:133])

pred.vector <- c("fung1","fung2","OakSeq")
best.qda <- qda(as.formula(paste("Oak ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Oak"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Oak==1)],collapse=";")

a  
a$habchar <- "OakSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)] #percentage correctly predicted improved, but not true positives. 

# Predicting using leave-one-out (re-run first model):
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Oak[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Oaclass[i] <- Pred$class
  PredictionBi$Oapost[i] <- Pred$posterior[2]
}

########### Willow ###########
PredictionBi$Willow <- Data$Willow
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Willow", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Willow ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Willow"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Willow==1)],collapse=";")
}

a #fung maximizes both corretly.predicted and true pos.  
a$habchar <- "Willow"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]

#does model improve from adding Salix sequences?
Data$WilSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Salix"),4:133])

pred.vector <- c("fung3","fung1","WilSeq")
best.qda <- qda(as.formula(paste("Willow ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Willow"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Willow==1)],collapse=";")

a  
a$habchar <- "WilSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)] #not improved!

# Predicting using leave-one-out (run model above again)
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Willow[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Wiclass[i] <- Pred$class
  PredictionBi$Wipost[i] <- Pred$posterior[2]
}

########### Heathland ###########
PredictionBi$Heathland <- Data$Heathland
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Heathland", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Heathland ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Heathland"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Heathland==1)],collapse=";")
}

a #no good model...  
a$habchar <- "Heathland"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]
#does model improve from adding Ericaceae sequences?
Data$EricaSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Erica"),4:133])

pred.vector <- c("fung1","EricaSeq")
best.qda <- qda(as.formula(paste("Heathland ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Heathland"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Heathland==1)],collapse=";")

a  
a$habchar <- "EriSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)] 
# is not good enough for prediction


########### Dwarfshrubs ###########
PredictionBi$Dwarfshrubs <- Data$Dwarfshrubs
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Dwarfshrubs", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Dwarfshrubs ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Dwarfshrubs"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Dwarfshrubs==1)],collapse=";")
}

a #fung maximizes both corretly.predicted and true pos.  
a$habchar <- "Dwarfshrubs"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]


# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Dwarfshrubs[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$DSclass[i] <- Pred$class
  PredictionBi$DSpost[i] <- Pred$posterior[2]
}

########### Alder ###########
PredictionBi$Alder <- Data$Alnus
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Alnus", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Alnus ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Alnus"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Alnus==1)],collapse=";")
}

a #fung maximizes both corretly.predicted and true pos.  
a$habchar <- "Alder"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]
#does model improve from adding Alnus sequences?
Data$AlnusSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Alnus"),4:133])

pred.vector <- c("fung1","fung3","AlnusSeq")
best.qda <- qda(as.formula(paste("Alnus ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Alnus"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Alnus==1)],collapse=";")

a  
a$habchar <- "AlnusSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)] #not improved!

# Predicting using leave-one-out (run model above again)
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Alnus[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Alclass[i] <- Pred$class
  PredictionBi$Alpost[i] <- Pred$posterior[2]
}

########### Reed swamp ###########
PredictionBi$Phragmites <- Data$Phragmites
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Phragmites", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Phragmites ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Phragmites"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Phragmites==1)],collapse=";")
}

a  
a$habchar <- "Phragmites"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]
#does model improve from adding Phragmites sequences?
Data$PhragSeq <- colSums(OTUplanter.frequencies2[which(OTUplanter.frequencies2$genus=="Phragmites"),4:133])
Data$PhragSeq 

pred.vector <- c("insect1","insect2","PhragSeq")
best.qda <- qda(as.formula(paste("Phragmites ~ ", paste(pred.vector, collapse = "+"))), data = Data, CV = TRUE)
best.cr <- table(best.qda$class, Data[,"Phragmites"])
a$best.vars[1] <- paste(pred.vector, collapse=" ")
a$correctly.predicted[1] <- sum(diag(prop.table(best.cr)))
a$true.positives[1] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
a$false.neg.sites[1] <- paste(Data$site_name[which(best.qda$class==0 & Data$Phragmites==1)],collapse=";")

a  
a$habchar <- "PhragSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)] #not improved!

# Predicting using leave-one-out (run model above again)
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Willow[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Phclass[i] <- Pred$class
  PredictionBi$Phpost[i] <- Pred$posterior[2]
}


########### Atlantic ###########
PredictionBi$Atlantic <- Data$Atlantic
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Atlantic", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Atlantic ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Atlantic"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Atlantic==1)],collapse=";")
}

a #no good models.  
a$habchar <- "Atlantic"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]

########### Jutland ###########
PredictionBi$Jutland <- Data$Jutland
# Model selection: 
for (i in seq_along(primers)) {
  pred.vector <- paste(primers[i],c(1:4),sep="")
  best.vars <- stepCV.qda("Jutland", pred.vector, Data)
  best.qda <- qda(as.formula(paste("Jutland ~ ", paste(best.vars, collapse = "+"))), data = Data, CV = TRUE)
  best.cr <- table(best.qda$class, Data[,"Jutland"])
  a$best.vars[i] <- paste(best.vars, collapse=" ")
  a$correctly.predicted[i] <- sum(diag(prop.table(best.cr)))
  a$true.positives[i] <- best.cr[2,2]/(best.cr[2,2]+best.cr[1,2])
  a$false.neg.sites[i] <- paste(Data$site_name[which(best.qda$class==0 & Data$Jutland==1)],collapse=";")
}

a #insect give highest correctly predicted  
a$habchar <- "Jutland"
ResultBi <- rbind(ResultBi,a)
ResultBi[,c(1:3,5)]
# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[2])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Jutland[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$Jyclass[i] <- Pred$class
  PredictionBi$Jypost[i] <- Pred$posterior[2]
}

write.table(PredictionBi, file=here::here("data","HabClassPredictions.txt"))
write.table(format(ResultBi, digits=2), file=here::here("data","HabClassModelResults.txt"))

############# Figures of Ellenberg ###############
library(ggplot2)

nul <- which(Data$stratum=="Field")
en <- which(Data$bio_stratum=="Agricultural")
to <- which(Data$bio_stratum=="EarlyDryPoor")
tre <- which(Data$bio_stratum=="EarlyDryRich")
fire <- which(Data$bio_stratum=="EarlyWetPoor")
fem <- which(Data$bio_stratum=="EarlyWetRich")
seks <- which(Data$bio_stratum=="LateDryPoor")
syv <- which(Data$bio_stratum=="LateDryRich")
otte <- which(Data$bio_stratum=="LateWetPoor")
ni <- which(Data$bio_stratum=="LateWetRich")
Data$new.stratum <- NA
Data$new.stratum[en] <- "Farmland"
Data$new.stratum[to] <- "Dry heath"
Data$new.stratum[tre] <- "Dry grassland"
Data$new.stratum[fire] <- "Mire"
Data$new.stratum[fem] <- "Fen"
Data$new.stratum[seks] <- "Moder forest"
Data$new.stratum[syv] <- "Mull forest"
Data$new.stratum[otte] <- "Forest mire"
Data$new.stratum[ni] <- "Swamp forest"

par(mar=c(5,8,2,2))

pdf(width=5.5, height=4,"figures/Ellenberg_RL.pdf")
a <- ggplot(Data[c(en,to,tre,seks,syv),], aes(ellenberg_r, ellenberg_l, color=new.stratum)) +
  scale_x_continuous(name = "pH") + 
  scale_y_continuous(name = "Light") +
  stat_ellipse(type="norm",size=1.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_blank(), legend.key = element_rect(fill = "white"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        text = element_text(size=16))+
  scale_colour_manual(values = c("Farmland" = "darkgrey", "Dry heath" = "darkorchid", 
                                 "Dry grassland" = "orange", "Moder forest" = "yellowgreen", 
                                 "Mull forest" = "darkgreen"))
print(a)
dev.off()

pdf(width=5.5, height=4,"figures/Ellenberg_RF.pdf")
b <- ggplot(Env[c(to,tre,fire,fem),], aes(ellenberg_r, ellenberg_f, color = new.stratum)) +
  scale_x_continuous(name = "pH") + 
  scale_y_continuous(name = "Moisture") +
  stat_ellipse(type="norm",size=1.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_blank(),legend.key = element_rect(fill = "white"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        text = element_text(size=16))+
  scale_colour_manual(values = c("Dry heath" = "darkorchid", "Dry grassland" = "orange", 
                                 "Mire" = "springgreen", "Fen" = "hotpink"))
print(b)
dev.off()

pdf(width=5.5, height=4,"figures/Ellenberg_FL.pdf")
c <- ggplot(Data[c(to,fire,syv,otte),], aes(ellenberg_f, ellenberg_l, color = new.stratum)) +
  scale_x_continuous(name = "Moisture") + 
  scale_y_continuous(name = "Light") +
  stat_ellipse(type="norm",size=1.5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_blank(),legend.key = element_rect(fill = "white"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        text = element_text(size=16))+
  scale_colour_manual(values = c("Dry heath" = "darkorchid", "Mire" = "springgreen", "Forest mire" = "lightseagreen","Mull forest" = "darkgreen"))
print(c)
dev.off()

Env$new.stratum[nul] <- "Rotational field"
pdf(width=5.5, height=4,"figures/Ellenberg_NL.pdf")
d <- ggplot(Data[c(nul,to,tre,syv,ni),], aes(ellenberg_n, ellenberg_l, color = new.stratum)) +
  scale_x_continuous(name = "Fertility") + 
  scale_y_continuous(name = "Light") +
  stat_ellipse(type="norm",size=1.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title=element_blank(),legend.key = element_rect(fill = "white"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        text = element_text(size=16))+
  scale_colour_manual(values = c("Rotational field" = "red", "Dry heath" = "darkorchid", 
                                 "Dry grassland" = "orange", "Mull forest" = "darkgreen",
                                 "Swamp forest"="royalblue"))
print(d)
dev.off()

# add predicted Ellenberg values and confidenceintervals for a target site: 
Ell <- read.table(here::here("data","PredictionEllenbergValuesGAM07032019.txt"))
Ell$Rmax <- Ell$R + Ell$RSE
Ell$Rmin <- Ell$R - Ell$RSE
Ell$Lmax <- Ell$L + Ell$LSE
Ell$Lmin <- Ell$L - Ell$LSE
Ell$Nmax <- Ell$N + Ell$NSE
Ell$Nmin <- Ell$N - Ell$NSE
Ell$Fmax <- Ell$M + Ell$MSE
Ell$Fmin <- Ell$M - Ell$MSE
colnames(Ell)[c(2,4,6,8)] <- c("ellenberg_n","ellenberg_f", "ellenberg_r", "ellenberg_l")

########### ...with predictions ###########
#sample(c(1:130),20)
samples <- c(67,104,40,70,94,13,84,34,12,114,128,61,124,113,46,81,115,74,45,5)

for (i in 1:length(samples)){
  print(a + annotate("point",size = 3, shape = 19, color = "blue", x = Data$ellenberg_r[samples[i]], y = Data$ellenberg_l[samples[i]]) +
          annotate("segment", x = Ell$Rmin[samples[i]], xend = Ell$Rmax[samples[i]], y = Ell$ellenberg_l[samples[i]], yend = Ell$ellenberg_l[samples[i]],
                   colour = "black") + 
          annotate("segment", x = Ell$ellenberg_r[samples[i]], xend = Ell$ellenberg_r[samples[i]], y = Ell$Lmin[samples[i]], yend = Ell$Lmax[samples[i]],
                   colour = "black"))
 here::here("plots",paste(paste("SpillekortGAM",i,"Ellenberg_RL", sep="_"),"pdf", sep=".")) 
 ggsave(here::here("plots",paste(paste("SpillekortGAM",i,"Ellenberg_RL", sep="_"),"pdf", sep=".")), width=5.5, height=4, device=cairo_pdf)
}

for (i in 1:length(samples)){
  print(b + annotate("point",size = 3, shape = 19, color = "blue", x = Data$ellenberg_r[samples[i]], y = Data$ellenberg_f[samples[i]]) +
          annotate("segment", x = Ell$Rmin[samples[i]], xend = Ell$Rmax[samples[i]], y = Ell$ellenberg_f[samples[i]], yend = Ell$ellenberg_f[samples[i]],
                   colour = "black") + 
          annotate("segment", x = Ell$ellenberg_r[samples[i]], xend = Ell$ellenberg_r[samples[i]], y = Ell$Fmin[samples[i]], yend = Ell$Fmax[samples[i]],
                   colour = "black"))
  ggsave(here::here("plots",paste(paste("SpillekortGAM",i,"Ellenberg_RF", sep="_"),"pdf", sep=".")), width=5.5, height=4, device=cairo_pdf)
  dev.off()
}

for (i in 1:length(samples)){
  print(c + annotate("point",size = 3, shape = 19, color = "blue", x = Env$ellenberg_f[samples[i]], y = Env$ellenberg_l[samples[i]]) +
          annotate("segment", x = Ell$Fmin[samples[i]], xend = Ell$Fmax[samples[i]], y = Ell$ellenberg_l[samples[i]], yend = Ell$ellenberg_l[samples[i]],
                   colour = "black") + 
          annotate("segment", x = Ell$ellenberg_f[samples[i]], xend = Ell$ellenberg_f[samples[i]], y = Ell$Lmin[samples[i]], yend = Ell$Lmax[samples[i]],
                   colour = "black"))
  ggsave(here::here("plots",paste(paste("SpillekortGAM",i,"Ellenberg_FL", sep="_"),"pdf", sep=".")), width=5.5, height=4, device=cairo_pdf)
  dev.off()
}
Data$new.stratum[nul] <- "Rotational field"
for (i in 1:length(samples)){
  print(d + annotate("point",size = 3, shape = 19, color = "blue", x = Env$ellenberg_n[samples[i]], y = Env$ellenberg_l[samples[i]]) +
          annotate("segment", x = Ell$Nmin[samples[i]], xend = Ell$Nmax[samples[i]], y = Ell$ellenberg_l[samples[i]], yend = Ell$ellenberg_l[samples[i]],
                   colour = "black") + 
          annotate("segment", x = Ell$ellenberg_n[samples[i]], xend = Ell$ellenberg_n[samples[i]], y = Ell$Lmin[samples[i]], yend = Ell$Lmax[samples[i]],
                   colour = "black"))
  ggsave(here::here("plots",paste(paste("SpillekortGAM",i,"Ellenberg_NL", sep="_"),"pdf", sep=".")), width=5.5, height=4, device=cairo_pdf)
  dev.off()
}

############ Figures of QDA results #################
# import predictions from QDA
Prediction <- read.table(here::here("data","HabClassPredictions.txt"))

#Highforest
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,3],x = Prediction[,3], y = Prediction[,5])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "HighForest", breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],5]), color="red"))
  ggsave(here::here("data",paste(paste("BinaryPred",i,"Highforest", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Forest
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,6],x = Prediction[,6], y = Prediction[,8])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Forest",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],8]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Forest", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Agriculture
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,9],x = Prediction[,9], y = Prediction[,11])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Agriculture",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],11]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Agriculture", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Beech
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,13],x = Prediction[,13], y = Prediction[,15])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Beech",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],15]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Beech", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Oak
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,16],x = Prediction[,16], y = Prediction[,18])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Oak",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],18]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Oak", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Willow
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,19],x = Prediction[,19], y = Prediction[,21])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Willow",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],21]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Willow", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Dwarfshrubs
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,23],x = Prediction[,23], y = Prediction[,25])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Dwarfshrubs",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],25]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Dwarfshrubs", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Alder
names(Prediction)
Prediction$Alder <- Data$Alnus
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,26],x = Prediction[,26], y = Prediction[,28])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Alder",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],28]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Alder", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Phragmites
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,29],x = Prediction[,29], y = Prediction[,31])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Reed",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],31]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Reed", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#Jutland
names(Prediction)
for (i in 1:length(samples)){ 
  print(ggplot(Prediction, aes(group=Prediction[,33],x = Prediction[,33], y = Prediction[,35])) +  
          geom_boxplot(fill = "grey", colour = "black") +
          scale_x_discrete(name = "Jutland",breaks = seq(0, 1,1),limits=c(0, 1)) + 
          scale_y_continuous(name = "Probability", breaks = seq(0, 1,1), limits=c(0,1)) +
          theme_bw() + theme(text=element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                             legend.title=element_blank(),legend.key = element_rect(fill = "white"),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"))+
          geom_hline(aes(yintercept=Prediction[samples[i],35]), color="red"))
  ggsave(here::here("plots",paste(paste("BinaryPred",i,"Jutland", sep="_"),"pdf", sep=".")), width=4, height=4, device=cairo_pdf)
  dev.off()
}

#### END ####