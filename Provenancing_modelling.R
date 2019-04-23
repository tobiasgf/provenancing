# Analysis SoilTracker: provenancing models
# Manuscript: Predicting provenance of forensic soil samples: 
# Linking soil to ecological habitats by metabarcoding and supervised classification
# Author: Camilla Fløjgaard
# Date: 23-04-2019
library(readxl)
library(here)
library(MASS)
library(mgcv)

#### Data: environmental data, NMS-ordination-axes and habitat characteristics
Data <- read.table(here::here("data","SoilTrackerData.txt")) # NMDS data added from the file "ordination_results.txt"

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
Annotated <- read.table(here::here("data","plant_table_final.txt"), header=TRUE, sep="\t")
#Annotated_v2 <- Annotated[which(Annotated$pident>=99),]
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

a #no true positives over 50%, don't ude this for modelling.  
a$habchar <- "ConiferousSeq"
ResultBi <- rbind(ResultBi,a[1,])
ResultBi[,c(1:3,5)]

# Predicting using leave-one-out:
pred.vector <- paste(a$best.vars[1])
for (i in 1:130){
  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Coniferous[-i])
  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
  PredictionBi$CoSeqclass[i] <- Pred$class
  PredictionBi$CoSeqpost[i] <- Pred$posterior[2]
}

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

a #fung maximizes both corretly.predicted. 
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

a #fung maximizes both corretly.predicted. 
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
ResultBi[,c(1:3,5)] #percentage correctly predicted improved, but true positives fell. 

# Predicting using leave-one-out (run first model selection):
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
#correctly predicted decreases and true positives is only 50%

#Predicting using leave-one-out:
#pred.vector <- paste(a$best.vars[1])
#for (i in 1:130){
#  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Heathland[-i])
#  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
#  PredictionBi$HLclass[i] <- Pred$class  
#  PredictionBi$HLpost[i] <- Pred$posterior[2]
#}

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
Abundance = data.frame(read_excel(here::here("data","BiowideAbundans.xlsx")))
Alnus.sites <- Abundance$plot[which(Abundance$id =="Alnus_glutinosa" & Abundance$abundans ==3)]
Data[Data$Site_nr %in% Alnus.sites,1:5]
# efter vurdering af billeder:
Alnus.sites2 <- Alnus.sites[-c(2,22,66,71,74)]
Data$Alnus <- "0"
Data$Alnus[which(Data$Site_nr %in% Alnus.sites2)] <-"1"

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
Phragmites.sites <- Abundance$plot[which(Abundance$id =="Phragmites_australis" & Abundance$abundans ==3)]
Data[Data$Site_nr %in% Phragmites.sites,1:5]
# efter vurdering af billeder:
Phragmites.sites2 <- Phragmites.sites[-c(37,50,91,92)]
Data$Phragmites <- "0"
Data$Phragmites[which(Data$Site_nr %in% Phragmites.sites2)] <-"1"

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
  PredictionBi$Wiclass[i] <- Pred$class
  PredictionBi$Wipost[i] <- Pred$posterior[2]
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
# Predicting using leave-one-out:
#pred.vector <- paste(b$best.vars[1])
#for (i in 1:130){
#  best.qda <- qda(as.matrix(Data[-i,paste(unlist(strsplit(pred.vector, split=" ")))]), Data$Atlantic[-i])
#  Pred <- predict(best.qda, newdata=Data[i,paste(unlist(strsplit(pred.vector, split=" ")))])
#  PredictionBi$Atclass[i] <- Pred$class
#  PredictionBi$Atpost[i] <- Pred$posterior[2]
#}

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
write.table(format(ResultBi, digits=2), file=here:here("data","HabClassModelResults.txt"))

