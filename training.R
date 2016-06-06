library(EBImage)
library(ForeCA)
library(amap)
library(e1071)

rm(list=ls())
source(file="/home/visionexpts/newsource3.R")

k<-500
unit<-16
s<-4
pNum<-150000
all<-7500

poster_list<-list()
origin_list<-list()
poster_train_list<-list()
poster_test_list<-list()
origin_train_list<-list()
origin_test_list<-list()
train_list<-list()

setwd("/protestify_data/images/addlexpts/poster")
poster_jpg <- list.files(pattern="*.jpg")
poster_list <- GrayResize(poster_jpg)
rm(poster_jpg)
setwd("/protestify_data/images/originals")
origin_jpg <- list.files(pattern = "*.jpg")
origin_list <- GrayResize(origin_jpg)
setwd("/home/visionexpts")
rm(origin_jpg)
poster_num <- length(poster_list)
origin_num <- length(origin_list)
poster_test_index <- sample(1:poster_num, round(0.20*poster_num))
origin_test_index <- sample(1:origin_num, round(0.20*origin_num))
rm(poster_num)
rm(origin_num)
poster_test_list <- poster_list[poster_test_index]       #108
poster_train_list <- poster_list[-poster_test_index]     #432
rm(poster_list)
rm(poster_test_index)
origin_test_list <- origin_list[origin_test_index]       #74
origin_train_list <- origin_list[-origin_test_index]     #295
rm(origin_list)
rm(origin_test_index)
train_list <- c(poster_train_list,origin_train_list)  #727

patch_Mat <- patch_extraction(train_list,pNum,unit)
rm(train_list)
normMat <- pre_process(patch_Mat)
rm(patch_Mat)
matZCA <- ZCAwhite(normMat,0.01)


kcenCo1<-mykmean(t(matZCA),k,"correlation",15)

kcentroids <- t(kcenCo1)

rm(matZCA)

k<-ncol(kcentroids)
poster_train__feature_mat <- feature_Generator(poster_train_list,kcentroids)
origin_train_feature_mat <- feature_Generator(origin_train_list,kcentroids)
rm(poster_train_list)
rm(origin_train_list)
poster_test__feature_mat <- feature_Generator(poster_test_list,kcentroids)
origin_test_feature_mat <- feature_Generator(origin_test_list,kcentroids)
rm(poster_test_list)
rm(origin_test_list)


poster_train <- rbind(poster_train__feature_mat,+1)
origin_train <- rbind(origin_train_feature_mat,-1)
poster_test <- rbind(poster_test__feature_mat,+1)
origin_test <- rbind(origin_test_feature_mat,-1)
train <- cbind(poster_train,origin_train)
test <- cbind(poster_test,origin_test)
rm(poster_train)
rm(origin_train)
rm(poster_test)
rm(origin_test)
rm(poster_train__feature_mat)
rm(origin_train_feature_mat)
rm(poster_test__feature_mat)
rm(origin_test_feature_mat)

train_x<-train[-nrow(train),]
train_y<-train[nrow(train),]
train_x<-t(train_x)

test_x<-test[-nrow(test),]
test_y<-test[nrow(test),]
test_x<-t(test_x)

model<-e1071::svm(train_x,train_y,type="C-classification",kernel="linear")
test_pred<-predict(model,test_x)
sum(test_pred==test_y)/length(test_pred)


rm(train)
rm(test)
rm(train_x)
rm(train_y)
rm(test_x)
rm(test_y)