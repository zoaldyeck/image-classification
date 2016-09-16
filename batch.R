library(EBImage)
library(ForeCA)
library(amap)
library(e1071)

testpath1<-"/protestify_data/images/addlexpts/poster"
testpath2<-"/protestify_data/images/originals"
offset<-8000
num<-50
label1<-+1
label2<--1
accuarcy1<-mymodel(testpath1,model,kcentroids,offset,num,label1)
accuarcy2<-mymodel(testpath2,model,kcentroids,offset,num,label2)

rm(testpath1)
rm(testpath2)
rm(offset)
rm(num)
rm(accuarcy1)
rm(accuarcy2)