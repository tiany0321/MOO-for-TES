setwd("C:\\Users\\12758\\Desktop\\GMOO-code\\feature selection  and performance\\hysteresis-final\\feature selection")
#source package
library(minerva)
library(corrplot)
library("gbm")

###########read training data
data<-read.csv("dataall_train.csv")
data0=data[,-c(1:14)]
data.hysteresis=data0[,-2]
#data.heat=data0[,-1]


#########Gradient boosting
set.seed(44)
boost.hysteresis = gbm(hysteresis ~ ., 
                 data = data.hysteresis, distribution = "gaussian", n.trees = 1000, interaction.depth = 9)
rel.inf.hysteresis = summary(boost.hysteresis)
inf.hysteresis=rel.inf.hysteresis[1:10,]
inf.hysteresis
write.csv(inf.hysteresis,"GB.result.hysteresis.csv")

feature.hysteresis=as.character(inf.hysteresis[,1])
feature.names.keep.hysteresis=data.hysteresis[,feature.hysteresis]
aa=as.data.frame(data.hysteresis$hysteresis)
colnames(aa)="hysteresis"
feature.names.keep.hysteresis=cbind(aa,feature.names.keep.hysteresis)

###########Maximal Information Coefficient(MIC)
mic<-mine(x=feature.names.keep.hysteresis,var.thr=1e-10,use="pairwise.complete.obs")$MIC

corrplot(mic)
corrplot(mic, order = 'AOE',type =  "upper")

############The final features for hysteresis prediction
keep.hysteresis=c("delta.r","VEN","Lambda")