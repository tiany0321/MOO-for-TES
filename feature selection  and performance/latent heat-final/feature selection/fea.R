setwd("C:\\Users\\12758\\Desktop\\GMOO-code\\feature selection  and performance\\latent heat-final\\feature selection")

#source package
library(minerva)
library(corrplot)
library(gbm)

###########read training data
data<-read.csv("dataall_train.csv")
data0=data[,-c(1:14)]
#data.hysteresis=data0[,-2]
data.heat=data0[,-1]


#########Gradient boosting
set.seed(484)
boost.heat = gbm(enthalpy.lin ~ ., 
               data = data.heat, distribution = "gaussian", n.trees = 1000, interaction.depth = 7)
rel.inf.heat = summary(boost.heat)
inf.heat=rel.inf.heat[1:10,]
inf.heat
write.csv(inf.heat,"GB.result.heat.csv")

feature.heat=as.character(inf.heat[,1])
feature.names.keep.heat=data.heat[,feature.heat]
aa=as.data.frame(data.heat$enthalpy.lin)
colnames(aa)="enthalpy.lin"
feature.names.keep.heat=cbind(aa,feature.names.keep.heat)

###########Maximal Information Coefficient(MIC)
mic<-mine(x=feature.names.keep.heat,var.thr=1e-10,use="pairwise.complete.obs")$MIC

corrplot(mic)
corrplot(mic, order = 'AOE',type =  "upper")

############The final features for latent heat prediction
feature.keep=c("Lambda","VEN","delta.x","ou")
