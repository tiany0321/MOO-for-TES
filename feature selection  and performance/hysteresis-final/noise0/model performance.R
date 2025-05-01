setwd("C:\\Users\\12758\\Desktop\\GMOO-code\\feature selection  and performance\\hysteresis-final\\noise0")
library("DiceKriging")

phi=0#########noise level prefactor

#######Kriging model function
fn.gp = function(training.data.x, training.data.y,virtual.data.x,phi)
{ set.seed(243)
  noiselevel=phi*(max(data.y0)-min(data.y0))
  noise=rep(noiselevel^2,nrow(training.data.x))
  m<-  km(formula<-~.,design=training.data.x, response=training.data.y , noise.var=noise ,covtype = "exp")
  set.seed(243)
  p <- predict(m, virtual.data.x, "UK")
  sd=as.data.frame(p$sd)
  mean=as.data.frame(p$mean)
  pr=cbind(mean,sd)
  colnames(pr)=c("mean","sd")
  return(pr)
}
#######Continuous Ranked Probability Score
CRPS=function(mean,sd,y){
  zz=(y-mean)/sd
  CRPS=sd*(1/(pi)^0.5-2*dnorm(zz)-zz*(2*pnorm(zz)-1))
  
  return(CRPS)
}
#####plot function
fn.plot.prevsmea = function(pre, meas){
  data = cbind(pre,meas)
  data = data.frame(data)
  colnames(data) = c("pre", "meas")
  if(min(data$meas) > min(data$pre)){mmin = min(data$pre)}else{mmin = min(data$meas)}
  if(max(data$meas) < max(data$pre)){mmax = max(data$pre)}else{mmax = max(data$meas)}
  plot(data$pre ~ data$meas ,tck=0.01,font=2,font.lab=2,
       ylab = "Predicted values", xlab = "Measured values",
       pch = 21, col = "darkblue", cex = 1.5, bg="gray",
       xlim= c(mmin - 5, mmax + 5), ylim = c(mmin - 5, mmax + 5))
  abline(0, 1, lwd =2, col = "red")
}
#####read data
data.training=read.csv("dataall_train.csv")
data=na.omit(data.training)

#############x
keep=c("delta.r","VEN","Lambda")
data0=data[,keep]
############y
data.y0=as.data.frame(data[,"hysteresis"])

############Define data frame for predicted results
heat.train.gauss=matrix(,nrow(data),3)
for(i in 1:nrow(data))
{
  #############leave one out method
  data.x=data0[-i,]#############x_-i as x of training data 
  data.vir.x=data0[i,]##########x_i as test data
  data.y=data.y0[-i,]#############dy_-i as y of training data 
  model.gp<-fn.gp(data.x,data.y,data.vir.x,phi)
  #######mse
  data.vir.y=data.y0[i,]#############y_i as y of test data 
  heat.train.gauss[i,1]=data.vir.y############y_i
  heat.train.gauss[i,2]=model.gp[,1]############prediction value of y_i
  heat.train.gauss[i,3]=model.gp[,2]############the uncertainty of the prediction
}

heat.train.gauss=as.data.frame(heat.train.gauss)
colnames(heat.train.gauss)=c("measured","predicted","sd")
heat.train.gauss=na.omit(heat.train.gauss)
write.csv(heat.train.gauss,"heat.train.gauss.noise0.csv")

data.error.frame=matrix(,1,3)
colnames(data.error.frame)=c("MSE","MAPE","CRPS_v")
#############calculate MSE
heat.train.gauss$abs_delta=abs(heat.train.gauss[,"predicted"]-heat.train.gauss[,"measured"])
data.error.frame[,"MSE"]=sum(heat.train.gauss$abs_delta^2)/nrow(heat.train.gauss)
#############calculate MAPE
heat.train.gauss$delta_ratio=abs((heat.train.gauss[,"predicted"]-heat.train.gauss[,"measured"])/heat.train.gauss[,"measured"])
data.error.frame[,"MAPE"]=sum(heat.train.gauss$delta_ratio)/nrow(heat.train.gauss)
#############calculate CRPS
heat.train.gauss$CRPS=CRPS(heat.train.gauss[,"predicted"],heat.train.gauss[,"sd"],heat.train.gauss[,"measured"])
CRPS_v=mean(abs(heat.train.gauss$CRPS))
CRPS_v
data.error.frame[,"CRPS_v"]=CRPS_v
write.csv(data.error.frame,"data.error.frame0.csv")


pdf("noise0.pdf")
fn.plot.prevsmea(heat.train.gauss[,2],heat.train.gauss[,1])
dev.off()


