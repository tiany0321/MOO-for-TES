
##################library
library(GPareto)
library(DiceDesign)
library(DiceKriging)
library(dplyr)
library(RColorBrewer)
display.brewer.all(type = "seq")

#####################plot function
fn_plot=function(selector,test.response.grid,response.grid,model.result,ref00,tar_point00,train_pareto,pareto_P,choose_p){
  
    x_lim=c(-38,2)
    y_lim=c(0,105)
    plot(test.response.grid[,"es1"],test.response.grid[,"es2"],cex=1.5, col=brewer.pal(9, "YlOrRd")[3],pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    par(new=T)
    plot(response.grid[,"es1"],response.grid[,"es2"], cex=1.5,col=brewer.pal(9, "Greens")[5],pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    #####reference point
    par(new=T)
    points(ref00[1],ref00[2],cex=1.5,pch=15,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    ######ideal point
    par(new=T)
    points(tar_point00[1],tar_point00[2],cex=1.5,pch=17,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    ###########pareto for training data
    par(new=T)
    plot(train_pareto[,"es1"],train_pareto[,"es2"],cex=1.5, col=brewer.pal(9, "YlOrRd")[6],type="p",pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    par(new=T)
    plot(train_pareto[,"es1"],train_pareto[,"es2"],cex=1.5,type="s",pch=16,col="red",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    ###########current best trade-off options
    par(new=T)
    plot(pareto_P[,"es1"],pareto_P[,"es2"],cex=1.5, col="cyan",type="p",pch=16,xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    par(new=T)
    plot(pareto_P[,"es1"],pareto_P[,"es2"],cex=1.5,type="s",pch=16,col="cyan",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    ##########constrain
    par(new=T)
    points(test.response.grid[choose_p,1],test.response.grid[choose_p,2],cex=1.5,pch=17,col="blue",xlab="Objective1",ylab="Objective2",ylim=y_lim,xlim=x_lim, lwd=3, tck=0.01,font=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
    
    text(x = tar_point00[1]+(ref00[1]-tar_point00[1])*0.05, y = ref00[2]-(ref00[2]-tar_point00[2])*0.05, labels = paste(selector,as.character(iter),sep="-iter"),cex = 1,font=1.7)
    par(new=T)
    arrows(ref00[1], ref00[2], tar_point00[1], tar_point00[2], col = "red")
}

#######pareto front function
pareto.front=function(data=response, obj){
    response_obj=as.data.frame(data[,obj])
  
    response_obj=response_obj[order(response_obj$es1,response_obj$es2,decreasing=FALSE),]
    num=which(!duplicated(cummin(response_obj$es2)))
    front = response_obj[num,]
    front_num=as.numeric(rownames(front))
  
  return(list(front=front,front_num=front_num))
}

###########decide the coordinates of the options which will construct the constrains
fn_pareto_P <- function(crit=NULL,vir_target,train_data,train_pareto,ref00,tar_point00,scale){

  if(crit=="3points"){
    ###########decide the current best trade-off option by maximizing cos_theta
    cos_theta=cos_theta_calculate(vir_point0=train_pareto,ref00,tar_point00,scale)
    opt_data=train_pareto[which(cos_theta[]==max(cos_theta[])),]
    ###########range the optimal options in training data according to cos_theta
    train_pareto_th=cbind(train_pareto,cos_theta)
    train_pareto_th_order=train_pareto_th[order(-train_pareto_th[,"cos_theta"]),]
    ###########choose three current best trade-off options
    point2=train_pareto_th_order[which(train_pareto_th_order[,"cos_theta"]==max(train_pareto_th_order[,"cos_theta"])),]
    train_pareto_th_order$abs=train_pareto_th_order[,"es1"]-point2[,"es1"]
    point3_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]<0),]
    point3=point3_c[1,]
    point1_c=train_pareto_th_order[which(train_pareto_th_order[,"abs"]>0),]
    point1=point1_c[1,]
    
    keep=c("es1","es2")
    pareto_P0=t(as.data.frame(c(point1[,"es1"],-Inf)))
    pareto_P1=point1[,keep]
    pareto_P2=point2[,keep]
    pareto_P3=point3[,keep]
    pareto_P4=t(as.data.frame(c(-Inf,point3[,"es2"])))
    
    colnames(pareto_P0)=keep
    colnames(pareto_P4)=keep
    ###########return the coordinates of the options which will construct the constrains
    pareto_P=as.data.frame(rbind(pareto_P0,pareto_P1,pareto_P2,pareto_P3,pareto_P4))
  }
  print(pareto_P)
  return(pareto_P)
}

###########cos_theta calculation function
cos_theta_calculate=function(vir_point0,ref00,tar_point00,scale){

  ref0=as.data.frame(t(as.data.frame(ref00)))
  colnames(ref0)=colnames(vir_point0)
  tar_point0=as.data.frame(t(as.data.frame(tar_point00)))
  colnames(tar_point0)=colnames(vir_point0)
  
  all=rbind(ref0,tar_point0,vir_point0)
  
  #######scale
  if(scale==TRUE){
    all[,1]=scale(all[,1])
    all[,2]=scale(all[,2])
  }
  ref=unlist(all[1,])
  tar_point=unlist(all[2,])
  vir_point=all[-c(1,2),]
  
  a1=(ref[1]-vir_point[,1])*(ref[1]-tar_point[1])+(ref[2]-vir_point[,2])*(ref[2]-tar_point[2])
  a2=(((ref[1]-vir_point[,1])^2+(ref[2]-vir_point[,2])^2)*((ref[1]-tar_point[1])^2+(ref[2]-tar_point[2])^2))^0.5
  cos_theta=a1/a2
  # angle=as.data.frame(EHI_grid*cos_theta)
  return(cos_theta)
}


####################psi function
fn_psi_ab=function(a_val,b_val,mean_val,sd_val){
  nor_val = (b_val-mean_val )/sd_val
  z = nor_val
  psi_ab= sd_val*dnorm(z)+(a_val-mean_val)*pnorm(z)   
  return(psi_ab)
}

###############Expected Hypervolume Improvement(EHVI)##########################################################

#########Expected Hypervolume Improvement(EHVI)
fn_EHVI=function(ref_poi=ref00,pareto_Data=pareto_Data,model.result){
  
  pareto_Data=pareto_Data[order(-pareto_Data[,"es1"]),]
  y_0=t(as.data.frame(c(ref_poi[1],-Inf)))
  colnames(y_0)=c("es1","es2")
  y_i1=t(as.data.frame(c(-Inf,ref_poi[2])))
  colnames(y_i1)=c("es1","es2")   
  pareto_Data=rbind(y_0,pareto_Data,y_i1)
  i_max=nrow(pareto_Data)-1
  EHVI=rep()
  
  for (i in 1:i_max){
    z_1=(pareto_Data[i+1,"es1"]-model.result["mean1"] )/model.result["sd1"]
    iterm1=(pareto_Data[i,"es1"]-pareto_Data[i+1,"es1"])*pnorm(z_1)*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    iterm2=(fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i,"es1"],model.result["mean1"],model.result["sd1"])-fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i+1,"es1"],model.result["mean1"],model.result["sd1"]))*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    EHVI[i]=iterm1+iterm2
  }
  EHVI=as.data.frame(EHVI)
  EHVI=na.omit(EHVI)
  EHVI_sum=sum(EHVI)
  return(EHVI_sum)
}
#########EHVI evaluation function
fn_EHI=function(test.grid,train_pareto,ref_poi=ref00,tar_point00,scale,opt_algri=fn_EHVI,model){
  #km prediction
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  
  mu <- as.data.frame(t(as.data.frame(pred$mean)))
  colnames(mu)=c("mean1","mean2")
  
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  colnames(sigma)=c("sd1","sd2")
  ##################new addition 
  model.result=cbind(mu,sigma)
  #########EHVI evaluation
  EHI_result=apply(model.result,1,opt_algri,ref_poi=ref_poi,pareto_Data=train_pareto)
  EHI_result=as.data.frame(EHI_result)
  EHI=cbind(model.result,EHI_result)
  colnames(EHI)=c("mean1","mean2","sd1","sd2","score")
  return(EHI)
}

################################################################################################################


####################################constrained EHVI###########################################################
####################constrained EHVI
fn_P_EHVI=function(ref_poi=ref00,pareto_Data=pareto_Data,model.result){
  
  pareto_Data=pareto_Data[order(-pareto_Data[,"es1"]),]
  pareto_Data=pareto_Data[c(1,3,5),]
  colnames(pareto_Data)=c("es1","es2")  
  i_max=nrow(pareto_Data)-1
  EHVI=rep()
  
  for (i in 1:i_max){
    z_1=(pareto_Data[i+1,"es1"]-model.result["mean1"] )/model.result["sd1"]
    iterm1=(pareto_Data[i,"es1"]-pareto_Data[i+1,"es1"])*pnorm(z_1)*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    iterm2=(fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i,"es1"],model.result["mean1"],model.result["sd1"])-fn_psi_ab(pareto_Data[i,"es1"],pareto_Data[i+1,"es1"],model.result["mean1"],model.result["sd1"]))*fn_psi_ab(pareto_Data[i+1,"es2"],pareto_Data[i+1,"es2"],model.result["mean2"],model.result["sd2"])
    EHVI[i]=iterm1+iterm2
  }
  EHVI=as.data.frame(EHVI)
  EHVI=na.omit(EHVI)
  EHVI_sum=sum(EHVI)
  return(EHVI_sum)
}

#########constrained EHVI evaluation function
fn_3points_EHI=function(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_P_EHVI,model){
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
  mu <- as.data.frame(t(as.data.frame(pred$mean))) 
  colnames(mu)=c("mean1","mean2")
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))
  colnames(sigma )=c("sd1","sd2")
  
  ##################new addition  
  model.result=cbind(mu,sigma)
  #the coordinates of the options which will construct the constrains
  pareto_P=fn_pareto_P(crit,vir_target,train_data,train_pareto,ref00,tar_point00,scale)
  
  EHI_result=apply(model.result,1,opt_algri,ref_poi=ref00,pareto_Data=pareto_P)
  EHI_result=as.data.frame(EHI_result)
  EHI=cbind(model.result,EHI_result)
  colnames(EHI)=c("mean1","mean2","sd1","sd2","score")
  return(EHI)
}
################################################################################################################


###############minimize the distance to the ideal point##########################################################

#########distance calculation function
fn_dist=function(model.result.option,tar_point00){
  dist_cal=((model.result.option["mean1"]-tar_point00[1])^2+(model.result.option["mean2"]-tar_point00[2])^2)^0.5
  return(dist_cal)
}

#########mean_dist evaluation function
fn_mean_dist=function(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_dist,model){
  pred <- predict_kms(model, newdata=test.grid, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
 
   mu <- as.data.frame(t(as.data.frame(pred$mean)))
  sigma <- as.data.frame(t(as.data.frame(pred$sd)))

  colnames(mu)=c("mean1","mean2")
  colnames(sigma )=c("sd1","sd2")
  ##################new addition  
  model.result=cbind(mu,sigma)
  
  dist_result=unlist(apply(model.result,1,opt_algri,tar_point00))
  dist_result=as.data.frame(-dist_result)
  dist=cbind(model.result,dist_result)
  colnames(dist)=c("mean1","mean2","sd1","sd2","score")
  return(dist)
}

################################################################################################################




############main function for bi-objective optimization #####################################################

fn_multi_obj=function(crit=NULL,vir_target,train_data,selector_num,test.grid,train_pareto,ref00,tar_point00,scale,model){
  ###### EHVI
  if(selector_num=="EHVI"){
    opt_algri=fn_EHVI
    score_calculate=fn_EHI(test.grid,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  ######constrained EHVI
  if(selector_num=="3points_EHI"){
    opt_algri=fn_P_EHVI
    score_calculate= fn_3points_EHI(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri,model)
  }
  #####minimize the distance to the ideal point
  if(selector_num=="mean_dist"){
    opt_algri=fn_dist
    score_calculate=fn_mean_dist(crit,vir_target,test.grid,train_data,train_pareto,ref00,tar_point00,scale,opt_algri=fn_dist,model)
  }
  return(score_calculate)
}

#######################################################################################################