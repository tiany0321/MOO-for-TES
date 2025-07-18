setwd("C:\\Users\\12758\\Desktop\\GMOO-code\\prediction\\GHVI-iter0")
set.seed(243)
source("function.R")
t1=proc.time()
#############preset ideal point
tar_point00=c(-35,10)

#############preset reference point
ref00=c(0,100)

crit="3points"
obj=c("es1","es2")
obj_opt="both_min"
scale=T
selector="3points_EHI"#######################selector

#############key features for latent heat prediction and hysteresis prediction
keep.heat=c("Lambda","VEN","delta.x","ou")
keep.hysteresis=c("delta.r","VEN","Lambda")
keep=c("Lambda","delta.x","VEN","ou","delta.r")
objective1.keep=keep.heat
objective2.keep=keep.hysteresis

#########read training data
data=read.csv("dataall_train.csv")
col_0=colnames(data)
#########x
design.grid0 <-data[,keep]
colnames(design.grid0)=keep
#########y
response.grid0<-data[,c("enthalpy.lin","hysteresis")]
colnames(response.grid0)=obj

#########read virtual space
virtualdata=read.csv("dataall_vir.csv")
virtualdata0=virtualdata[,]
vir.grid_org <- virtualdata0[,keep]
colnames(vir.grid_org)=keep
vir.grid= anti_join(vir.grid_org,design.grid0)

################model construction
#######latent heat model
formula1<- es1~ Lambda + VEN + delta.x + ou
noiselevel_heat=0.07*(max(response.grid0[,"es1"])-min(response.grid0[,"es1"]))
noise_heat=rep(noiselevel_heat^2,nrow(response.grid0))
set.seed(243)
mf1_0 <- km(formula1, design = design.grid0, response = response.grid0[,1],noise.var=noise_heat,covtype = "exp")
#######hysteresis model
formula2<- es2~ delta.r + VEN + Lambda
noiselevel_hysteresis=0.03*(max(response.grid0[,"es2"])-min(response.grid0[,"es2"]))
noise_hysteresis=rep(noiselevel_hysteresis^2,nrow(response.grid0))
set.seed(243)
mf2_0 <- km(formula2, design = design.grid0, response = response.grid0[,2],noise.var=noise_hysteresis,covtype = "exp")
########combine two model
model0 = list(mf1_0, mf2_0)

#############the pareto frontier for the current training data
train_pareto0=pareto.front(data=response.grid0, obj)$front
#the coordinates of the options which will construct the constrains
pareto_P0=fn_pareto_P(crit,vir_target,train_data=response.grid0,train_pareto0,ref00,tar_point00,scale)

############The hypervolume between the Pareto front of the existing training set and the reference point
hv <- dominated_hypervolume(points=t(train_pareto0), ref=ref00)

##############################Prediction and Evaluation of Virtual Spaces using Multi-Objective Optimization Algorithms
vir_target=tar_point00
multi_obj_calculate0=fn_multi_obj(crit,vir_target,train_data=response.grid0,selector_num=selector,vir.grid,train_pareto0,ref00,tar_point00,scale, model0)
multi_obj_score0=multi_obj_calculate0$score
model.result0=multi_obj_calculate0[,c("mean1","mean2")]
model.result0.sd=multi_obj_calculate0[,c("sd1","sd2")]
vir.element=virtualdata[,c(2:14)]

##########recommended point
multi_obj_num0=as.data.frame(multi_obj_score0)
vir.result=cbind(vir.element,multi_obj_calculate0,vir.grid)

multi_obj_order0=rep(1:nrow(multi_obj_num0))
multi_obj_num0=cbind(multi_obj_num0,multi_obj_order0)
choose0=multi_obj_num0[order(-multi_obj_num0[,1]),]
choose_p0=choose0[1:4,"multi_obj_order0"]
iter=0
test.response.grid=model.result0
colnames(test.response.grid)=obj
iter=0
fn_plot(selector,test.response.grid,response.grid=response.grid0,model.result=model.result0,ref00,tar_point00,train_pareto=train_pareto0,pareto_P=pareto_P0,choose_p=choose_p0)
vir.result[choose_p0,]
write.csv(vir.result,"vir.result-3point.csv")

