setwd("C:\\Users\\12758\\Desktop\\GMOO-code\\virtual space construction")
##############define a data frame
vir.data=matrix(,1,13)
colnames(vir.data)=c("Ti","Ni","Cu","Hf","Zr","Nb","Co","Cr","Fe","Mn","Pd","V","Al")
#Elements, limits and steps
#Ni–Ti–Cu-M(N) shape memory alloys 
Co=seq(0,6,0.4)#
Fe=seq(0,2,0.4)#
Cr=seq(0,2,0.2)#
Cu=seq(5,15,0.4)#
Hf=seq(0,15,0.4)#
V=seq(0,6,0.4)#
Al=seq(0,4.8,0.4)#
Nb=seq(0,4.8,0.4)#
Zr=seq(0,4.8,0.4)#

######## add elements: Hf,Cu,Cr,Fe,Co,zr
temp = data.frame(Co = rep(Co, each = length(Zr)), Zr = rep(Zr, length(Co)))
temp1 = do.call("rbind", replicate(length(Fe), temp, simplify = F))
temp = data.frame(Fe = rep(Fe, each = dim(temp)[1]), temp1)
rm(temp1)
temp1 = do.call("rbind", replicate(length(Cr), temp, simplify = F))
temp = data.frame(Cr = rep(Cr, each = dim(temp)[1]), temp1)
rm(temp1)
temp1 = do.call("rbind", replicate(length(Cu), temp, simplify = F))
temp = data.frame(Cu = rep(Cu, each = dim(temp)[1]), temp1)
rm(temp1)
temp1 = do.call("rbind", replicate(length(Hf), temp, simplify = F))
temp = data.frame(Hf = rep(Hf, each = dim(temp)[1]), temp1)
rm(temp1)

f<-function(x) sum(x==0)
aa=as.data.frame(apply(temp,1,f))
colnames(aa)="order"
data.all1=cbind(temp,aa)
data.all1=data.all1[-which(data.all1[,"order"]==0),]######## Delete rows that contain 6 elements doped in the alloy
data.all1=data.all1[-which(data.all1[,"order"]==1),]######## Delete rows that contain 5 elements doped in the alloy
data.all1=data.all1[-which(data.all1[,"order"]==2),]######## Delete rows that contain 4 elements doped in the alloy
temp=data.all1[,-ncol(data.all1)]

# add element V
temp1 = do.call("rbind", replicate(length(V), temp, simplify = F))
temp = data.frame(V = rep(V, each = dim(temp)[1]), temp1)
rm(temp1)
aa=as.data.frame(apply(temp,1,f))
colnames(aa)="order"
data.all1=cbind(temp,aa)
data.all1=data.all1[-which(data.all1[,"order"]==3),]######## Delete rows that contain 4 elements doped in the alloy
temp=data.all1[,-ncol(data.all1)]


# add element Al
temp1 = do.call("rbind", replicate(length(Al), temp, simplify = F))
temp = data.frame(Al = rep(Al, each = dim(temp)[1]), temp1)
rm(temp1)
aa=as.data.frame(apply(temp,1,f))
colnames(aa)="order"
data.all1=cbind(temp,aa)
data.all1=data.all1[-which(data.all1[,"order"]==4),]######## Delete rows that contain 4 elements doped in the alloy
temp=data.all1[,-ncol(data.all1)]

# add element Nb
temp1 = do.call("rbind", replicate(length(Nb), temp, simplify = F))
temp = data.frame(Nb = rep(Nb, each = dim(temp)[1]), temp1)
rm(temp1)
aa=as.data.frame(apply(temp,1,f))
colnames(aa)="order"
data.all1=cbind(temp,aa)
data.all1=data.all1[-which(data.all1[,"order"]==5),]######## Delete rows that contain 4 elements doped in the alloy
data.all1=data.all1[-which(data.all1[,"order"]==8),]######## Delete rows that only contain cu element doped in the alloy
temp=data.all1[,-ncol(data.all1)]

temp11=temp
temp0=temp

temp0=matrix(0,nrow=nrow(temp),ncol=13)
colnames(temp0)=c("Ti","Ni","Cu","Hf","Zr","Nb","Co","Cr","Fe","Mn","Pd","V","Al")
temp0[,"Hf"]=temp[,"Hf"]
temp0[,"Cu"]=temp[,"Cu"]
temp0[,"Zr"]=temp[,"Zr"]
temp0[,"Nb"]=temp[,"Nb"]
temp0[,"Co"]=temp[,"Co"]
temp0[,"Cr"]=temp[,"Cr"]
temp0[,"Fe"]=temp[,"Fe"]
temp0[,"Mn"]=0
temp0[,"Pd"]=0
temp0[,"V"]=temp[,"V"]
temp0[,"Al"]=temp[,"Al"]

temp=temp0

Ti.site=temp[,"Hf"]+0.5*temp[,"Nb"]+temp[,"Zr"]+temp[,"Al"] ###############Ti-site doping
Ni.site=temp[,"Cu"]+0.5*temp[,"Nb"]+temp[,"Co"]+temp[,"Cr"]+temp[,"Fe"]+temp[,"V"]###############Ni-site doping

########Ni-site
seq_j=seq(47,50.6,0.4)

temp[,"Ti"]=100-seq_j[1]-Ti.site #Ti
temp[,"Ni"]=seq_j[1]-Ni.site #Ni
#constrain Ti,Ni content
data=temp[which(temp[,"Ti"]>=20),]
data.all=data[which(data[,"Ni"]>=35),]

for(j in 2:length(seq_j)){
  temp0[,"Ti"]=100-seq_j[j]-Ti.site
  temp0[,"Ni"]=seq_j[j]-Ni.site
  data=temp0[which(temp0[,"Ti"]>=20),]
  data2=data[which(data[,"Ni"]>=35),]
  data.all=rbind(data.all,data2)
}
data.all=as.data.frame(data.all)



f<-function(x) sum(x==0)
aa=as.data.frame(apply(data.all,1,f))
colnames(aa)="order"
data.all1=cbind(data.all,aa)
bb=as.data.frame(apply(data.all,1,sum))

data.all.u=unique(data.all)
data.all.u=as.data.frame(data.all.u)
write.csv(data.all.u,"vir.data0.4.csv")

