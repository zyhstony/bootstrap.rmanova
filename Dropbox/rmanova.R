require("nlme")
require("lme4")
require("robust")
library(foreach)
library(doMC)
registerDoMC(10)
options(contrasts=c('contr.sum','contr.poly'))
########################defined function#################################
resample<-function(data,change=0)
{
  #resample subject while restricting repeated measures within subject 
  residue<-rep(0,nrow(data))
  f<-tapply(data[,1],factor(paste(data$Race,data$Tumor)),mean)
  for(i in 1:nrow(data))
  {
    name<-paste(data$Race[i],data$Tumor[i])
    temp<-data[i,1]-f[name]
    residue[i]<-unlist(temp)
  }
  index<-data$id
  if(change==1)
  {
    #bootstrap residual without restricting on between-subject groups
    index.rs<-sample(1:length(unique(index)),replace=T)
    residue.rs<-residue
    for(i in 1:(length(residue)/2))
    {
      residue.rs[2*i-1]<-residue[index.rs[i]*2-1]
      residue.rs[2*i]<-residue[index.rs[i]*2]
    }
    
  }
  if(change==0)
  {
    #bootstrap residual with restricting on between-subject groups
    
    index.rs<-sample(1:(length(index)/2))
    index.rs[1:(length(index)/4)]<-sample(1:(length(index)/4),replace=T)
    index.rs[1:(length(index)/4)+(length(index)/4)]<-sample((length(index)/4)+1:(length(index)/4),replace=T)
    residue.rs<-residue
    for(i in 1:(length(residue)/2))
    {
      residue.rs[2*i-1]<-residue[index.rs[i]*2-1]
      residue.rs[2*i]<-residue[index.rs[i]*2]
    }
  }
  #print(index.rs)
  data[,1]<-residue.rs
  return(data)
}
quant<-function(value,dis)
{
  #get bootstrap pvalue
  return(sum(dis>value)/length(dis))
}

get.p<-function(data)
{
  #get F value for each bootstrap run
  varlist<-names(data)[1:(ncol(data)-3)]
  colnames(data)[1]<-"miRNA"
  gen.formula<- function(x) {
    ezANOVA(data=data,within=.(Tumor),wid=.(id),between=.(Race),dv=miRNA,type=3)
  }
  get.F<-function(x)
  {
    x$ANOVA[,4]
  }
  
  models <- mclapply(varlist,gen.formula)  
  f<-unlist(lapply(models,get.F))
  p<-array(f,c(3,ncol(data)-3))
  return(p)
}

get.p1<-function(data)
{
  #get p value for each bootstrap run
  varlist<-names(data)[1:(ncol(data)-3)]
  colnames(data)[1]<-"miRNA"
  gen.formula<- function(x) {
  ezANOVA(data=data,within=.(Tumor),wid=.(id),between=.(Race),dv=miRNA,type=3)
  }
  get.F<-function(x)
  {
    x$ANOVA[,5]
  }
  
    models <- mclapply(varlist,gen.formula)  
f<-unlist(lapply(models,get.F))
  p<-array(f,c(3,ncol(data)-3))
  return(p)
}
###########################################data preparition##############################
setwd("~/MicroRNA")
data.all<-read.csv("ANOVA/Genespring_5_27.csv",stringsAsFactors = FALSE,header=T)[,1:123]
rownames(data.all)<-data.all[,1]
colnames(data.all)<-gsub(".gTotalP.*","",colnames(data.all))
data.all<-data.all[,-1]
SBU<-data.all[,1:60]
WU<-data.all[,61:122]
SBU<-t(SBU)
rownames(SBU)<-gsub("US84103569_","",rownames(SBU))
rownames(SBU)<-gsub("S01_miRNA_107_Sep09_","",rownames(SBU))
rownames(SBU)<-gsub(".txt.*","",rownames(SBU))
rownames(SBU)<-gsub(".*miRNA/","",rownames(SBU))
WU<-t(WU)
rownames(WU)<-gsub("US84103569_","",rownames(WU))
rownames(WU)<-gsub("S01_miRNA_107_Sep09_","",rownames(WU))
rownames(WU)<-gsub(".txt.*","",rownames(WU))
rownames(WU)<-gsub(".*miRNA/","",rownames(WU))
#############ANOVA #####################################################
metadata<-read.csv("ANOVA/Copy of microarray table1.csv",stringsAsFactors=F,sep="\t")
metadata$Race[metadata$Race=="Caucasian"]<-"White"
metadata<-metadata[paste(metadata[,1],metadata[,3],sep="_") %in% c(rownames(SBU)),]
metadata$Race<-factor(as.character(metadata$Race))
SBU.data<-data.frame(SBU,Race=metadata[,4],Tumor=metadata[,5])
SBU.data<-SBU.data[order(SBU.data$Race),]
SBU.data.save.save<-data.frame(SBU.data,id=factor(rep(1:30,each=2)))
num.col<-ncol(SBU.data.save.save)
out<-array(0,c(3,num.col-3))

rownames(out)<-c("Race","Tumor","Int")
colnames(out)<-colnames(SBU)[1:(num.col-3)]
pas<-readLines("pas.SBU") #detectable list of miRNA names

####bootstrap rmANOVA pvalues#################3333
for(i in 1:961)
{
  if(colnames(WU)[i] %in% pas)
  {
    data<-SBU.data.save.save[,c(i,(num.col-2):num.col)]
    v<-tapply(data[,1],data[,2],var)
    change=sum(abs(log(v[1]/v[2]))<0.05)
    r<-tapply(data[,1],data$Race,mean)
    t<-tapply(data[,1],data$Tumor,mean)
    fil<-max(abs(r[1]-r[2]),abs(t[1]-t[2]))
    if(fil>log2(1.5))
    {
      data.list<-vector("list",1000)
      data.list<-mclapply(data.list,function(x){x<-data})
      data.list<-mclapply(data.list,resample,change=change)
      F<-get.p(data)
      pv<-get.p1(data)
      Fv<-mclapply(data.list,get.p)
      p<-array(0,c(3,1000))
      for(j in 1:1000)
      {
        p[,j]<-Fv[[j]]
      }
      pvalue<-pv
      for(j in 1:3)
      {
        
        pvalue[j]<-quant(F[j],p[j,])
        
      }
    }else{out[,i]=c(3,3,3)}
    out[,i]<-pvalue
  }else{out[,i]=c(2,2,2)}
  cat(c("SBU",i))
}
write.csv(t(out),"SBU.anova.pvalue.csv",row.names=T)