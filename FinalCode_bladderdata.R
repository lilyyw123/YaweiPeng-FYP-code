#############################################################
# Final code for bladder data analysis
# create on 2022/12/17
# By Yawei Peng
############################################################
setwd("E:/Lily University/FYP/MyProgress/Joint_frailty_Model_code/Code_and_Data_lily_mod/FYP_implement")

library(survival)
library(matrixcalc)
library(ggplot2)


############################################################
## Data prepocessing
############################################################

data(cancer,package = 'survival')
bld1 = bladder1

head(bld1)
names(bld1)
################################################
# meaning of the data

# id:	Patient id
# treatment:	Placebo, pyridoxine (vitamin B6), or thiotepa

# number:	Initial number of tumours (8=8 or more)
# size:	Size (cm) of largest initial tumour
# recur:	Number of recurrences
# start,stop:	The start and end time of each time interval

# status:	End of interval code, 0=censored, 1=recurrence,
#     2=death from bladder disease, 3=death other/unknown cause
# rtumor:	Number of tumors found at the time of a recurrence
# rsize:	Size of largest tumor at a recurrence

# enum:	Event number (observation number within patient)
########################################################
##########################################################
### data prepocessing#####
### step 1: deal with char.

### step 2: select previous 4 gap times's individual only -->
## assume the bladder cancer only happen 4
## if recurrent more than 4--> regard as censor

### step 3: decide terminal event: "event"
## status = 0 -- censor
## status = 1 -- recurrence
## status = 2-- death from cancer-- assume recurrence
## status = 3 -- death from other reason --- assume censor

### step4: standardlize continuous data

bld1$rtumor[bld1$rtumor=='.']='0'
bld1$rtumor = as.integer(bld1$rtumor)
bld1$rsize[bld1$rsize=='.']='0'
bld1$rsize = as.integer(bld1$rsize)

summary(bld1)
View(bld1)

# ## find terminal event
nn = dim(bld1)[1]
rank = seq(1,nn,by=1)
bld1$rank = rank
bld1[1:10,]

last = aggregate(bld1$rank~bld1$id,FUN = max)
colnames(last) = c("id","rank")
# event = bld1$status[last$rank]

# last = data.frame(last,event = event)
# 
bld1 = merge(bld1,last,by = "id",all.x=T)
bld1[1:10,]

### update status

bld1$status1 = bld1$status
bld1$status1[which(bld1$status==2)]=1 # die from bladder cancer--> treat as recurrent
bld1$status1[which(bld1$status==3)]=0 # die from other reason--> treat as recurrent


# # ser revent indicator
# bld1$revent=0
# bld1$revent[which(bld1$status==1)]=1 # recurrent event

# treatment
bld1$treatmentPY = 0
bld1$treatmentTH=0
bld1$treatmentPY[bld1$treatment=='pyridoxine']=1
bld1$treatmentTH[bld1$treatment=='thiotepa']=1

# incoperate covariate bldnumber and size

datamid = bld1
datamid$newnum = bld1$number
datamid$newsize = bld1$size

dim(bld1[which(bld1$enum==1),])

patient = length(which(bld1$enum==1)) # total number
i=1
for(i in 1:patient){
  
  group = datamid[which(datamid$id==i),]
  groupsize = dim(group)[1] # total number of occurrence
  
  if (groupsize>1){
    firstrec = which.min(group$rank.x) # find the start recurr
    while(firstrec<groupsize){
      group$newnum[firstrec+1] = group$rtumor[firstrec]
      group$newsize[firstrec+1] = group$rsize[firstrec]
      firstrec = firstrec+1
    }
  }
  datamid[which(datamid$id==i),] = group
}

# View(datamid)
# write.csv(datamid,"datamid.csv")

### use new size & new num tumor as covariate

bld1 = datamid

max(bld1$newnum)

## change scale--> on the scale 0/1
## devide by 8
bld1$newnum = bld1$newnum/8
bld1$newsize = bld1$newsize/8

length(which(bld1$enum==1))
# View(bld1)
# write.csv(bld1,"newbladder.csv")

#### 



####################################################
## data visualization
####################################################

fstatus = as.factor(bld1$status)
M = table(fstatus)
M
pie(M)
b = barplot(M,legend.text = c("censored ","recurrence","death from bladder disease","death from other reason"), col = c(2,4,1,7),
        main= "Bar chart of event status",ylab ="frequency"  )

barplot(M,main = 'Bar chart of number of recurrence',xlab = 'number of recurrence',ylab = 'relative frequency',col = c(rep(2,5),rep(7,5)))



#######################################################
#      combine previous 4 gaps

#########################################################

bld4gap = bld1[which(bld1$enum<5),]
# bld4gap = bld4gap[which((bld4gap$stop-bld4gap$start)>0),]  # delete the patient with zero gaptime
dim(bld4gap)
nn2 = dim(bld4gap)[1]

bld4gap$gap1 = 0
bld4gap$gap2 = 0
bld4gap$gap3 = 0

bld4gap$mi = 1*(bld4gap$recur>0)
head(bld4gap)
stop = bld4gap$stop/10
start =  bld4gap$start/10
gapp = stop-start

bld4gap$start = start
bld4gap$stop = stop
bld4gap$gapp = gapp

head(bld4gap)

################################
## first 4 gap's patients
###############################

#################################
#### store data--> look############
# names(bld4gap)
j = 1
for(i in 1:118){
  # every indivivdual
  J = length(which(bld4gap$id==i))
  # if(J==1){
  #   bld4gap$gap1[j] = gapp[j]
  # }
  if(J==2){
    # second time
    bld4gap$gap1[j+1] = gapp[j]
  }
  if(J==3){
    # second time
    bld4gap$gap1[j+1] = gapp[j]
    # third time
    bld4gap$gap1[j+2] = gapp[j]
    bld4gap$gap2[j+2] = gapp[j+1]
  }
  if(J==4){
    # second time
    bld4gap$gap1[j+1] = gapp[j]
    # third time
    bld4gap$gap1[j+2] = gapp[j]
    bld4gap$gap2[j+2] = gapp[j+1]
    # fourth time
    bld4gap$gap1[j+3] = gapp[j]
    bld4gap$gap2[j+3] = gapp[j+1]
    bld4gap$gap3[j+3] = gapp[j+2]
  }
  j = j+J
}
head(bld4gap)
bld4gap[80:90,]
# write.csv(bld4gap,"look.csv")

# bld4gap = read.csv("look.csv",header = T)

x1 = bld4gap$treatmentPY
x2 = bld4gap$treatmentTH
x3 = bld4gap$newnum
x4 = bld4gap$newsize

g1 = bld4gap$gap1
g2 = bld4gap$gap2
g3 = bld4gap$gap3

mi = bld4gap$mi
status1 = bld4gap$status1


######

firstgap = bld4gap[which(bld4gap$enum==1),]
dim(firstgap)
# write.csv(firstgap,'111.csv')
summary(firstgap)
summary(firstgap[which(firstgap$mi==1),])

#########################################################
## ## Proportional Hazards Cured Model
#########################################################

### start MLE
ll1 = L1 = 0
### loglikelihood 
loglik2 = function(x){
  
  lambda = abs(x[1])
  alpha = abs(x[2])  # if cannot solve, change to exp
  b1 = x[3]
  b2 = x[4]
  b3 = x[5]
  b4 = x[6]
  ## coeff for gaptime
  b5 = x[7]
  b6 = x[8]
  b7 = x[9] # previous 3 is enough
  
  d0 = x[10]
  d1 = x[11]  # trts
  d2 = x[12]
  d3 = x[13]  # gap time
  
  loghaz = log(alpha)+alpha*log(lambda)+(alpha-1)*log(gapp) +(b1*x1) + (b2*x2) +(b3*x3) +(b4*x4) + (b5*g1)+(b6*g2)+(b7*g3)
  cumhaz = (lambda*gapp)^alpha * exp((b1*x1) + (b2*x2) +(b3*x3) +(b4*x4) + (b5*g1)+(b6*g2)+(b7*g3))
  
  pi = 1/(1+exp(d0+(d1*x1)+(d2*x2)+(d3*gapp)))
  
  for(i in 1:nn2){
    
    if(mi[i]==1){
      # not a zero-recurrence
      if(status1[i]==1){
        # recurrence
        L1 = log(pi[i])+loghaz[i]-cumhaz[i]
      }else{
        # censored
        L1 = log(pi[i])-cumhaz[i]
      }
    }else{
      L1 = log(1-pi[i]+pi[i]*exp(-cumhaz[i]))
    }
    ll1 = ll1+L1
  }
  f1 = -ll1
  return(f1)
}

## can find solution!!

x = rep(0.1,13)
loglik2(x)
rrr2.1<- optim(x,loglik2, hessian = T, control = list(maxit = 20000),method="BFGS")
rrr2.1


##################################################################
######### MLE coefficients & statistical inference ###############
theta = rrr2.1$par  ## MLE: theta
hessian_theta = rrr2.1$hessian # since -loglike, var no need -
var_theta = diag(solve(hessian_theta))
se_theta = sqrt(var_theta)
loglike_theta = -rrr2.1$value

wald.theta  = theta/se_theta
pvalue.theta = pnorm(wald.theta,lower.tail = F)

result = data.frame(theta = round(theta,2), se_theta = round(se_theta,2), waldstat=round(wald.theta,2),pvalue= round(pvalue.theta ,2))
write.csv(result,"model1MLE_v3.csv")


#########################################################
## ## Proportional Hazards Model without Cured proportion
#########################################################
### start MLE

x1 = bld4gap$treatmentPY
x2 = bld4gap$treatmentTH
x3 = bld4gap$newnum
x4 = bld4gap$newsize

g1 = bld4gap$gap1
g2 = bld4gap$gap2
g3 = bld4gap$gap3

mi = bld4gap$mi
status1 = bld4gap$status1

ll2 = L2 = 0

### loglikelihood 
loglik3 = function(x){
  
  lambda = abs(x[1])
  alpha = abs(x[2])  # if cannot solve, change to exp
  b1 = x[3]
  b2 = x[4]
  b3 = x[5]
  b4 = x[6]
  ## coeff for gaptime
  b5 = x[7]
  b6 = x[8]
  b7 = x[9] # previous 3 is enough
  
  loghaz = log(alpha)+alpha*log(lambda)+(alpha-1)*log(gapp) +(b1*x1) + (b2*x2) +(b3*x3) +(b4*x4) + (b5*g1)+(b6*g2)+(b7*g3)
  cumhaz = (lambda*gapp)^alpha * exp((b1*x1) + (b2*x2) +(b3*x3) +(b4*x4) + (b5*g1)+(b6*g2)+(b7*g3))
  
  for(i in 1:nn2){
      # not a zero-recurrence
    if(status1[i]==1){
      # recurrence
      L2 = loghaz[i]-cumhaz[i]
    }else{
        # censored
      L2 = -cumhaz[i]
    }
    ll2 = ll2+L2
  }
  f1 = -ll2
  return(f1)
}

## can find solution!!

x = rep(0.1,9)
loglik3(x)
rrr3.1<- optim(x,loglik3, hessian = T, control = list(maxit = 20000),method="BFGS")
rrr3.1


##################################################################
######### MLE coefficients & statistical inference ###############
theta = rrr3.1$par  ## MLE: theta
hessian_theta = rrr3.1$hessian # since -loglike, var no need -
var_theta = diag(solve(hessian_theta))
se_theta = sqrt(var_theta)
loglike_theta = -rrr3.1$value

wald.theta  = theta/se_theta
pvalue.theta = pnorm(wald.theta,lower.tail = F)
result = data.frame(theta = round(theta,2), se_theta = round(se_theta,2), waldstat=round(wald.theta,2),pvalue= round(pvalue.theta ,2))

# result = data.frame(theta = theta, waldstat=wald.theta,pvalue = pvalue.theta)
write.csv(result,"model_nocure_MLE.csv")

##### remark! same as AFT (survreg) result!!!






