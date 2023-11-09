#*********   Markov Bridges **********#
# Reversed
# Bisection
# Uniformization
# Rejection
# Modified Rejection
# Direct
#************************************
# R Libraries
if(!require(ReIns))
{install.packages("ReIns", dependencies = TRUE)}
if(!require(Matrix))
{install.packages("Matrix", dependencies = TRUE)}
if(!require(matrixcalc))
{install.packages("matrixcalc", dependencies = TRUE)}
if(!require(expm))
{install.packages("expm", dependencies = TRUE)}
if(!require(dplyr))
{install.packages("dplyr", dependencies = TRUE)}
if(!require(ggplot2))
{install.packages("ggplot2", dependencies = TRUE)}
if(!require(patchwork))
{install.packages("patchwork", dependencies = TRUE)}
#***********************************#
set.seed(19)

#***********************************#
#_______ STATIONARY TIMES______#
#***********************************#
#dimension of the matrices
n=3:20

nest<-NULL
for(i in 1:length(n))
{
  Q=gen_inf(n[i])
  nest<-c(nest,timeest(Q))
}
time_stat<-data.frame(n,nest)
colnames(time_stat)<-c("n","Stationary times")
# graphic
gts<-time_stat %>%
  ggplot(aes(x=n,y=`Stationary times`)) 
gts + geom_line() +geom_point()+
  scale_x_continuous(limits = c(min(n), max(n)),breaks = n)+
  My_Theme

#***********************************#
#_______ JUMPS & TIMES _____________#
#***********************************#
ts=seq(1,10,by=0.01)
tablafinalEdos<-tablafinalTiempos<-NULL
test<-NULL

# Analysis for all matrices
for(i in 1:length(n))
{
  print(i)
  # Generador infinitesimal
    Q=gen_inf(n[i])
    a=1
    b=2
    t<-timeest(Q)
    test<-c(test,t)
  mc=1000
  SS=RVSS(Q,ceiling(t)+4,a,b)
  SSR=ESSREV(SS,a,b,Q,ceiling(t)+4,mc)
#  SSR
 tablafinalEdos<-cbind(tablafinalEdos,SSR[[1]])
 tablafinalTiempos<-cbind(tablafinalTiempos,SSR[[2]])
}

#_____________Results ____________#
metodos<-c("Bisection","Direct","Reverse",
           "Rejection", "Modified","Uniformization")
# States
finaledos<-data.frame(rep(n,6),as.vector(c(tablafinalEdos[4,],
                                       tablafinalEdos[3,],
                                       tablafinalEdos[1,],
                                       tablafinalEdos[5,],
                                       tablafinalEdos[6,],
                                       tablafinalEdos[2,])),
                      gl(6,length(n),label=metodos)) 
colnames(finaledos)<-c("Time","Mean","Method")

gedos<-finaledos %>%
  ggplot(aes(x=Time,y=Mean,color=Method)) 

gedos + geom_line() +geom_point()+
  ylab("Norm")+xlab("n")+labs(title = "Number of Jumps") + 
  scale_x_continuous(limits = c(min(n), max(n)),breaks = n)+
  My_Theme


# Stay times
finaltime<-data.frame(rep(n,6),as.vector(c(tablafinalTiempos[4,],
                                           tablafinalTiempos[3,],
                                           tablafinalTiempos[1,],
                                           tablafinalTiempos[5,],
                                           tablafinalTiempos[6,],
                                           tablafinalTiempos[2,])),
                      gl(6,length(n),label=c("Bisection","Direct","Reverse",
                                             "Rejection", "Modified","Uniformization"))) 
colnames(finaltime)<-c("Time","Mean","Method")

gtime<-finaltime %>%
  ggplot(aes(x=Time,y=Mean,color=Method)) 

gtime + geom_line() +geom_point()+
  ylab("Norm")+xlab("n")+
  labs(title = "Stay times") + 
  scale_x_continuous(limits = c(min(n), max(n)),breaks = n)+
  My_Theme

#***********************************#
#__________ ESTIMATION______________#
#***********************************#
# Model 1
n<-3
Q=gen_inf(n)
te<-timeest(Q) # tiempo estacionario

mc=1000
a=1
b=2
ite<-10
matl<-NULL
vecl<-NULL
for(i in 1:ite)
{
  resest<-ELU(a,b,Q,te,mc)
  matl<-cbind(matl,as.vector(resest[[1]]))
  vecl<-cbind(vecl,as.vector(resest[[2]]))
}
matlem<-apply(matl,1,mean)
matlesd<-apply(matl,1,sd)

veclem<-apply(vecl,1,mean)
veclesd<-apply(vecl,1,sd)

# reales
SS=RVSS(Q,te,a,b)
SS$NJ

# compara edos
metodos<-c("Reverse","Uniformization","Direct","Bisection", 
           "Rejection", "Modified Rejection")
ordm<-c(rep(metodos[1],n),rep(metodos[2],n),rep(metodos[3],n),
  rep(metodos[4],n),rep(metodos[5],n),rep(metodos[6],n))

lu<-data.frame(c(rep(SS$NJ[,1],6),rep(SS$NJ[,2],6),rep(SS$NJ[,3],6)),
               rep(ordm,n),matlem,matlesd)
colnames(lu)<-c("Real","Method","Estimation", "IC")
luj<-lu[lu$Real!=0,]
luj

# Comparison of times
luz<-data.frame(rep(SS$TS,6),ordm,veclem,veclesd)
colnames(luz)<-c("Real","Method","Estimation", "IC")
luz
  
  
#************************************#
#_____________ EXECUTION TIME________#
#************************************#
set.seed(2023)
# EXAMPLE 1
n<-3
Q=gen_inf(n)
te<-timeest(Q) # tiempo estacionario

# for t>t^* 
T2= seq(ceiling(te),ceiling(te)+5,by=1)
TiemposEjecucion<-checar_eje(n,Q,T2)
TiemposEjecucion

finalTE1<-data.frame(rep(T2,6),as.vector(c(TiemposEjecucion[4,],
                                        TiemposEjecucion[3,],
                                        TiemposEjecucion[1,],
                                        TiemposEjecucion[5,],
                                        TiemposEjecucion[6,],
                                        TiemposEjecucion[2,])),
                  gl(6,6,label=c("Bisection","Direct","Reverse",
                                 "Rejection", "Modified","Uniformization"))) 
colnames(finalTE1)<-c("Time","Mean","Method")

gTE1<-finalTE1 %>%
  ggplot(aes(x=Time,y=Mean,color=Method)) 

gTE1 + geom_line() +geom_point()+
  ylab("CPU Time")+xlab("Time between end-points")+
  scale_x_continuous(limits = c(min(T2), max(T2)),breaks = T2)+
  My_Theme



#********** EXAMPLE 2 ********#
n=3
dim<-n
Q<-matrix(0,nrow=dim,ncol=dim)
Q<-matrix(c(-2, 1,1, 
            0,-10,10,
            4,1,-5),ncol=dim,byrow=T)

 # stationary time
te<-timeest(Q)
# t>t^*
T2= seq(ceiling(te),ceiling(te)+5,by=1)
TiemposEjecucionE<-checar_eje(n,Q,T2)
TiemposEjecucionE<-as.matrix(TiemposEjecucionE)

final<-data.frame(rep(T2,6),as.vector(c(TiemposEjecucionE[4,],
                        TiemposEjecucionE[3,],
                        TiemposEjecucionE[1,],
                        TiemposEjecucionE[5,],
                        TiemposEjecucionE[6,],
                        TiemposEjecucionE[2,])),
                  gl(6,6,label=c("Bisection","Direct","Reverse",
                                 "Rejection", "Modified","Uniformization"))) 
colnames(final)<-c("Time","Mean","Method")

# Graphic   
g<-final %>%
  ggplot(aes(x=Time,y=Mean,color=Method)) 

g + geom_line() +geom_point()+
  ylab("CPU Time")+xlab("Time between end-points")+
  scale_x_continuous(limits = c(min(T2), max(T2)),breaks = T2)+
  My_Theme





#__________________________________________#
# ************ FUNCTIONS *********
#__________________________________________#

#************************#
#### Enlarged legend
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16))

#************************#
#### Stationary time
timeest<-function(Q)
{
  n<-dim(Q)[1]
  da=rep(1,n)
  a1<-rbind(t(Q),da)
  b1<-c(rep(0,n),1)
  DistEst<-qr.solve(a1,b1)
  matdisest<-t(matrix(rep(DistEst,n),nrow=n))
  
  ban<-0
  k<-1
  tol<-0.01
  ts=seq(1,10,by=0.01)
  while(ban==0)
  {
    Xn<-expm(ts[k]*Q)
    if(norm(Xn-matdisest)<tol) ban<-1
    k<-k+1
  }
  tiempoesta=ts[k-1]
  return(tiempoesta)
}

#************************#
#### Execution time
checar_eje<-function(dim,Q,T2)
{
  ite<-100
  T1=0
  TiemposR<-TiemposU<-TiemposD<-TiemposB<-TiemposRE<-TiemposREM<-NULL
  
  # loop
  for(i in 1:ite)
  { 
    ayutiR<-ayutiU<-ayutiD<-ayutiB<-ayutiRE<-ayutiREM<-NULL
    
   for(j in 1:length(T2))
    {
    ctiR<-ctiU<-ctiD<-ctiB<-ctiRE<-ctiREM<-NULL 
    for(a in 1:dim)
    { for (b in 1:dim)
     { 
    
     if(a!=b)
      {
       # REVERSED
      t<-proc.time()
      Puente= bridge_MJP_rev(a,b,Q,T2[j]-T1)
      ti<-list(proc.time()-t)
      ctiR<-c(ctiR,ti[[1]][1])
      
      #UNIFORMIZATION
      t<-proc.time()
      Puente=Uniformization(a,b,Q,T2[j]-T1)  
      ti<-list(proc.time()-t)
      ctiU<-c(ctiU,ti[[1]][1])
      
      #DIRECT
      t<-proc.time()
      Puente=dir_bri(a,b,Q,T2[j]-T1) 
      ti<-list(proc.time()-t)
      ctiD<-c(ctiD,ti[[1]][1])
      
      #BISECTION
      t<-proc.time()
      Puente=Bisection(a,b,T1,T2[j],Q,itera=0)
      ti<-list(proc.time()-t)
      ctiB<-c(ctiB,ti[[1]][1])
     
      #REJECTED
      t<-proc.time()
      Puente= bridge_MJP_REJ(a,b,Q,T2[j]-T1)
      ti<-list(proc.time()-t)
      ctiRE<-c(ctiRE,ti[[1]][1])
      
      # MODIFIED REJECTED
      t<-proc.time()
      Puente= bridge_MJP_MR(a,b,Q,T2[j]-T1)
      ti<-list(proc.time()-t)
      ctiREM<-c(ctiREM,ti[[1]][1])
     }#if
     }#for b
    }# for a
    ayutiR<-cbind(ayutiR,mean(ctiR))
    ayutiU<-cbind(ayutiU,mean(ctiU))
    ayutiD<-cbind(ayutiD,mean(ctiD))
    ayutiB<-cbind(ayutiB,mean(ctiB))
    ayutiRE<-cbind(ayutiRE,mean(ctiRE))
    ayutiREM<-cbind(ayutiREM,mean(ctiREM))
  }#  j
    TiemposR<-rbind(TiemposR,ayutiR)
    TiemposU<-rbind(TiemposU,ayutiU)
    TiemposD<-rbind(TiemposD,ayutiD)
    TiemposB<-rbind(TiemposB,ayutiB)
    TiemposRE<-rbind(TiemposRE,ayutiRE)
    TiemposREM<-rbind(TiemposREM,ayutiREM)
  
  }# i
  
  MediaTiemposR<-apply(TiemposR,2,mean)
  MediaTiemposU<-apply(TiemposU,2,mean)
  MediaTiemposD<-apply(TiemposD,2,mean)
  MediaTiemposB<-apply(TiemposB,2,mean)
  MediaTiemposRE<-apply(TiemposRE,2,mean)
  MediaTiemposREM<-apply(TiemposREM,2,mean)
  
 Tablafinal<- rbind(MediaTiemposR,MediaTiemposU,MediaTiemposD,
                    MediaTiemposB,MediaTiemposRE,MediaTiemposREM)
 colnames(Tablafinal)<-T2
 return(Tablafinal)
}


#************************#
####  Infinitesimal generador
gen_inf=function(n)
{
  Q=matrix(1/(n-1),nrow = n,ncol=n)
  diag(Q)=-1
  return(Q)
}

#************************#
#### Integrals Iabcd
I_abcd=function(Q,t,a,b,c,d)
{
  n=dim(Q)[1]
  if(a!=c&&d!=b){
    Int=(1/n^2)*(t+t*exp(-(n*t)/(n-1))-((2*(n-1))/n)*(1-exp(-(n*t)/(n-1))))
  }else if(a==c&&d==b){
    Int=(1/n^2)*(t+((n-1)^2)*t*exp(-(n*t)/(n-1))+((2*(n-1)^2)/n)*(1-exp(-(n*t)/(n-1))))
  }else{
    Int=(1/n^2)*(t-(n-1)*t*exp(-(n*t)/(n-1))+(((n-2)*(n-1))/n)*(1-exp(-(n*t)/(n-1))))
  }
  
  return(Int)
}

#************************#
#### Integrals
RVSS=function(Q,t,a,b)
{
  n=dim(Q)[1]
  MInt=matrix(0,nrow=n,ncol=n)
  SS=matrix(0,nrow=n,ncol=n)
  Pab=expm(t*Q)[a,b]
  for(c in 1:n)
  {
    for(d in 1:n)
    {
      MInt[c,d]=I_abcd(Q,t,a,b,c,d)
      if(c==d){
        SS[c,d]=MInt[c,d]/Pab
      }else{
        SS[c,d]=Q[c,d]*MInt[c,d]/Pab
      }
    }
  }
  NJ=SS
  diag(NJ)=0
  TS=diag(SS)
  return(list(NJ=NJ,TS=TS))
}


#__________________________________________#
#### Generator of the inverse process
#__________________________________________#
QRF=function(Q)
{
  n=dim(Q)[1]
  da=rep(1,n)
  a1<-rbind(t(Q),da)
  b1<-c(rep(0,n),1)
  pi<-qr.solve(a1,b1)
  QR=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      QR[i,j]=pi[j]*Q[j,i]/pi[i]
    }
  }
  return(QR)
}


#__________________________________________#
# Generate the transition rate matrix
# @param Q: Intensity matrix
# @return P: Jump matrix
#__________________________________________#
JumpMat<-function (Q)
{
  q<-dim(Q)[1]
  P<-matrix(0,q,q)  
  for(i in 1:q)
  {if (Q[i,i]<0)
  {
    P[i,] <- -Q[i,]/Q[i,i]
    P[i,i]<-0}
    else
    {P[i,i]<-1}
  }
  return(P)
}


#__________________________________________#
#Simulate a path of a MJP.

#@param pi: initial distribution.
#@param Q: infinitesimal generator matrix.
#@param t: time horizon of path.
#@return (time,tray): a tuple containing:
#  length(time)=number of jumps + 2
#time: time_1=initial time, time_i= time of jump i (i in {2,3,...,length(time)-1),time_end= time horizon
#tray: tray_i= states visited with arrival time time_i for in {1,2,...,length(time)-1} and state_end= same that end-1 (do no jump)
#__________________________________________#
PathSim<-function (pi,Q,t)
{
  m<-length(pi)
  #Sampling Space
  E <-seq(1,m,1)
  #Jump matrix
  P<-JumpMat(Q)
  #Path
  tray<-sample(E,size=1,prob=pi)
  #Jumping times
  time<-0
  while (time[length(time)] <t)
  { end<-length(tray)
  newtray<-tray[length(tray)]
  if (Q[newtray,newtray]==0)
  {
    tray<-c(tray,newtray) 
    time<-c(time,t)
  }
  else
  {
    s<-rexp(1,-Q[newtray,newtray])
    time<-c(time,time[length(time)]+s)
    tray<-c(tray,sample(E,size=1,prob=P[newtray,]))
  }
  }
  time[length(time)]<-t
  tray[length(tray)]<-tray[length(time)-1]
  return(list(tray=tray,time=time))
}


#__________________________________________#
#Simulate a path of a MJP reverse time

#@param pi: initial distribution.
#@param Q: infinitesimal generator matrix.
#@param t: time horizon of path.
#@return (time,tray): a tuple containing:
#  length(time)=number of jumps + 2
#time: time_1=initial time,
#      time_i= time of jump i (i in {2,3,...,length(time)-1),
#      time_end= time horizon
#tray: tray_i= states visited with arrival time time_i for in {1,2,...,length(time)-1} 
#      state_end= same that end-1 (do no jump)
#__________________________________________#
RevPathSim<-function(pi,Q,t)
{
  QR=QRF(Q)
  path=PathSim(pi,QR,t)
  tray=c(rev(path$tray[-length(path$tray)]),path$tray[1])
  
  time= time=c(0,cumsum(rev(diff(path$time))))
  return(list(tray=tray,time=time))
}


#__________________________________________#
### Simulate a Reverse Markov Bridge 
# input
# a, b
#@param Q: infinitesimal generator matrix.
#@param t: time horizon of path.
#__________________________________________#

bridge_MJP_rev=function(a,b,Q,t)
{
  
  pi1=numeric(dim(Q)[1]) # inicializar
  pi1[a]=1
  pi2=numeric(dim(Q)[1]) # inicializar
  pi2[b]=1

  # Generate a forward path using a
  dpath=PathSim(pi1,Q,t)
  de=dpath$tray
  dt=dpath$time
  
  # Flag
  crit=0
  
  while(crit==0)
  {
    tray=NULL
    time=NULL
    LU<-NULL
    
    # Reversed Path using b
    rpath=RevPathSim(pi2,Q,t)
    re=rpath$tray
    rt=rpath$time
   
    # first jump coincides in both chains    
    if(de[1]==re[1])
    {
      tray=re
      time=rt
      LU="R"
      crit=1
      ban<-1
      break
    }
    
    # More jumps   
    jumpr=length(re) # jumps reversed 
    jumpd=length(de) # jumps forward
    
    nr=2 #counting jumps reversed
    nd=2 #counting jumps forward
    
    tray=de[1] #initializes tray
    time=dt[1] #initializes time
    LU<-c(LU,"F")
    
    while(nd<jumpd && nr<jumpr)
    {
      while(dt[nd]<rt[nr]) # Movemos la forward
      {  
        if(re[nr-1]==de[nd]) #R[2]=D[1]
        { 
          tray=c(tray,de[nd],re[(nr):jumpr])
          time=c(time,dt[nd],rt[(nr):jumpr])
          LU<-c(LU,"F",rep("R",length(rt[(nr):jumpr])))
          nr=jumpr
          nd=jumpd
          crit=1 # terminamos
          ban<-1
          break
        } # if 
        tray=c(tray,de[nd])
        time=c(time,dt[nd])
        LU<-c(LU,"F")
        nd=nd+1 
      } 
      
      ban<-1
      if(abs(dt[nd]-rt[nr])<0.0001)
      {  
        ban<-0
        break
      }
      
      # the other case
      if(ban==1)
      {
        while(rt[nr]<dt[nd])
      {  
        if(re[nr]==de[nd-1])
        {
          if(nr<(jumpr-1))
          {
            tray=c(tray,re[(nr+1):jumpr])
            time=c(time,rt[(nr+1):jumpr])
            LU<-c(LU,rep("R",length(rt[(nr+1):jumpr])))
            nr=jumpr
            nd=jumpd
            crit=1
            ban<-1
            break
          }#if 

          tray=c(tray,re[jumpr])
          time=c(time,rt[jumpr])
          LU<-c(LU,"R")
          nr=jumpr
          nd=jumpd
          crit=1
          ban<-1
          break
        } #if
        nr=nr+1
        ban<-1
      } # while (rt[nr]<dt[nd])
      } #if ban  
    } # while(nd<jumpd && nr<jumpr)
    
  }  # while(crit==0)
  return(list(tray=tray,time=time))
}

#__________________________________________#
#Calculate the sufficient statistics
#@param m: number of states
#@param time: the start time, jump times and horizon time
#@param tray: state visited and final state
#@return (R,MS): a tuple containing:
#  R: total time spent in each state,
#MS: number of observed jumps from one state to
#another.
#__________________________________________#
MJPSS<-function (m,time,tray)
{
  jmp<-length(tray)-2
  R<-rep(0,m)
  MS<-matrix(0,m,m)
  for (i in 1:jmp)
  {
    ini<-tray[i]
    fin<-tray[i+1]
    MS[ini,fin]<-MS[ini,fin]+1
    R[ini]<-R[ini]+time[i+1]-time[i]
  }
  R[tray[jmp+1]]<-R[tray[jmp+1]]+time[jmp+2]-time[jmp+1]
  return(list(vector=R,matrix=MS))
}


#**************************************** #
NN<-function(a,b,gen_inf,HT){
  u=runif(1)
  z=0
  c=gen_inf
  m=-min(c)
  R=diag(length(c[1,]))+c/m
  Pab=expm(gen_inf*HT)
  
  if(a!=b){
    i=1
  }else{
    i=0
  }
  Z=0
  aux=0
  if(a!=b){
    i=1
  }else{
    i=0
  }
  po=(exp(-m*HT)*((m*HT)^i)/factorial(i))
  con=((matrix.power(R,i)[a,b])/Pab[a,b])
  aux=po*con
  z=aux+z
  
  while(z<u){
    i=i+1
    po=(exp(-m*HT)*((m*HT)^i)/factorial(i))
    con=((matrix.power(R,i)[a,b])/Pab[a,b])
    aux=po*con
    z=aux+z
  }
  return(i)
}


#**************************************** #
markovchain2<-function(a,b,R,n){
  c=a
  if(n<2){
    return(c(a,b))
  }
  estados=a
  while( length(estados)<=n ){ #tail(estados,1)!=b ||
    c=sample(1:length(R[,1]),size=1,prob= R[c,])
    estados=c(estados,c)
  }
  return(c(estados))
}

#**************************************** #
unifor<-function(a,b,gen_inf,HT){
  c=gen_inf
  m=-min(c)
  R=diag(length(c[1,]))+c/m
  n=NN(a,b,gen_inf,HT)
  edos=markovchain2(a,b,R,n)
  while(tail(edos,1)!=b){
    edos=markovchain2(a,b,R,n)
  }
  n=length(edos)-1
  if((n==0 && a==b) || (n==1 && a==b) ){
    edos=a
    return(list(edos,c(0,HT),gen_inf))
  }
  if((n==1 && a!=b) || (n==0 && a!=b) ){
    edos=c(a,b)
    return(list(edos,c(0,HT*runif(n),HT),gen_inf))
  }
  if(n>=2){
    saltos=HT*runif(n)
    salt=sort(saltos)
    tiempos=c(0,salt,HT)
    E=1:length((gen_inf[1,]))
    gen_inf
    return(list(edos,tiempos,gen_inf,R))
  }
}

#**************************************** #
Uniformization<-function(a,b,gen_inf,HT){
  options(warn=-1)
  t<-proc.time()
  w=unifor(a,b,gen_inf ,HT)
  i=1
  while(i <=length(w[[1]])){
    i=i+1
    while(w[[1]][i]==w[[1]][i-1] && i <=length(w[[1]])){
      w[[1]]=w[[1]][-i]
      x=w[[2]][i-1]+(w[[2]][i]-w[[2]][i-1])
      x
      w[[2]]=w[[2]][-i]
      w[[2]][i]=x
    }
  }
  if(tail(w[[2]],1)!=HT){
    w[[2]][length(w[[2]])]=HT
  }
  return(list(w,proc.time()-t))
}

#**************************************** #
# Probabilities of no jumping
pr_nj=function(Q,delta)
{
  pr=exp(diag(Q)*delta)/diag(expm(Q*delta))
  return(pr)
}
pr_tj=function(a,b,U,D,UI,Q,delta)
{
  n=dim(Q)[1]
  ja=Jaj(a,Q,D,delta)
  Pab=expm(Q*delta)[a,b]
  pr=numeric(n)
  for(i in 1:n)
  {
    if(i==a){
      pr[i]=0
    }else{
      pr[i]=Q[a,i]/Pab
      aux=0
      for(j in 1:n){
        aux=aux+U[i,j]*UI[j,b]*ja[j]
      }
      pr[i]=pr[i]*aux
    }
  }
  pr[a]=0
  pr=pr/sum(pr)
  return(pr)
}

#**************************************** #
fipi=function(t,i,a,b,Q,pi,Dd,Ui,UIb,delta)
{
  Qa=-Q[a,a]
  Pab=expm(Q*delta)[a,b]
  
  fi=(Q[a,i]/Pab)*(sum(Ui*UIb*exp(delta*Dd)*exp(-t*(Dd+Qa))))
  fi=fi/pi
  return(fi)  
}


#**************************************** #
# simulation of tau
sim_tau=function(i,a,b,Q,pi,Dd,Ui,UIb,delta)
{
  
  crit=0
  while(crit==0)
  {
    candi=rtexp(1,-Q[a,a],delta)
    fit=fipi(candi,i,a,b,Q,pi,Dd,Ui,UIb,delta)
    evaexpt=dtexp(candi,-Q[a,a],delta)
    u=runif(1)
    
    acep=fit/evaexpt
    if(acep>u&candi<delta){
      tau=candi
      crit=1
    }
  }
  return(tau)
}


#**************************************** #
Jaj=function(a,Q,D,delta)
{
  Qa=-Q[a,a]
  n=dim(Q)[1]
  Dd=diag(D)
  ja=numeric(n)
  for(j in 1:n){
    if((Dd[j]+Qa)==0){
      ja[j]=delta*exp(Dd[j]*delta)
    }else{
      ja[j]=(exp(Dd[j]*delta)-exp(-Qa*delta))/(Dd[j]+Qa)
    }
  }
  return(ja)
}

#**************************************** #
dir_bri=function(a,b,Q,delta)
{
  tray=a
  time=0
  crit=0
  
  n=dim(Q)[1]
  E=seq(1:n)
  U <- eigen(Q)$vectors
  D <- diag(eigen(Q)$values)
  UI=solve(U)
  UIb=UI[,b]
  Dd=diag(D)
  del_aux=delta
  
  while(crit==0)
  {
    pr=pr_nj(Q,del_aux)
    Z=rbinom(1,1,pr[a])
    
    if(a==b&Z==1){
      tray=c(tray,a)
      time=c(time,delta)
      crit=1
    }
    else{
      pis=pr_tj(a,b,U,D,UI,Q,del_aux)
      i=sample(E,1,prob=pis)
      Ui=U[i,]
      pi=pis[i]
      
      tau=sim_tau(i,a,b,Q,pi,Dd,Ui,UIb,del_aux)
      tray=c(tray,i)
      time=c(time,tail(time,1)+tau)
      del_aux=delta-tail(time,1)
      a=i
    }
  }
  return(list(tray=tray,time=time))
}

#********************************************#
Rab<-function(a,b,Q,HT){
  if(-Q[a,a]==-Q[b,b]){
    return(Q[a,b]*HT*exp(Q[a,a]*HT))
  }else{
    return(Q[a,b]*(exp(Q[a,a]*HT)-exp(Q[b,b]*HT))/(-Q[b,b]+Q[a,a]))
  }
}

#********************************************#
#********************************************#
Bisaa<-function(a,b,T1,T2,Q1,itera){
  HT=T2-T1
  HT0=HT/2
  n=length(Q1[1,])
  estados=1:n
  Po=expm(Q1*HT0)
  C=estados[-a]
  probab=NULL
  opciones1=NULL
  probba=c(exp(Q1[a,a]*HT0)*exp(Q1[a,a]*HT0),
           exp(Q1[a,a]*HT0)*(Po[a,a]-exp(Q1[a,a]*HT0)), ##2
           exp(Q1[a,a]*HT0)*(Po[a,a]-exp(Q1[a,a]*HT0)), ##3
           (Po[a,a]-exp(Q1[a,a]*HT0))*(Po[a,a]-exp(Q1[a,a]*HT0)))##4
  probab=c(probab,probba)
  opciones1=list(c(0,0),c(0,2),c(2,0),c(2,2))
  
  for(i in C ){
    probbc=c(Rab(a,i,Q1,HT0)*Rab(i,a,Q1,HT0), ##5
             Rab(a,i,Q1,HT0)*(Po[i,a]-Rab(i,a,Q1,HT0)), ##6
             Rab(i,a,Q1,HT0)*(Po[a,i]-Rab(a,i,Q1,HT0)), ##7
             (Po[a,i]-Rab(a,i,Q1,HT0))*(Po[i,a]-Rab(i,a,Q1,HT0)))##8
    probab=c(probab,probbc)
    opcionesc=list(c(1,1),c(1,2),c(2,1),c(2,2))
    opciones1=cbind(opciones1,opcionesc)
  }
  complete=c(a,C)
  if(itera==1){
    probab[1] <- 0
  }
  opcion=sample(1:length(probab),size=1,prob=abs(probab)/sum(abs(probab),na.rm=T))
  estado=floor(opcion/4)+1
  if(opcion%%4==0){
    estado=opcion/4
  }
  if(estado>n){
    estado=n
  }
  # interval 1
  if(opciones1[[opcion]][1]==0){
    estad1=c(a,a)
    tempos1=c(T1,T1+HT0)
    media1=cbind(estad1,tempos1)
    
    # interval 2
    if(opciones1[[opcion]][2]==0){
      estad2=c(a,a)
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      
      if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
      }
      estad2=c(complete[estado],b,b)
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      a=a
      b=a
      media2=Bisection(a,b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
      a=a
      b=b
    }
  }else if(opciones1[[opcion]][1]==1){
    
    if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
      time1=T1+HT0*runif(1)
    }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
      time1=T1+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
    }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
      time1=T1+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
    }
     estad1=c(a,complete[estado],complete[estado])
    tempos1=c(T1,time1,T1+HT0)
    media1=cbind(estad1,tempos1)
    # Intervalo 2
    
    if(opciones1[[opcion]][2]==0){
      estad2=c(complete[estado],b)
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
      }
       #time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]],endpoint=HT0)
      estad2=c(complete[estado],b,b)
      T2=T1+2*HT0
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      a=complete[estado]
      media2=Bisection(a,b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
    }
  }else if(opciones1[[opcion]][1]==2){
    T1=T1
    T2=T1+HT0
    a=a
    media1=Bisection(a,complete[estado],T1,T2,Q1,itera=1)
     
    ####intervalo2
    if(opciones1[[opcion]][2]==0){
      estad2=c(a,a)
      T2=T1+2*HT0
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      
      if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
      }
      #time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]],endpoint=HT0)
      estad2=c(complete[estado],b,b)
      T2=T1+2*HT0
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      media2=Bisection(complete[estado],b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
    }
  }
}

#********************************************#
#********************************************#
Bisab<-function(a,b,T1,T2,Q1,itera){
  
  HT=T2-T1
  HT0=HT/2
  n=length(Q1[1,])
  estados=1:n
  Po=expm(Q1*HT0)
  indi=which(Q1[a,]>0)
  estados=1:n
  C=estados[-c(a,b)]
  probab=NULL
  opciones1=NULL
  opciones2=NULL
  probb=c(exp(Q1[a,a]*HT0)*Rab(a,b,Q1,HT0),
          exp(Q1[a,a]*HT0)*(Po[a,b]-Rab(a,b,Q1,HT0)), ##2
          Rab(a,b,Q1,HT0)*(Po[a,a]-exp(Q1[a,a]*HT0)), ##3
          (Po[a,a]-exp(Q1[a,a]*HT0))*(Po[a,b]-Rab(a,b,Q1,HT0)))##4
  probab=c(probab,probb)
  opciones1=list(c(0,1),c(0,2),c(2,1),c(2,2))
  
  probb2=c(exp(Q1[b,b]*HT0)*Rab(a,b,Q1,HT0),
           Rab(a,b,Q1,HT0)*(Po[b,b]-exp(Q1[b,b]*HT0)), ##2
           exp(Q1[b,b]*HT0)*(Po[a,b]-Rab(a,b,Q1,HT0)), ##3
           (Po[b,b]-exp(Q1[b,b]*HT0))*(Po[a,b]-Rab(a,b,Q1,HT0)))##4
  probab=c(probab,probb2)
  
  opciones2=list(c(1,0),c(1,2),c(2,0),c(2,2))
  opciones1=cbind(opciones1,opciones2)
  
  for(i in C ){
    probbc=c(Rab(a,i,Q1,HT0)*Rab(i,b,Q1,HT0), ##5
             Rab(a,i,Q1,HT0)*(Po[i,b]-Rab(i,b,Q1,HT0)), ##6
             Rab(i,b,Q1,HT0)*(Po[a,i]-Rab(a,i,Q1,HT0)), ##7
             (Po[a,i]-Rab(a,i,Q1,HT0))*(Po[i,b]-Rab(i,b,Q1,HT0)))##8
    
    probab=c(probab,probbc)
    opcionesc=list(c(1,1),c(1,2),c(2,1),c(2,2))
    opciones1=cbind(opciones1,opcionesc)
  }
  complete=c(a,b,C)
  
  if(itera==1){
    probab[1] <- 0
    probab[5]<- 0
  }
  opcion=sample(1:length(probab),size=1,prob=abs(probab)/sum(abs(probab),na.rm=T))
  estado=floor(opcion/4)+1
  if(opcion%%4==0){
    estado=opcion/4
  }
  if(estado>n){
    estado=n
  }
  if(opciones1[[opcion]][1]==0){
    estad1=c(a,a)
    tempos1=c(T1,T1+HT0)
    media1=cbind(estad1,tempos1)
    
    ####interval 
    if(opciones1[[opcion]][2]==0){
      estad2=c(a,a)
      T2=T1+2*HT0
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
      }
      estad2=c(a,b,b)
      T2=T1+2*HT0
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
      
    }
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      a=a
      b=b
      media2=Bisection(a,b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
    }
    
    #####elif primero
  }else if(opciones1[[opcion]][1]==1){
    if(Q1[complete[estado],complete[estado]]==Q1[complete[1],complete[1]]){
      time1=T1+HT0*runif(1)
    }else if(Q1[complete[estado],complete[estado]]<Q1[complete[1],complete[1]]){
      time1=T1+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[1],complete[1]],endpoint=HT0)
    }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[1],complete[1]])){
      time1=T1+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[1],complete[1]],endpoint=HT0)
    }
    #time1=T1+rtexp(1,rate=-Q1[complete[estado],complete[estado]],endpoint=HT0)
    estad1=c(a,complete[estado],complete[estado])
    tempos1=c(T1,time1,T1+HT0)
    media1=cbind(estad1,tempos1)
    ###Intervalo 2
    
    T2=T1+2*HT0
    
    ####intervalo2
    if(opciones1[[opcion]][2]==0){
      estad2=c(complete[estado],complete[estado])
      
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      ####-----version anterior  1 en lugar de 2 complete[2] ahora antes, complete[1]]
      if(Q1[complete[estado],complete[estado]]==Q1[complete[2],complete[2]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[2],complete[2]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[2],complete[2]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[2],complete[2]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[2],complete[2]],endpoint=HT0)
      }
      estad2=c(complete[estado],b,b)
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      a=complete[estado]
      
      media2=Bisection(complete[estado],b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
    }
  }else if(opciones1[[opcion]][1]==2){
    
    T1=T1
    T2=T1+HT0
    media1=Bisection(a,complete[estado],T1,T2,Q1,itera=1)
    T2=T1+2*HT0
    ####intervalo2
    if(opciones1[[opcion]][2]==0){
      estad2=c(b,b)
      tempos2=c(T1+HT0,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==1){
      #####-------version anterior 1 en lugar de 2
      if(Q1[complete[estado],complete[estado]]==Q1[complete[2],complete[2]]){
        time2=T1+HT0+HT0*runif(1)
      }else if(Q1[complete[estado],complete[estado]]<Q1[complete[2],complete[2]]){
        time2=T1+HT0+rtexp(1,rate=-Q1[complete[estado],complete[estado]]+Q1[complete[2],complete[2]],endpoint=HT0)
      }else if((Q1[complete[estado],complete[estado]])>(Q1[complete[2],complete[2]])){
        time2=T1+HT0+rtexp(1,rate=+Q1[complete[estado],complete[estado]]-Q1[complete[2],complete[2]],endpoint=HT0)
      }
      estad2=c(complete[estado],b,b)
      tempos2=c(T1+HT0,time2,T2)
      media2=cbind(estad2,tempos2)
      total=rbind(media1,media2)
      return(total)
    }
    if(opciones1[[opcion]][2]==2){
      T1=T1+HT0
      T2=T1+HT0
      
      media2=Bisection(complete[estado],b,T1,T2,Q1,itera=1)
      media2=rbind(media1,media2)
      return(media2)
    }
  }
}

#********************************************#
Bisection<-function(a,b,T1,T2,Q1,itera){
  n=length(Q1[,1])
  if(a==b){
    return(Bisaa(a,b,T1,T2,Q1,itera))
  }else{
    return(Bisab(a,b,T1,T2,Q1,itera))
  }
}

#**************************************** #
prob_salto=function(gen_inf)
{
  options(warn=-1)
  dim=dim(gen_inf)[1]
  P=matrix(0,dim,dim)
  for(i in 1:dim)
  {
    if(gen_inf[i,i]==0){
      P[i,i]=1
    }
    else{
      for(j in 1:dim)
      {
        if(i!=j)
        {
          P[i,j]=-gen_inf[i,j]/gen_inf[i,i]
        }
      }
    }
  }
  return((P))
}


#**************************************** #
tray_psm_a=function(a,gen_inf,HT)
{
  edos=NULL
  tiempos_saltos=NULL
  edos=a
  E=1:length(gen_inf[1,])
  tiempos_saltos=0
  tpaso=rexp(1,-gen_inf[tail(edos,1),tail(edos,1)])
  tiempos_saltos=c(tiempos_saltos,tpaso)
  P=prob_salto(gen_inf)
  while(tail(tiempos_saltos,1)< HT)
  {
    edos=c(edos,sample(x=E,size=1,prob=P[tail(edos,1),]))
    aux=tail(tiempos_saltos,1)
    tpaso=rexp(1,-gen_inf[tail(edos,1),tail(edos,1)])
    tiempos_saltos=c(tiempos_saltos,aux+tpaso)
  }
  tiempos_saltos[length(tiempos_saltos)]=HT
  return(list(edos,tiempos_saltos,gen_inf))
}

#**************************************** #
Forward.Samp=function(a,b,gen_inf,HT){ 
  t <- proc.time() # Inicia el cron?metro
  v=tray_psm_a(a,gen_inf,HT)
  while(tail(v[[1]],1)!=b){
    v=tray_psm_a(a,gen_inf,HT)
  }
  return(list(v,proc.time()-t ))
}


#**************************************** #
Rejtsmaux<-function(a,gen_inf,HT){
  edos=NULL
  tiempos_saltos=NULL
  edos=a
  E=1:length(gen_inf[1,])
  tiempos_saltos=0
  Qa=-gen_inf[tail(edos,1),tail(edos,1)]
  t=-log(1-runif(1)*(1-exp(-HT*Qa)))/Qa
  tiempos_saltos=c(tiempos_saltos,tail(tiempos_saltos,1)+t)
  P=prob_salto(gen_inf)
  
  while(tail(tiempos_saltos,1)< HT)
  {
    edos=c(edos,sample(x=E,size=1,prob=P[tail(edos,1),]))
    Qa=-gen_inf[tail(edos,1),tail(edos,1)]
    t=rexp(1,-gen_inf[tail(edos,1),tail(edos,1)])
    tiempos_saltos=c(tiempos_saltos,tail(tiempos_saltos,1)+t)
  }
  tiempos_saltos[length(tiempos_saltos)]=HT
  return(list(edos,tiempos_saltos,gen_inf))
}

#**************************************** #
Mod.Rej.Sam<-function(a,b,gen_inf,HT){
  t<-proc.time()
  if(a==b){
    return(Forward.Samp(a,b,gen_inf,HT))
  }else{
    v=Rejtsmaux(a,gen_inf,HT)
    while(tail(v[[1]],1)!=b){
      v=Rejtsmaux(a,gen_inf,HT)
    }
    return(list(v,proc.time()-t))
  }
}

#**************************************** #
Gen.Mod.Rej.Samp<-function(camino,tempos,gen_inf){
  t=proc.time()
  completo=NULL
  tiempos=NULL
  n=length(camino)
  for(i in 1:(n-1)){
    HT=-tempos[i]+tempos[i+1]
    a=camino[i]
    b=camino[i+1]
    Mod=Mod.Rej.Sam(a,b,gen_inf,HT)
    if(i==1){
      tiempos=c(tiempos,Mod[[1]][[2]])
      completo=c(completo,c(Mod[[1]][[1]],b))
      cambio=c(length(tiempos))
      
    }else{
      
      tiempos1=Mod[[1]][[2]]+tail(tiempos,1)[[1]]
      completo1=c(Mod[[1]][[1]],b)
      tiempos=c(tiempos,Mod[[1]][[2]]+tail(tiempos,1)[[1]])
      completo=c(completo,c(Mod[[1]][[1]],b))
      cambio=c(cambio,length(tiempos))
    }
  }
  return(list(completo,tiempos,proc.time()-t,cambio))
}

#******************************************
#### REJECTED
bridge_MJP_REJ=function(a,b,Q,t)
{
  n=dim(Q)[1]
  pi=numeric(n)
  pi[a]=1
  
  crit=0
  int=0
  while(crit<1)
  {
    path=PathSim(pi,Q,t)
    end=length(path$tray)
    if(b==path$tray[end]){
      tray=path$tray
      time=path$time
      crit=1
    }
    int=int+1
  }
  return(list(tray=tray,time=time,int=int))
}
#******************************************
#### MODIFIED REJECTED
bridge_MJP_MR=function(a,b,Q,t)
{
  n=dim(Q)[1]
  pi=numeric(n)
  
  E <-seq(1,n,1)
  P=JumpMat(Q)
  Qa=-Q[a,a]
  crit=0
  int=0
  
  while(crit<1)
  {
    if(a==b)
    {
      path=bridge_MJP_REJ(a,b,Q,t)
      tray=path$tray
      time=path$time
      int=path$int
      crit=1
    }else{
      u=runif(1)
      
      tau=-log(1-u*(1-exp(-t*Qa)))/Qa
      c=sample(E,size=1,prob=P[a,])
      
      pi[c]=1
      path=PathSim(pi,Q,t-tau)
      end=length(path$tray)
      if(b==path$tray[end]){
        tray=c(a,path$tray)
        time=c(0,tau+path$time)
        crit=1
      }
      int=int+1
    }
  }
  return(list(tray=tray,time=time,int=int))
}

#*****************************************#
# Estimation of the sufficient statistics
#### ELU
#*****************************************#
ELU=function(a,b,Q,t,mc)
{
  n=dim(Q)[1]
  MJ<-MJU<-MJD<-MJB<-MJRE<-MJREM<-matrix(0,nrow=n,ncol=n)
  ST<-STU<-STD<-STB<-STRE<-STREM<-numeric(n)
  
  for(i in 1:mc)
    {
    # REVERSE
    path=bridge_MJP_rev(a,b,Q,t)
    SS=MJPSS(n,path$time,path$tray)
    MJ=MJ+SS$matrix/mc
    ST=ST+SS$vector/mc
    
    #UNIFORMIZATION
    pathU= Uniformization(a,b,Q,t)  
    SSU=MJPSS(n,pathU[[1]][[2]],c(pathU[[1]][[1]],b))
    MJU=MJU+SSU$matrix/mc
    STU=STU+SSU$vector/mc

  # DIRECTO     
    pathD=dir_bri(a,b,Q,t)
    SSD=MJPSS(n,pathD$time,pathD$tray)
    MJD=MJD+SSD$matrix/mc
    STD=STD+SSD$vector/mc
    
  #Bisection
    pathB= Bisection(a,b,0,t,Q,itera=0)
    SSB=MJPSS(n,pathB[,2],pathB[,1])
    diag(SSB$matrix)<-0
    MJB=MJB+SSB$matrix/mc
    STB=STB+SSB$vector/mc
    
    #RECHAZO
    pathRE= bridge_MJP_REJ(a,b,Q,t) #Forward.Samp(a,b,Q,t) 
    SSRE=MJPSS(n,pathRE$time,pathRE$tray) #MJPSS(n,pathRE[[1]][[2]],c(pathRE[[1]][[1]],b))
    MJRE=MJRE+SSRE$matrix/mc
    STRE=STRE+SSRE$vector/mc
    
    #RECHAZO_MODIFICADO
    pathREM= bridge_MJP_MR(a,b,Q,t) # Mod.Rej.Sam(a,b,Q,t) 
    SSREM=MJPSS(n,pathREM$time,pathREM$tray) #MJPSS(n,pathREM[[1]][[2]],c(pathREM[[1]][[1]],b))
    MJREM=MJREM+SSREM$matrix/mc
    STREM=STREM+SSREM$vector/mc
  }# for
  MAT<-rbind(MJ,MJU,MJD,MJB,MJRE,MJREM)
  TIEM<-rbind(ST,STU,STD,STB,STRE,STREM)
  return(list(MAT,TIEM))
}


#*****************************************
# Count Sufficient Statistics for all methods
#### ESSREV
ESSREV=function(RESS,a,b,Q,t,mc)
{
  n=dim(Q)[1]
  MJ<-MJU<-MJD<-MJB<-MJRE<-MJREM<-matrix(0,nrow=n,ncol=n)
  ST<-STU<-STD<-STB<-STRE<-STREM<-numeric(n)
  
  for(i in 1:mc)
  {
    # reversed  
    path=bridge_MJP_rev(a,b,Q,t)
    SS=MJPSS(n,path$time,path$tray)
    MJ=MJ+SS$matrix/mc
    ST=ST+SS$vector/mc
    
    #UNIFORMIZATION
    pathU= Uniformization(a,b,Q,t)  
    SSU=MJPSS(n,pathU[[1]][[2]],c(pathU[[1]][[1]],b))
    MJU=MJU+SSU$matrix/mc
    STU=STU+SSU$vector/mc
    
    # DIRECTO     
    pathD=dir_bri(a,b,Q,t)
    SSD=MJPSS(n,pathD$time,pathD$tray)
    MJD=MJD+SSD$matrix/mc
    STD=STD+SSD$vector/mc
    
    #Bisection
    pathB= Bisection(a,b,0,t,Q,itera=0)
    SSB=MJPSS(n,pathB[,2],pathB[,1])
    diag(SSB$matrix)<-0
    MJB=MJB+SSB$matrix/mc
    STB=STB+SSB$vector/mc
    
    #RECHAZO
    pathRE= bridge_MJP_REJ(a,b,Q,t) #Forward.Samp(a,b,Q,t) 
    SSRE=MJPSS(n,pathRE$time,pathRE$tray) #MJPSS(n,pathRE[[1]][[2]],c(pathRE[[1]][[1]],b))
    MJRE=MJRE+SSRE$matrix/mc
    STRE=STRE+SSRE$vector/mc
    
    
    #RECHAZO_MODIFICADO
    pathREM= bridge_MJP_MR(a,b,Q,t) # Mod.Rej.Sam(a,b,Q,t) 
    SSREM=MJPSS(n,pathREM$time,pathREM$tray) #MJPSS(n,pathREM[[1]][[2]],c(pathREM[[1]][[1]],b))
    MJREM=MJREM+SSREM$matrix/mc
    STREM=STREM+SSREM$vector/mc
  }# for
  
  # Reverse  
  Redos<-norm(RESS$NJ-MJ)
  Rtime<-sum(abs(as.matrix(RESS$TS)-ST))
  #Uniformization  
  Uedos<-norm(RESS$NJ-MJU)
  Utime<-sum(abs(as.matrix(RESS$TS)-STU)) 
  #DIRECTO
  Dedos<-norm(RESS$NJ-MJD)
  Dtime<-sum(abs(as.matrix(RESS$TS)-STD)) 
  #BISECTION
  Bedos<-norm(RESS$NJ-MJB)
  Btime<-sum(abs(as.matrix(RESS$TS)-STB))  
  #RECHAZO
  REedos<-norm(RESS$NJ-MJRE)
  REtime<-sum(abs(as.matrix(RESS$TS)-STRE))  
  #RECHAZO MOD
  REMedos<-norm(RESS$NJ-MJREM)
  REMtime<-sum(abs(as.matrix(RESS$TS)-STREM))  
  #todos juntos 
  Estados<-rbind(Redos,Uedos,Dedos,Bedos,REedos,REMedos)
  Tiempos<-rbind(Rtime,Utime,Dtime,Btime,REtime,REMtime)
  
  return(list(Estados,Tiempos))
}
