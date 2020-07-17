## 潜在クラスロジットモデル
source("MNLmodel.R")
load("logitdata.Rdata")
set.seed(555)
options(digits=10) #表示桁数を10桁に

## パラメータの設定
beta1<- c(-5, -10, -2)
beta2<- c(2, 3, 1)
beta01<- c(2, 0, 0)
beta02<- c(1, 2, 1)
beta03<- c(0, 3, 2)
BETAs<-cbind(beta1,beta2,beta01,beta02,beta03)

## セグメント別購買データの発生

## 家計数
hh<-logitdata$hh

### 購買機会の回数
pt<-logitdata$pt

hhpt<-hh*pt

ID<-logitdata$ID
PRICE<-logitdata$PRICE
DISP<-logitdata$DISP
BUY<-matrix(0,hhpt,4)
CLS<-matrix(0,hh,1) ## 所属クラス
HTC<-matrix(0,hhpt,hh)

## 所属クラスの決定
rCLS<-runif(hh) 
CLS[(rCLS<0.2),1]<-1
CLS[((rCLS>=0.2) & (rCLS<0.7)),1]<-2
CLS[(rCLS>=0.7),1]<-3

## マーケティング変数と購買データの発生
for(i in 1:hh){
  ### 選択確率の計算
  PP<-MNLmodel(BETAs[CLS[i], ], PRICE[ID[,1]==i, ], DISP[ID[,1]==i,])
  
  for(j in 1:pt){
    r1<-(i-1)*pt+j
    ## 選択ブランドの決定
    CUMSPP<-cumsum(PP[j,])
    rn2<-runif(1)
    PPM<-which.max(CUMSPP>=rn2)
    BUY[r1,PPM]<- 1
  }
}


### 潜在クラスロジットモデルによる推定

cll <- function(x) {   ## 完全データのロジットモデル部分の尤度
  b1 <- x[1:3]
  b2 <- x[4:6]
  b01<- x[7:9]
  b02<- x[10:12]
  b03<- x[13:15]
  
  ones<-rep(1,hhpt)
  
  # 完全データでの購買機会×セグメント別の尤度を計算して合計
  U1<- outer(PRICE[,1],b1) + outer(DISP[,1],b2)+outer(ones,b01)
  U2<- outer(PRICE[,2],b1) + outer(DISP[,2],b2)+outer(ones,b02)
  U3<- outer(PRICE[,3],b1) + outer(DISP[,3],b2)+outer(ones,b03)
  U4<- outer(PRICE[,4],b1) + outer(DISP[,4],b2)
  d<- exp(U1)+exp(U2)+exp(U3)+exp(U4)
  LLC<-(matrix(BUY[,1],hhpt,3)*U1+matrix(BUY[,2],hhpt,3)*U2
        +matrix(BUY[,3],hhpt,3)*U3+matrix(BUY[,4],hhpt,3)*U4-log(d))
  LL<-sum(zpt * LLC)
  return(LL)
}


ollz <- function(x) {   ## 観測データでの尤度とzの計算
  b1 <- x[1:3]
  b2 <- x[4:6]
  b01<- x[7:9]
  b02<- x[10:12]
  b03<- x[13:15]
  r  <- x[16:18]
  
  LCo<- matrix(0,hhpt,3)
  LLho<- matrix(0,hh,3)
  
  for(s in 1:3){
    betas<-c(b1[s], b2[s], b01[s], b02[s], b03[s])
    PP<-MNLmodel(betas, PRICE, DISP)
    LCo[,s]<-apply(PP^BUY,1,prod)
  }
  
  for(i in 1:hh){
    LLho[i,]<-apply(LCo[ID[,1]==i,],2,prod)
  }
    
  #観測データでの対数尤度	
  LLo<-sum(log(apply(matrix(r,hh,3,byrow=T) * LLho,1,sum)))
  
  z0<-matrix(r,hh,3,byrow=T) * LLho
  z1<-z0 / matrix(apply(z0,1,sum),hh,3)
  
  rval<-list(LLo=LLo,z1=z1)
  return(rval)
}

## EM アルゴリズムの初期値
iter<-0
beta<-c(-5,-10,-2, rep(0,12))
r<-c(0.3, 0.3, 0.4)
obsllz<-ollz(c(beta,r))
LL1<-obsllz$LLo
z<- obsllz$z1
zpt<-matrix(0,hhpt,3)

dl <- 100 # EMステップでの対数尤度の差の初期値を設定
tol<- 10^(-10)


## EM アルゴリズムによる推定
while (abs(dl) >= tol) {  ### dlがtol以上の場合は繰り返す
  
  for(i in 1:hhpt){
    zpt[i,]<-z[ID[i,1],]
  }
  
  #完全データでのロジットモデルの推定    
  res<-optim(beta, cll, method = "BFGS", hessian = TRUE, control=list(fnscale=-1))
  beta<-res$par
  r<- apply(z,2,sum)/hh
  
  obsllz<-ollz(c(beta,r))
  LL<-obsllz$LLo
  z<- obsllz$z1
  
  iter<-iter+1
  dl<-LL-LL1
  print(LL)
  LL1<-LL
}

options(digits=7)
B<-matrix(beta,5,3,byrow=TRUE)
B
r