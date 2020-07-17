# MNLmodelの読み出し
source("./MNLmodel.R") 

# 乱数のseedを設定
set.seed(555)

#パラメータの設定
beta1<- -5
beta2<- 2
beta01<- 2
beta02<- 1
beta03<- 0
betas<-c(beta1,beta2,beta01,beta02,beta03)

hh<-100 #家計数
pt<-20 #購買回数

hhpt<-hh*pt


ID<-matrix(0,hhpt,2) #個人ID
BUY<-matrix(0,hhpt,4) #購買ダミー行列
PRICE<-matrix(0,hhpt,4) #価格
DISP<-matrix(0,hhpt,4) #エンド陳列の有無

for(i in 1:hh){  #i 家計
  for(j in 1:pt){  #j 購買機会
    r<-(i-1)*pt+j 
    ID[r,1]<-i
    ID[r,2]<-j
    
    # ブランド1の販売価格，特別陳列の有無の発生
    rn<-runif(2)
    # 確率0.8で価格は1, 確率0.15で価格は0.9, 確率0.05で価格は0.8
    if (rn[1]<0.8) SP<-1 else
    {if (rn[1]<0.95) SP<-0.9 else SP<-0.8}
    PRICE[r,1]<-SP
    # 確率0.2で特別陳列あり
    DISP[r,1]<-(rn[2]>0.8)
    
    # ブランド2の販売価格，特別陳列の有無の発生
    rn<-runif(2)
    # 確率0.5で価格は1, 確率0.3で価格は0.8, 確率0.2で価格は0.6
    if (rn[1]<0.5) SP<-1 else
    {if (rn[1]<0.8) SP<-0.8 else SP<-0.6}
    PRICE[r,2]<-SP
    # 確率0.1で特別陳列あり
    DISP[r,2]<-(rn[2]>0.9)
    
    # ブランド3の販売価格，特別陳列の有無の発生
    rn<-runif(2)
    # 確率0.7で価格は1, 確率0.1で価格は0.8, 確率0.2で価格は0.6
    if (rn[1]<0.7) SP<-1 else
    {if (rn[1]<0.8) SP<-0.8 else SP<-0.6}
    PRICE[r,3]<-SP
    # 確率0.4で特別陳列あり
    DISP[r,3]<-(rn[2]>0.6)
    
    # ブランド4の販売価格，特別陳列の有無の発生
    rn<-runif(2)
    # 確率0.5で価格は1, 確率0.3で価格は0.8, 確率0.2で価格は0.6
    if (rn[1]<0.5) SP<-1 else
    {if (rn[1]<0.8) SP<-0.8 else SP<-0.6}
    PRICE[r,4]<-SP
    # 確率0.4で特別陳列あり
    DISP[r,4]<-(rn[2]>0.6)
  }
}

#選択確率の計算
PPr<- MNLmodel(betas,PRICE,DISP)

#購買ブランドを決定
for(i in 1:hhpt){
  CSPPr<-cumsum(PPr[i,])  #累積確率を計算
  rn2<-runif(1)           #0‾1の一様分布を発生
  PPM<-which.max(CSPPr>=rn2) #乱数より大きな累積確率の値を持つ対象を選択
  BUY[i,PPM]<- 1 
}

logitdata<-list(betas=betas,hh=hh,pt=pt,ID=ID, PRICE=PRICE,DISP=DISP,BUY=BUY)
save(logitdata,file="logitdata.Rdata")