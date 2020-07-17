# ロジットモデルの関数 --------------------------------------------------------------
MNLmodel <- function(x, PRICE, DISP){
  b1 <- x[1]
  b2 <- x[2]
  b01 <- x[3]
  b02 <- x[4]
  b03 <- x[5]
  
  # 効用の計算
  U1 <- b1*PRICE[, 1]+b2*DISP[, 1]+b01
  U2 <- b1*PRICE[, 2]+b2*DISP[, 2]+b02
  U3 <- b1*PRICE[, 3]+b2*DISP[, 3]+b03
  U4 <- b1*PRICE[, 4]+b2*DISP[, 4]
  
  d <- exp(U1)+exp(U2)+exp(U3)+exp(U4)
  
  # 選択確率
  P1 <- exp(U1)/d
  P2 <- exp(U2)/d
  P3 <- exp(U3)/d
  P4 <- exp(U4)/d
  
  Pr <- cbind(P1,P2,P3,P4)
  return(Pr)
}

