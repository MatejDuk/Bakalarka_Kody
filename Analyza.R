library(quantmod)
library(quadprog)
library(openxlsx)
library(readxl)
library(lpSolve)
library(Rglpk)


#Pre potreby spustenia tohto skriptu musime najskor spustit Model2_Simplex a Model3_IPM a Model3_Simplex
###################################################################################################
#Nami zvolene akcie
Symbols <- c("VRTX", "AWK", "DD", "AON", "ON", "EMR", "SWK", "OMC", "AME", "DHI", "IDXX", "WTW", "OKE", "EMN", "PEG", "VTRS", "BF-B", "GWW", "LDOS", "FI", "AEP", "JKHY", "CSX", "WFC", "REGN")

#Nacitanie dat o cenach akcii z Yahoo finance
Data <- new.env()
getSymbols(Symbols, from = "2021-11-13", to ="2023-12-13", env = Data)

#Prevedenie cien akcii do formatu aby sme s nim vedeli pracovat 
Akcie <- NULL
for (i in Symbols){
  Akcie <- cbind(Akcie, Data[[i]][,6])
}
Akcie <- as.matrix(Akcie) 

#Casova perioda
casova_perioda <- seq(1, nrow(Akcie), by = 20)
Akcie <- Akcie[casova_perioda,]

#Nemali by tam byt 0 ale pre istotu to tu spustim (nuly by nemali byt lebo akcie su vsetky z rovnakeho trhu, v rovnake dni je otvoreny a zatvoreny)
ind <- which(rowSums(is.na(Akcie)) == 0)
Akcie <- Akcie[ind,]

Akcie_vynosy <- Akcie[2:nrow(Akcie),]/Akcie[1:(nrow(Akcie)-1),] - 1


#Prepocet vynosov na rocny vynos
r <- apply(X = Akcie_vynosy  ,
           MARGIN = 2 ,
           FUN = mean)
r <- r * 252/20

#Vypocet kovariancnej matice s rocnymi udajmi
V <- cov(Akcie_vynosy)
V <- V * 252/20

############################################################################################################################################################
#Markowitz
y <- seq(from = max(r)-0.01, to = 0, by = -0.01)
x <- rep(0, length(y))
n <- length(Symbols)
W <- NULL

for (i in 1:length(y)){
  r_p <- y[i]
  sol_numericke <- solve.QP(Dmat =  V,
                            dvec = rep(0, n) ,
                            Amat =  cbind(r, rep(1,n), diag(n)),
                            bvec =  c(r_p, 1, rep(0,n)),
                            meq =   2)
  
  W <-cbind(W, sol_numericke$solution)
  
  x[i] <- sqrt(sol_numericke$solution %*% V %*% sol_numericke$solution)
}

#Vykreslenie efektivnej hranice, pre lepsiu  predstavu 
plot(x = x, y = y,col = "red", type = "l", xlab = "sigma_p", ylab = "r_p")

#Vykreslenie iba portfolii, ktore su relevantne pre vyber 
index <- which.min(x)
new_x <- x[1:index]
new_y <- y[1:index]
new_W <- W[,1:index]

lines(x = new_x, y=new_y, col='black', type = 'l')
Port_vynosy <- as.matrix(Akcie_vynosy %*% new_W)


#########################################################################################################
#Model2_Simplex - aj sa vykreslia, ktore portfolia z efektivnej Markowitzovej hranice su SSD efektivne a neefektivne

neefektivne <- NULL
hodnoty <- NULL
test2_cas <- NULL
for (i in 1:ncol(new_W)){
  R <- Akcie_vynosy
  ref <- as.matrix(Port_vynosy[,i])
  test2 <- Model2_Simplex(R,ref)
  if (round(test2$optimalna_hodnota, digits = 10)){
    neefektivne <- c(neefektivne, i)
  }
  hodnoty <- c(hodnoty, test2$optimalna_hodnota)
  test2_cas <- c(test2_cas, test2$cas)
}
if (length(neefektivne)>0){
  neefektivne <- as.matrix(neefektivne)
}

for (k in 1:ncol(new_W)) {
  if(k %in% neefektivne){
    lines(new_x[k], new_y[k], type = "p", pch = 19, col = "red", cex = 2.2)
  }
  else{
    lines(new_x[k], new_y[k], type = "p", pch = 19, col = "green", cex = 2.2)
  }
}
legend("bottomright", legend = c("Markowitzova efektívna hranica", "SSD efektívne", "SSD neefektívne"), 
       col = c("black", "green", "red"), 
       lty = c(1, NA, NA),  
       pch = c(NA, 19, 19),  
       pt.cex = c(1, 2, 2),  
       title = "Legenda")

#########################################################################################################
#Model3_IPM - pre kazde neefektivne portfolio sa najde to ktore mu dominuje na zaklade modelu3, aj sa vykresli graf s povodnym portfoliom a k nemu dominantnym
plot(x = x[8:31], y = y[8:31],col = "red", type = "l", xlab = "sigma_p", ylab = "r_p")
lines(x = new_x, y=new_y, col='black', type = 'l')
for(j in 1:length(neefektivne)){
  cp <- neefektivne[j]
  ref <- as.matrix(Port_vynosy[,cp])
  lines(new_x[cp],new_y[cp], type = "p", pch =19, col = "red", cex = 3.4)
  text(new_x[cp],new_y[cp], labels = j, cex = 1.2, font = 2)
  R <- Akcie_vynosy
  vysledok <- Model3_Simplex(R,ref)
  vys <- Model3_IPM(R,ref)
  rs <- vys$optimalne_riesenie[1:n]
  x1 <- sqrt(rs %*% V %*% rs)
  y1 <- r %*% rs
  lines(x1, y1, type = "p",pch = 19, col = "blue", cex = 3.4)
  text(x1,y1, labels = j, cex = 1.2, col = "white", font = 2)
}
legend("bottomright", legend = c("Markowitzova efektívna hranica", "SSD neefektívne", "Dominantné portfólio"), 
       col = c("black", "red", "blue"), 
       lty = c(1, NA, NA),  
       pch = c(NA, 19, 19), 
       pt.cex = c(1, 2, 2),  
       title = "Legenda")
