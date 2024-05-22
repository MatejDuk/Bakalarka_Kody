library(quantmod)
library(quadprog)
library(openxlsx)
library(readxl)
library(lpSolve)
library(Rglpk)

#Pre potreby spustenia tohto skriptu musime najskor spustit  Model3_div_IPM
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
y <- seq(from = 0.195, to = 0.06, by = -0.005)
x <- rep(0, length(y))
n <- length(Symbols)
W <- NULL

for (i in 1:length(y)){
  r_p <- y[i]
  sol_numericke <- solve.QP(Dmat =  V,
                            dvec = rep(0, n) ,
                            Amat =  cbind(r, rep(1,n), diag(n), -diag(n)),
                            bvec =  c(r_p, 1, rep(0,n), -rep(0.1,n)),
                            meq =   2)
  
  W <-cbind(W, sol_numericke$solution)
  
  x[i] <- sqrt(sol_numericke$solution %*% V %*% sol_numericke$solution)
}

#Vykreslenie efektivnej hranice, pre lepsiu  predstavu 
plot(x = x[1:(length(x)-4)], y = y[1:(length(y)-4)],col = "red", type = "l", xlab = "sigma_p", ylab = "r_p")

#Vykreslenie iba portfolii, ktore su relevantne pre vyber 
index <- which.min(x)
new_x <- x[1:index]
new_y <- y[1:index]
new_W <- W[,1:index]

lines(x = new_x, y=new_y, col='black', type = 'l')
Port_vynosy <- as.matrix(Akcie_vynosy %*% new_W)

#########################################################################################################
#Model3_div_IPM 
#Ku kazdemu portfoliu z Markowitzovej efektivnej hranice sa najde domianntne, obe sa vykreslia
#K tomu sa aj vypocita chyba potrebna na zistenie, ci sa jedna o efektivne alebo neefektivne portfolio
chyba <- NULL
for (l in 1:ncol(new_W)){
  ref <- as.matrix(Port_vynosy[,l])
  R <- Akcie_vynosy
  model <- Model3_div_IPM(R,ref)
  w_povodna <- new_W[,l]
  w_nova <- model$optimalne_riesenie[1:n,]
  chyba <- c(chyba, norm(w_povodna-w_nova, type = "2"))
  lines(new_x[l], new_y[l], type = "p", pch =19, col = "red", cex = 3.4)
  text(new_x[l], new_y[l], labels = l, cex = 1.2, font = 2)
  x1 <- sqrt(w_nova %*% V %*% w_nova)
  y1 <- r %*% w_nova
  lines(x1,y1, new_y[l], type = "p", pch =19, col = "blue", cex = 3.4)
  text(x1,y1, labels = l, cex = 1.2, col = "white", font = 2)
}
legend("bottomright", legend = c("Markowitzova efektívna hranica", "Pôvodné portfólio", "Neostro dominantné portfólio"), 
       col = c("black", "red", "blue"), 
       lty = c(1, NA, NA),  # NA for no line (only for points)
       pch = c(NA, 19, 19),  # NA for no points (only for lines)
       pt.cex = c(1, 2, 2),  # Point size
       title = "Legenda")

chyba <- as.matrix(round(chyba,16))
