library(quantmod)
library(quadprog)
library(openxlsx)
library(readxl)
library(lpSolve)
library(Rglpk)
library(lubridate)


#Pre potreby spustenia tohto testu musime najskor spustit Model1_FSD 
#################################################################################################

#Vstupne parametre (mozu sa menit pouzivatelom)
number_of_shares <- 25
pocet_mesiacov <- 25
step <- 20 #Pri dennych datach treba zadat 1, pri mesacnych 20
start_date <- as.Date("2018-01-01")
end_date <- as.Date("2024-04-01") %m-% months(pocet_mesiacov)
epsilon <- 0.01

#Precitanie dat o firmach z S&P 500 
current_dir <- getwd()
setwd(current_dir)
excel_file <- "S&P_500_Company_List.xlsx"
Companies <- read_excel(excel_file)

#Definovanie pocitadiel 
Pocet <- NULL
Pocet_FSD <- NULL
Cas_FSD <- NULL
Info <- NULL
Pocet_dni <- NULL

#Testovanie 
for (j in 1:10){
  #Vyber nahodnych akcii 
  rand <- sample(1:nrow(Companies), number_of_shares, replace = FALSE)
  Symbols <- Companies[rand,1]
  Symbols <- as.vector(Symbols[[1]])
  
  #Vyber nahodneho casoveho obdobia
  from <- sample(seq(start_date, end_date, by = "day"), 1)
  to <- from %m+% months(pocet_mesiacov)
  from <- as.character(from)
  to <- as.character(to)
  
  
  #Ukladanie informacii o nahodnych vyberoch
  A <- c(from, to, Symbols)
  Info <- cbind(Info, A)
  
  
  #Nacitanie dat o cenach akcii z Yahoo finance
  Data <- new.env()
  getSymbols(Symbols, from = from, to = to, env = Data)
  
  #Prevedenie cien akcii do formatu aby sme s nim vedeli pracovat 
  Akcie <- NULL
  for (i in Symbols){
    Akcie <- cbind(Akcie, Data[[i]][,6])
  }
  Akcie <- as.matrix(Akcie) 
  
  #Casova perioda
  casova_perioda <- seq(1, nrow(Akcie), by = step)
  Akcie <- Akcie[casova_perioda,]
  
  #Nemali by tam byt 0 ale pre istotu to tu spustim (nuly by nemali byt lebo akcie su vsetky z rovnakeho trhu, v rovnake dni je otvoreny a zatvoreny)
  ind <- which(rowSums(is.na(Akcie)) == 0)
  Akcie <- Akcie[ind,]
  
  Akcie_vynosy <- Akcie[2:nrow(Akcie),]/Akcie[1:(nrow(Akcie)-1),] - 1
  
  Pocet_dni <- c(Pocet_dni, nrow(Akcie_vynosy))
  
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
  y <- seq(from = max(r)-epsilon, to = min(r)+epsilon, by = -epsilon)
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
  
  Pocet <- c(Pocet, ncol(new_W))
  
  #########################################################################################################
  #########################################################################################################
  #Model 1 FSD test
  test1_FSD <- Model1_FSD(index, new_y, Akcie_vynosy, new_W)
  test1_FSD_cas <- test1_FSD$cas
  Pocet_FSD <- c(Pocet_FSD, length(test1_FSD$neefektivna_mnozina))
  Cas_FSD <- c(Cas_FSD, round(test1_FSD_cas, digits = 4))
  
  #Vypisovanie poctu iteracii a uspatie systemu na 7 sekund kvoli eliminacii preplnenia pamate
  print(j)
  Sys.sleep(7)
}
rm(Akcie, Akcie_vynosy, Companies, Data, new_W, sol_numericke, V, W, test1_FSD, Port_vynosy)
Info <- rbind(Info, Pocet_dni, Pocet, Pocet_FSD, Cas_FSD)

Pocet <- as.matrix(Pocet)
Pocet_FSD <- as.matrix(Pocet_FSD)



