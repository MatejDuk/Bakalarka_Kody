library(lpSolve)
library(Rglpk)

#V tomto skripte si postupne definujeme vsetky modely z bakalarskej prace ako funkcie
##############################################################################################################################
#Model 1 FSD, tak ako je algoritmus naznaceny v bakalarke

Model1_FSD <- function(index, new_y, Akcie_vynosy, new_W){
  #Vstupne parametre:
  #index = pocet pozorovanych portfolii, v nasom pripade pocet portfolii, ktore su na Markowitzovej efektivnej hranici
  #new_y = ocakavany vynos (os y) pre kazde portfolio z efeektivnej Markowitzovej hranice
  #Akcie_vynosy = vynosy jednotlivych akcii, v bakalarke matica R
  #new_W = vaha pre jednotlive portfolia z Markowitzovej hranice
  
  #Vytvorenie vsetkych moznych dvojic
  st <- Sys.time()
  dvojice <- expand.grid(c(1:index), c(1:index))
  dvojice <- dvojice[dvojice[, 1] !=  dvojice[, 2],]
  r_port <- as.matrix(new_y)
  Port_vynosy <- as.matrix(Akcie_vynosy %*% new_W)
  
  
  #Redukcia 1 z bakalarky 
  red1 <- dvojice[dvojice[,1] < dvojice[,2],]
  
  #Redukcia 2 z bakalarky 
  min_Port_vynosy <- as.matrix(apply(Port_vynosy, 2, min))
  cisla <- NULL
  
  for (i in 1:nrow(red1)){
    a <- red1[i,1]
    b <- red1[i,2]
    if(min_Port_vynosy[a] >= min_Port_vynosy[b]){
      cisla <- cbind(cisla, i)
    }
  }
  red2 <-red1[cisla,]
  
  
  #Zoradenie portfolii od najvacsieho vynosu po najmensi
  zor_Port_vynosy <- apply(Port_vynosy, 2, sort)
  n <- nrow(red1)
  m <- nrow(zor_Port_vynosy)
  
  #Test na FSD
  neefektivne <- NULL
  
  #Podmienka na velkost red2
    for (i in 1:n){
      a <- red1[i,1]
      b <- red1[i,2]
      j <- 1
      test <- 0 
      while (zor_Port_vynosy[j,a] >= zor_Port_vynosy[j,b]){
        if (zor_Port_vynosy[j,a] > zor_Port_vynosy[j,b]){
          test <- 1
        }
        j <- j + 1
        if (j-1 == m){
          break
        }
      }
      if(j-1 == m & test == 1){
        neefektivne <- cbind(neefektivne, b)
      }
    }
    neefektivne <- as.vector(neefektivne)
    if (length(neefektivne) != 0){
      neefektivne <- as.matrix(neefektivne[!duplicated(neefektivne)])
    }
  et <- Sys.time()
  cas <- et - st
  
  return(list(dvojice = dvojice, redukcia1 = red1, redukcia2 = red2, neefektivna_mnozina = neefektivne, cas = cas))
  #Vystupy z funkcie su vsetky dvojice, dvojice po redukcii 1, dvojice po redukcii 2, cisla portfolii, ktore sa zahrnu do neefektivnej mnoziny a cas, ktory trva posudenie vsetkych portfolii
}

##############################################################################################################################
#Model 1 SSD

Model1_SSD <- function(index, new_y, Akcie_vynosy, new_W){
  #Vstupne parametre:
  #index = pocet pozorovanych portfolii, v nasom pripade pocet portfolii, ktore su na Markowitzovej efektivnej hranici
  #new_y = ocakavany vynos (os y) pre kazde portfolio z efeektivnej Markowitzovej hranice
  #Akcie_vynosy = vynosy jednotlivych akcii, v bakalarke matica R
  #new_W = vaha pre jednotlive portfolia z Markowitzovej hranice
  
  #Vytvorenie vsetkych moznich dvojic
  st <- Sys.time()
  dvojice <- expand.grid(c(1:index), c(1:index))
  dvojice <- dvojice[dvojice[, 1] !=  dvojice[, 2],]
  r_port <- as.matrix(new_y)
  Port_vynosy <- as.matrix(Akcie_vynosy %*% new_W)
  
  
  #Redukcia 1
  red1 <- dvojice[dvojice[,1] < dvojice[,2],]
  
  #Redukcia 2
  min_Port_vynosy <- as.matrix(apply(Port_vynosy, 2, min))
  cisla <- NULL
  
  for (i in 1:nrow(red1)){
    a <- red1[i,1]
    b <- red1[i,2]
    if(min_Port_vynosy[a] >= min_Port_vynosy[b]){
      cisla <- cbind(cisla, i)
    }
  }
  red2 <-red1[cisla,]
  
  #Zoradenie portfolii od najvacsieho vynosu po najmensi
  zor_Port_vynosy <- apply(Port_vynosy, 2, sort)
  n <- nrow(red1)
  m <- nrow(zor_Port_vynosy)
  
  #SSD
  cs <- apply(zor_Port_vynosy, 2, cumsum)
  neefektivne_SSD <- NULL
  
  for (i in 1:n){
    a <- red1[i,1]
    b <- red1[i,2]
    j <- 1
    test <- 0 
    while (cs[j,a] >= cs[j,b]){
      if (cs[j,a] > cs[j,b]){
        test <- 1
      }
      j <- j + 1
      if (j-1 == m){
        break
      }
    }
    if(j-1 == m & test == 1){
      neefektivne_SSD <- cbind(neefektivne_SSD, b)
    }
  }
  neefektivne_SSD <- as.vector(neefektivne_SSD)
  if (length(neefektivne_SSD) != 0){
    neefektivne_SSD <- as.matrix(neefektivne_SSD[!duplicated(neefektivne_SSD)])
  }
  
  et <- Sys.time()
  cas <- et - st
  
  return(list(dvojice = dvojice, redukcia1 = red1, redukcia2 = red2, neefektivna_mnozina = neefektivne_SSD, cas = cas))
  #Vystupy z funkcie su vsetky dvojice, dvojice po redukcii 1, dvojice po redukcii 2, cisla portfolii, ktore sa zahrnu do neefektivnej mnoziny a cas, ktory trva posudenie vsetkych portfolii
}

##############################################################################################################################
#Model 2 - Postov test - pri rieseni ulohy linearneho programovania vyuzita Simplexova metoda

Model2_Simplex <- function(R, ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  st <- Sys.time()
  #Usporiadanie vynosov referencneho aktiva vzostupne 
  ref <- apply(ref, 2, sort)
  
  #Usporiadanie vynosov vsetkych akcii podla referencneho aktiva, podla bakalarky je to vlastne matica X
  desired_order <- rownames(ref)
  R <- R[desired_order, ,drop=FALSE]
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Pri odvolavani sa na bakalarku myslime ulohu (4)
  
  #Prva nerovnica z bakalarky
  a <- as.vector(ref)
  a <- t(replicate(n,  a))-t(R)
  b <- as.matrix(rep(m,n))
  
  A <- cbind(b,a)
  rm(a,b)
  
  #Druha nerovnica z bakalarky
  a <- NULL
  for (i in 1:(m-1)){
    aa <- c(rep(0,i),1,-1,rep(0,(m-i-1)))
    a <- rbind(a, aa)
  }
  A <- rbind(A, a)
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m-1,ncol = 1), diag(m-1), matrix(0, nrow = m-1,ncol = 1))
  A <- rbind(A,a)           
  
  #Prva rovnica z bakalarky
  a <- c(rep(0,m),1)
  A <- rbind(A,a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  B <- c(rep(0,(m-1)+(m-1)+n),1)
  
  #Definicia ucelovej funkcie
  f <- c(1, rep(0,m))
  
  #Pocitanie Simplexovou metodou
  vysledok <- lp("min", f, A, c(rep(">=",(m-1)+(m-1)+n), "=="), B)
  
  #Prvy vystup = optimalna hodnota ucelovej funkcie
  optimalna_hodnota <- vysledok$objval
  et <- Sys.time()
  
  #Druhy vystup = optimalne riesenie
  optimalne_riesenie <- vysledok$solution
  
  #Treti vystup = dlzka trvania jedneho Postovho testu
  cas <- et - st
  
  #Stvrty vystup = existencia riesenia
  existencia_riesenia <- NULL
  if (vysledok$status == 0){
    existencia_riesenia <- "Solution found."
  }
  else{
    existencia_riesenia <- "Error occured during solving."
  }
  
  return(list(optimalna_hodnota = optimalna_hodnota, optimalne_riesenie = optimalne_riesenie, cas = cas, existencia_riesenia = existencia_riesenia))
}

##############################################################################################################################
#Model 2 - Postov test - pri rieseni ulohy linearneho programovania vyuzita metoda vnutorneho bodu

Model2_IPM <- function(R, ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  st <- Sys.time()
  #Usporiadanie vynosov referencneho aktiva vzostupne 
  ref <- apply(ref, 2, sort)
  
  #Usporiadanie vynosov vsetkych akcii podla referencneho aktiva
  desired_order <- rownames(ref)
  R <- R[desired_order, ,drop=FALSE]
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Pri odvolavani na bakalrku myslime ulohu (4)
  
  #Prva nerovnica z bakalarky
  a <- as.vector(ref)
  a <- t(replicate(n,  a))-t(R)
  b <- as.matrix(rep(m,n))
  
  A <- cbind(b,a)
  rm(a,b)
  
  #Druha nerovnica z bakalarky
  a <- NULL
  for (i in 1:(m-1)){
    aa <- c(rep(0,i),1,-1,rep(0,(m-i-1)))
    a <- rbind(a, aa)
  }
  A <- rbind(A, a)
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m-1,ncol = 1), diag(m-1), matrix(0, nrow = m-1,ncol = 1))
  A <- rbind(A,a)           
  
  #Prva rovnica z bakalarky
  a <- c(rep(0,m),1)
  A <- rbind(A,a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  B <- c(rep(0,(m-1)+(m-1)+n),1)
  
  #Definicia ucelovej funkcie
  f <- c(1, rep(0,m))
  
  #Pocitanie Simplexovou metodou
  vysledok <- Rglpk_solve_LP(f, A, c(rep(">=",(m-1)+(m-1)+n), "=="), B, max = FALSE, method = "interior")
  
  #Prvy vystup = optimalna hodnota ucelovej funkcie
  optimalna_hodnota <- vysledok$optimum
  et <- Sys.time()
  
  #Druhy vystup = optimalne riesenie
  optimalne_riesenie <- vysledok$solution
  
  #Treti vystup = dlzka trvania jedneho Postovho testu
  cas <- et - st
  
  #Stvrty vystup = existencia riesenia
  existencia_riesenia <- NULL
  if (vysledok$status == 0){
    existencia_riesenia <- "Solution found."
  }
  else{
    existencia_riesenia <- "Error occured during solving."
  }
  
  return(list(optimalna_hodnota = optimalna_hodnota, optimalne_riesenie = optimalne_riesenie, cas = cas, existencia_riesenia = existencia_riesenia))
}

##############################################################################################################################
#Model 3 - Dentcheva-Ruszczynsky test - pri rieseni ulohy linearneho programovania vyuzita Simplexova metoda

Model3_Simplex <- function(R, ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Definovanie pravdepodobnosti vyskytu daneho vynosu pri jednotlivych akciach
  p <- 1/m
  
  #Rozmer nasej neznámej bude m*m+n
  A <- NULL
  
  #Pri odvolavani na bakalrku myslime ulohu (7)
  
  #Prva nerovnica z ulohy
  for (i in 1:m){
    a <- c(rep(0,n), rep(0,m*(i-1)), rep(p,m), rep(0,m*(m-i)))
    A <- rbind(A,a)
  }
  
  
  #Druha nerovnica z bakalarky
  rep_a <- do.call(rbind, replicate(m, -R, simplify = FALSE))
  a <- cbind(rep_a, -diag(m*m))
  A <- rbind(A,a)
  
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m*m, ncol = n), -diag(m*m))
  A <- rbind(A, a)
  
  #Stvrta nerovnica z bakalarky 
  a <- cbind(-diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A,a)
  
  #Rovnica z bakalarky 
  a <- c(rep(1,n), rep(0,(m*m)))
  A <- rbind(A,a)
  rm(a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  
  bb <- NULL
  B <- NULL
  for(k in 1:m){
    for (j in 1:m){
      b <- max(ref[k,1]-ref[j,1],0)
      bb <- cbind(bb,b)
    }
    B <- cbind(B,bb%*%rep(p,m))
    bb<- NULL
  }
  
  for (k in 1:m){
    b <- c(-rep(ref[k,1],m))
    B <- c(B,b)
  }
  
  B <- c(B,rep(0,(m*m)))
  B <- c(B,rep(0,n))
  B <- c(B,1)
  
  f <- apply(X = R  ,
             MARGIN = 2 ,
             FUN = mean)
  f <- c(f,rep(0,m*m))
  
  #Pocitanie simplexovou metodou
  st <- Sys.time()
  vysledok_simplex <- lp("max", f, A, c(rep("<=",2*m^2+m+n), "=="), B)
  
  #Prvy vystup
  optimalne_riesenie <- vysledok_simplex$solution
  optimalne_riesenie <- as.matrix(optimalne_riesenie)
  et <- Sys.time()
  
  #Druhy vystup
  cas <- et - st
  
  #Treti vystup
  riesenie <- NULL
  if (vysledok_simplex$status == 0){
    riesenie <- "Solution found."
  }
  else{
    riesenie <- "Error occured during solving."
  }
  
  #Stvrty vystup 
  optimalna_hodnota <- vysledok_simplex$objval
  
  return(list(optimalne_riesenie =optimalne_riesenie, cas = cas, riesenie = riesenie, optimalne_hodnota = optimalna_hodnota))
}

##############################################################################################################################
#Model 3 - Dentcheva-Ruszczynsky test - pri rieseni ulohy linearneho programovania vyuzita metoda vnutorneho bodu

Model3_IPM <- function(R,ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Definovanie pravdepodobnosti vyskytu daneho vynosu pri jednotlivych akciach
  p <- 1/m
  
  #Rozmer nasej neznámej bude m*m+n
  A <- NULL
  
  #Pri odvolavani na bakalraku myslime ulohu (7)
  
  #Prva nerovnica z ulohy
  for (i in 1:m){
    a <- c(rep(0,n), rep(0,m*(i-1)), rep(p,m), rep(0,m*(m-i)))
    A <- rbind(A,a)
  }
  
  
  #Druha nerovnica z bakalarky
  rep_a <- do.call(rbind, replicate(m, -R, simplify = FALSE))
  a <- cbind(rep_a, -diag(m*m))
  A <- rbind(A,a)
  
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m*m, ncol = n), -diag(m*m))
  A <- rbind(A, a)
  
  #Stvrta nerovnica z bakalarky 
  a <- cbind(-diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A,a)
  
  #Rovnica z bakalarky 
  a <- c(rep(1,n), rep(0,(m*m)))
  A <- rbind(A,a)
  rm(a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  
  bb <- NULL
  B <- NULL
  for(k in 1:m){
    for (j in 1:m){
      b <- max(ref[k,1]-ref[j,1],0)
      bb <- cbind(bb,b)
    }
    B <- cbind(B,bb%*%rep(p,m))
    bb<- NULL
  }
  
  for (k in 1:m){
    b <- c(-rep(ref[k,1],m))
    B <- c(B,b)
  }
  
  B <- c(B,rep(0,(m*m)))
  B <- c(B,rep(0,n))
  B <- c(B,1)
  
  f <- apply(X = R  ,
             MARGIN = 2 ,
             FUN = mean)
  f <- c(f,rep(0,m*m))
  
  #Pocitanie simplexovou metodou
  st <- Sys.time()
  vysledok_IPM <- Rglpk_solve_LP(f, A, c(rep("<=",2*m^2+m+n), "=="), B, max = TRUE, method = "interior")
  
  #Prvy vystup
  optimalne_riesenie <- vysledok_IPM$solution
  optimalne_riesenie <- as.matrix(optimalne_riesenie)
  et <- Sys.time()
  
  #Druhy vystup
  cas <- et - st
  
  #Treti vystup
  riesenie <- NULL
  if (vysledok_IPM$status == 0){
    riesenie <- "Solution found."
  }
  else{
    riesenie <- "Error occured during solving."
  }
  
  #Stvrty vystup 
  optimalna_hodnota <- vysledok_IPM$optimum
  
  return(list(optimalne_riesenie =optimalne_riesenie, cas = cas, riesenie = riesenie, optimalne_hodnota = optimalna_hodnota))
}

##############################################################################################################################
#Model 3 - upraveny pre ohranicenie maximalnych vah (vacsiu diverzifikaciu) 
#Pri rieseni ulohy linearneho programovania vyuzita Simplexova metoda

Model3_div_Simplex <- function(R, ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Definovanie pravdepodobnosti vyskytu daneho vynosu pri jednotlivych akciach
  p <- 1/m
  
  #Rozmer nasej neznámej bude m*m+n
  A <- NULL
  #Prva nerovnica z ulohy
  for (i in 1:m){
    a <- c(rep(0,n), rep(0,m*(i-1)), rep(p,m), rep(0,m*(m-i)))
    A <- rbind(A,a)
  }
  
  
  #Druha nerovnica z bakalarky
  rep_a <- do.call(rbind, replicate(m, -R, simplify = FALSE))
  a <- cbind(rep_a, -diag(m*m))
  A <- rbind(A,a)
  
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m*m, ncol = n), -diag(m*m))
  A <- rbind(A, a)
  
  #Stvrta nerovnica z bakalarky 
  a <- cbind(-diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A,a)
  
  #Piata nerovnica z bakalarky = toto je cast, kde sa ohranicuje maximalne vaha jednotlivych aktiv 
  a <- cbind(diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A, a)
  
  #Rovnica z bakalarky 
  a <- c(rep(1,n), rep(0,(m*m)))
  A <- rbind(A,a)
  rm(a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  
  bb <- NULL
  B <- NULL
  for(k in 1:m){
    for (j in 1:m){
      b <- max(ref[k,1]-ref[j,1],0)
      bb <- cbind(bb,b)
    }
    B <- cbind(B,bb%*%rep(p,m))
    bb<- NULL
  }
  
  for (k in 1:m){
    b <- c(-rep(ref[k,1],m))
    B <- c(B,b)
  }
  
  B <- c(B,rep(0,(m*m)))
  B <- c(B,rep(0,n))
  B <- c(B,rep(0.1,n))
  B <- c(B,1)
  
  f <- apply(X = R  ,
             MARGIN = 2 ,
             FUN = mean)
  f <- c(f,rep(0,m*m))
  
  #Pocitanie simplexovou metodou
  st <- Sys.time()
  vysledok_simplex <- lp("max", f, A, c(rep("<=",2*m^2+m+n+n), "=="), B)
  
  #Prvy vystup
  optimalne_riesenie <- vysledok_simplex$solution
  optimalne_riesenie <- as.matrix(optimalne_riesenie)
  et <- Sys.time()
  
  #Druhy vystup
  cas <- et - st
  
  #Treti vystup
  riesenie <- NULL
  if (vysledok_simplex$status == 0){
    riesenie <- "Solution found."
  }
  else{
    riesenie <- "Error occured during solving."
  }
  
  #Stvrty vystup 
  optimalna_hodnota <- vysledok_simplex$objval
  
  return(list(optimalne_riesenie =optimalne_riesenie, cas = cas, riesenie = riesenie, optimalne_hodnota = optimalna_hodnota))
}

##############################################################################################################################
#Model 3 - upraveny pre ohranicenie maximalnych vah (vacsiu diverzifikaciu) 
#Pri rieseni ulohy linearneho programovania vyuzita metoda vnutorneho bodu

Model3_div_IPM <- function(R,ref){
  #Vstupy do tejto funkcie:
  #R = matica vynosov jednotlivych aktiv
  #ref = vynos refernecneho portfolia
  
  #Definovanie rozmerov
  m <- nrow(R)
  n <- ncol(R)
  
  #Definovanie pravdepodobnosti vyskytu daneho vynosu pri jednotlivych akciach
  p <- 1/m
  
  #Rozmer nasej neznámej bude m*m+n
  A <- NULL
  #Prva nerovnica z ulohy
  for (i in 1:m){
    a <- c(rep(0,n), rep(0,m*(i-1)), rep(p,m), rep(0,m*(m-i)))
    A <- rbind(A,a)
  }
  
  
  #Druha nerovnica z bakalarky
  rep_a <- do.call(rbind, replicate(m, -R, simplify = FALSE))
  a <- cbind(rep_a, -diag(m*m))
  A <- rbind(A,a)
  
  
  #Tretia nerovnica z bakalarky
  a <- cbind(matrix(0, nrow = m*m, ncol = n), -diag(m*m))
  A <- rbind(A, a)
  
  #Stvrta nerovnica z bakalarky 
  a <- cbind(-diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A,a)
  
  #Piata nerovnica z bakalarky = toto je cast, kde sa ohranicuje maximalne vaha jednotlivych aktiv
  a <- cbind(diag(n), matrix(0, nrow = n, ncol = m*m))
  A <- rbind(A, a)
  
  #Rovnica z bakalarky 
  a <- c(rep(1,n), rep(0,(m*m)))
  A <- rbind(A,a)
  rm(a)
  
  
  #Nasleduje definicia b (vektor na pravej strane)
  
  bb <- NULL
  B <- NULL
  for(k in 1:m){
    for (j in 1:m){
      b <- max(ref[k,1]-ref[j,1],0)
      bb <- cbind(bb,b)
    }
    B <- cbind(B,bb%*%rep(p,m))
    bb<- NULL
  }
  
  for (k in 1:m){
    b <- c(-rep(ref[k,1],m))
    B <- c(B,b)
  }
  
  B <- c(B,rep(0,(m*m)))
  B <- c(B,rep(0,n))
  B <- c(B, rep(0.1, n))
  B <- c(B,1)
  
  f <- apply(X = R  ,
             MARGIN = 2 ,
             FUN = mean)
  f <- c(f,rep(0,m*m))
  
  #Pocitanie simplexovou metodou
  st <- Sys.time()
  vysledok_IPM <- Rglpk_solve_LP(f, A, c(rep("<=",2*m^2+m+n+n), "=="), B, max = TRUE, method = "interior")
  
  #Prvy vystup
  optimalne_riesenie <- vysledok_IPM$solution
  optimalne_riesenie <- as.matrix(optimalne_riesenie)
  et <- Sys.time()
  
  #Druhy vystup
  cas <- et - st
  
  #Treti vystup
  riesenie <- NULL
  if (vysledok_IPM$status == 0){
    riesenie <- "Solution found."
  }
  else{
    riesenie <- "Error occured during solving."
  }
  
  #Stvrty vystup 
  optimalna_hodnota <- vysledok_IPM$optimum
  
  return(list(optimalne_riesenie =optimalne_riesenie, cas = cas, riesenie = riesenie, optimalne_hodnota = optimalna_hodnota))
}

##############################################################################################################################