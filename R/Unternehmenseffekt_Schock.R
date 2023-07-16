#' Standardfehler bei Korrelation über die Unternehmenscluster hinaus
#'
#' Die Funktion berechnet die geschätzten und wahren Standardfehler bei erklärenden Variablen und Fehlervariablen,
#' die über die Unternehmenscluster hinaus korreliert sind.
#' @details Da die Anzahl der Regressoren nicht die Fehlschätzung beeinflusst,
#' simuliert diese Funktion nur Regressionen mit einem Regressor
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T T bestimmt die Anzahl der Perioden im Paneldatensatz
#' @param Anz_Sim Anzahl der Simulationen. Anz_Sim bestimmt, wieviele Paneldatensätze
#' mit N Unternehmen und T Perioden simuliert werden sollen.
#' @param Korr_mu Höhe der Korrelation innerhalb der Unternehmenscluster der erklärenden Variable
#' (Anteil der Varianz des Unternehmenseffekts an der gesamten Varianz der erklärenden Variable)
#' @param Korr_u Höhe der Korrelation innerhalb der Unternehmenscluster der Fehlervariable
#' (Anteil der Varianz des Unternehmenseffekts an der gesamten Varianz der Fehlervariable)
#' @param D_X Distanz |i-j| für die, die erklärenden Variablen verschiedener Unternehmen korreliert sind
#' @param D_eta Distanz |i-j| für die, die Fehlervariablen verschiedener Unternehmen korreliert sind
#' @param L Höhe des L des Discroll-Kraay Schätzers
#' @returns Vektor mit Länge 8. Durchschnittliche Schätzung der Standardfehler nach OLS, Fama-MacBeth, Cluster,
#' Discroll-Kraay, Cluster-Schock und Newey-West sowie
#' der durchschnittlichen wahren Standardfehler der OLS Regression und der Fama-MacBeth Regression
#' @importFrom plm plm
#' @export
Unternehmenseffekt_Schock <- function(N = 500, T = 10, Anz_Sim = 100,
                                      Korr_mu = 0.25, Korr_u = 0.25,
                                      D_X = 0, D_eta = 0, L = 10){

  D_X <- D_X + 1
  D_eta <- D_eta + 1

  beta <- 1

  Var_X <- 1
  Var_eta <- 2^2

  Cov_Cluster_X <- Korr_mu*Var_X
  Cov_Cluster_eta <- Korr_u*Var_eta

  Var_mu <- Cov_Cluster_X/D_X
  Var_u <- Cov_Cluster_eta/D_eta

  # Fehlermeldungen

  if((N %% D_X) != 0 | (N %% D_eta) != 0){
    stop(paste(N, "ist nicht durch", D_eta, "teilbar"))
  }

  if(Var_X < (D_eta*Var_mu) | Var_eta < (D_eta*Var_u)){
    stop("Die Varianz der Unternehmenseffekte ist größer als die Varianz des Regressors oder der Fehlervariable")
  }

  # Unternehmens- und Zeitvariable erstellen

  unternehmen <- vector(mode = "double", length = N*T)
  for(i in 1:N){
    for(t in 1:T){
      unternehmen[T*(i-1)+t] <- i
    }
  }
  jahr <- vector(mode = "double", length = N*T)
  for(i in 1:N){
    for(t in 1:T){
      jahr[T*(i-1)+t] <- 2000+t
    }
  }


  # Kovarianzmarix erstellen

  # OLS
  Cov_Matrix_OLS <- matrix(rep(0, times = N*T*N*T), ncol = N*T)
  for(i in 1:N){
    for(t in 1:T){
      for(j in 1:N){
        for(s in 1:T){
          x <- abs(i-j)
          if(T*(i-1)+t==T*(j-1)+s){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- Var_eta
          }
          else if(i==j){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- D_eta*Var_u
          }
          else if(D_eta > 1 & x > 0 & x <= (D_eta-1)){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- (D_eta-x)*Var_u
          }
          else{
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- 0
          }
        }
      }
    }
  }

  # Fama
  Cov_Matrix_Fama1 <- matrix(rep(0, times = N*N), ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      x <- abs(i-j)
      if(i==j){
        Cov_Matrix_Fama1[i,j] <- Var_eta
      }
      else if(D_eta > 1 & x > 0 & x <= (D_eta-1)){
        Cov_Matrix_Fama1[i,j] <- (D_eta-x)*Var_u
      }
      else{
        Cov_Matrix_Fama1[i,j] <- 0
      }
    }
  }
  Cov_Matrix_Fama2 <- matrix(rep(0, times = N*N), ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      x <- abs(i-j)
      if(i==j){
        Cov_Matrix_Fama2[i,j] <- D_eta*Var_u
      }
      else if(D_eta > 1 & x > 0 & x <= (D_eta-1)){
        Cov_Matrix_Fama1[i,j] <- (D_eta-x)*Var_u
      }
      else{
        Cov_Matrix_Fama2[i,j] <- 0
      }
    }
  }


  # Vorbereitung Simulation

  Ergeb_Matrix <- matrix(rep(0, times = Anz_Sim*8), ncol = 8)
  colnames(Ergeb_Matrix) <- c("OLS", "Fama-MacBeth", "Cluster", "Driscoll-Kraay", "Cluster, Schock",
                              "Newey-West", "wahre Standardabweichung OLS", "wahre Standardabweichung Fama")


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    # Simulation

    # X

    for(r in 1:D_X){
      assign(paste0("mu", r), stats::rnorm((N/D_X)+1, mean = 0, sd = sqrt(Var_mu)))
      assign("y", get(paste0("mu", r)))
      assign(paste0("mu", r, "_neu"), vector(mode = "double", length = N*T))
      assign("x", get(paste0("mu", r, "_neu")))

      for(i in 1:((N/D_X)-1)){
        for(t in 1:(T*D_X)){
          x[(r-1)*T+(i-1)*(T*D_X)+t] <- y[i]
        }
      }
      if(r == 1){
        for(z in 1:(T*D_X)){
          x[(r-1)*T+((N/D_X)-1)*(T*D_X)+z] <- y[(N/D_X)]
        }
      }else{
        for(l in 1:(T*(r-1))){
          x[l] <- y[(N/D_X)]
        }
        for(f in 1:(T*(D_X-(r-1)))){
          x[(r-1)*T+((N/D_X)-1)*(T*D_X)+f] <- y[(N/D_X)+1]
        }
      }

      assign(paste0("mu", r, "_neu"), x)
    }


    nu <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_X-D_X*Var_mu))
    X <- nu
    for(r in 1:D_X){
      assign("X", X + get(paste0("mu", r, "_neu")))
    }



    # eta

    for(r in 1:D_eta){
      assign(paste0("u", r), stats::rnorm((N/D_eta)+1, mean = 0, sd = sqrt(Var_u)))
      assign("y", get(paste0("u", r)))
      assign(paste0("u", r, "_neu"), vector(mode = "double", length = N*T))
      assign("x", get(paste0("u", r, "_neu")))

      for(i in 1:((N/D_eta)-1)){
        for(t in 1:(T*D_eta)){
          x[(r-1)*T+(i-1)*(T*D_eta)+t] <- y[i]
        }
      }
      if(r == 1){
        for(z in 1:(T*D_eta)){
          x[(r-1)*T+((N/D_eta)-1)*(T*D_eta)+z] <- y[(N/D_eta)]
        }
      }else{
        for(l in 1:(T*(r-1))){
          x[l] <- y[(N/D_eta)]
        }
        for(f in 1:(T*(D_eta-(r-1)))){
          x[(r-1)*T+((N/D_eta)-1)*(T*D_eta)+f] <- y[(N/D_eta)+1]
        }
      }

      assign(paste0("u", r, "_neu"), x)
    }


    epsilon <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_eta-D_eta*Var_u))
    eta <- epsilon
    for(r in 1:D_eta){
      assign("eta", eta + get(paste0("u", r, "_neu")))
    }




    # y Variable erstellen

    y <- X + eta

    # Paneldaten Regression durchführen

    daten <- data.frame(X, y, unternehmen, jahr)
    colnames(daten) <- c("X", "y", "unternehmen", "jahr")

    paneldaten <- plm::pdata.frame(daten, index = c("unternehmen", "jahr"))
    Reg <- plm::plm(y ~ X, data = paneldaten, model = "pooling")

    paneldaten_Fama <- plm::pdata.frame(daten, index = c("jahr", "unternehmen"))
    Reg2 <- plm::plm(y ~ X, data = paneldaten_Fama, model = "pooling")

    # Standardfehler berechnen

    # OLS Standardfehler

    Sd_OLS <- lmtest::coeftest(Reg)[2,2]
    Ergeb_Matrix[D,1] <- Sd_OLS


    # Fama-MacBeth Standardfehler

    Reg_Fama <- plm::pmg(y ~ X, data = paneldaten_Fama, model = "mg")

    Sd_Fama <- lmtest::coeftest(Reg_Fama)[2,2]
    Ergeb_Matrix[D,2] <- Sd_Fama


    # Cluster Standardfehler

    Var_Cluster <- plm::vcovHC(Reg)
    Sd_Cluster <- sqrt(Var_Cluster[2,2])
    Ergeb_Matrix[D,3] <- Sd_Cluster


    # Driscoll-Kraay

    #maxlag <- L
    #wj <-  function(j, maxlag) 1 - j/(maxlag + 1)
    #Var_Schock <- matrix(rep(0, times = 4), ncol = 2)
    #Var_ClusterSchock <- plm::vcovHC(Reg, cluster = "group", type = "sss")
    #if(maxlag > 0){
    #  for(i in 1:maxlag){
    #    Vctl <- plm::vcovG(Reg, type = "sss", cluster = "group",
    #                       l = i, inner = "cluster")
    #    Var_Schock <- Var_Schock + wj(i, maxlag) * (Vctl + t(Vctl))
    #  }
    #}
    #Var_ClusterSchock <- Var_ClusterSchock + Var_Schock

    Var_D_K <- plm::vcovSCC(Reg, cluster = "group", type = "sss",
                            inner = "cluster", maxlag = L)
    Sd_D_K <- sqrt(Var_D_K[2,2])
    Ergeb_Matrix[D,4] <- Sd_D_K


    # Cluster und Schock

    w1 <- function(j, maxlag) 1
    Var_Cluster_S <- plm::vcovSCC(Reg, cluster = "group", wj = w1, type = "sss",
                                  inner = "cluster", maxlag = (D_eta-1))
    Sd_Cluster_S <- sqrt(Var_Cluster_S[2,2])
    Ergeb_Matrix[D,5] <- Sd_Cluster_S


    # Newey-West Standardfehler

    Var_Newey <- plm::vcovNW(Reg, maxlag = (T-1))
    Sd_Newey <- sqrt(Var_Newey[2,2])
    Ergeb_Matrix[D,6] <- Sd_Newey


    # wahre Varianz OLS
    # Var(b|X)=(X'X)^{-1}X'\Omega X(X'X)^{-1}

    intercept <- as.matrix(rep(1, times = N*T))
    X <- cbind(intercept, as.matrix(as.vector(paneldaten[,1])))

    wahre_Varianz_OLS <- solve(t(X)%*%X)%*%t(X)%*%Cov_Matrix_OLS%*%X%*%solve(t(X)%*%X)

    wahre_Sd_OLS <- sqrt(wahre_Varianz_OLS[2,2])
    Ergeb_Matrix[D,7] <- wahre_Sd_OLS


    # wahre Varianz Fama

    X_Mat_Fama <- cbind(intercept, as.matrix(as.data.frame(paneldaten_Fama[,1])))

    X_t <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_t[[t]] <- data.matrix(as.data.frame(X_Mat_Fama[(N*(t-1)+1):(N*t),]))
    }
    X_sum1 <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_sum1[[t]] <- solve(t(X_t[[t]])%*%X_t[[t]])%*%
        t(X_t[[t]])%*%Cov_Matrix_Fama1%*%X_t[[t]]%*%solve(t(X_t[[t]])%*%X_t[[t]])
    }
    X_sum2 <- vector(mode = "list", length = ((T-1)*T)/2)
    index <- 0
    for(t in 1:(T-1)){
      for(s in (t+1):T){
        index <- index + 1
        X_sum2[[index]] <- solve(t(X_t[[t]])%*%X_t[[t]])%*%t(X_t[[t]])%*%
          Cov_Matrix_Fama2%*%X_t[[s]]%*%solve(t(X_t[[s]])%*%X_t[[s]])
      }
    }

    wahre_Var_Fama <- (1/(T^2))*(Reduce("+", X_sum1))+2*(1/(T^2))*(Reduce("+", X_sum2))
    wahre_Sd_Fama <- sqrt(wahre_Var_Fama[2,2])
    Ergeb_Matrix[D,8] <- wahre_Sd_Fama

  }

  Ergeb_Mat <- colMeans(Ergeb_Matrix, na.rm = TRUE)

  return(Ergeb_Mat)

}
