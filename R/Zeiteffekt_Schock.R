#' Standardfehler bei Korrelation über die Zeitcluster hinaus
#'
#' Die Funktion berechnet die geschätzten und wahren Standardfehler bei erklärenden Variablen und Fehlervariablen,
#' die über die Zeitcluster hinaus korreliert sind.
#' @details Da die Anzahl der Regressoren nicht die Fehlschätzung beeinflusst,
#' simuliert diese Funktion nur Regressionen mit einem Regressor
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T T bestimmt die Anzahl der Perioden im Paneldatensatz
#' @param Anz_Sim Anzahl der Simulationen. Anz_Sim bestimmt, wieviele Paneldatensätze
#' mit N Unternehmen und T Perioden simuliert werden sollen.
#' @param Korr_zeta Höhe der Korrelation innerhalb der Zeitcluster der erklärenden Variable
#' (Anteil der Varianz des Zeiteffekts an der gesamten Varianz der erklärenden Variable)
#' @param Korr_delta Höhe der Korrelation innerhalb der Zeitcluster der Fehlervariable
#' (Anteil der Varianz des Zeiteffekts an der gesamten Varianz der Fehlervariable)
#' @param D_X Distanz |t-s| für die, die erklärenden Variablen verschiedener Perioden korreliert sind
#' @param D_eta Distanz |t-s| für die, die Fehlervariablen verschiedener Perioden korreliert sind
#' @param L Höhe des L des Discroll-Kraay Schätzers
#' @returns Vektor mit Länge 8. Durchschnittliche Schätzung der Standardfehler nach OLS, Fama-MacBeth, Cluster,
#' Discroll-Kraay, Cluster-Schock und Newey-West sowie
#' der durchschnittlichen wahren Standardfehler der OLS Regression und der Fama-MacBeth Regression
#' @importFrom plm plm
#' @export
Zeiteffekt_Schock <- function(N = 500, T = 10, Anz_Sim = 100,
                              Korr_zeta = 0.25, Korr_delta = 0.25,
                              D_X = 0, D_eta = 0, L = 9){

  D_X <- D_X + 1
  D_eta <- D_eta + 1

  beta <- 1

  Var_X <- 1
  Var_eta <- 2^2

  Cov_Cluster_X <- Korr_zeta*Var_X
  Cov_Cluster_eta <- Korr_delta*Var_eta

  Var_zeta <- Cov_Cluster_X/D_X
  Var_delta <- Cov_Cluster_eta/D_eta

  # Fehlermeldungen

  if((T %% D_X) != 0 | (T %% D_eta) != 0){
    stop(paste(N, "ist nicht durch", D_eta, "teilbar"))
  }

  if(Var_X < (D_eta*Var_zeta) | Var_eta < (D_eta*Var_delta)){
    stop("Die Varianz der Unternehmenseffekte ist größer als die Varianz des Regressors oder der Fehlervariable")
  }


  # Unternehmens- und Zeitvariable erstellen

  unternehmen <- vector(mode = "double", length = N*T)
  for(t in 1:T){
    for(i in 1:N){
      unternehmen[N*(t-1)+i] <- i
    }
  }
  jahr <- vector(mode = "double", length = N*T)
  for(t in 1:T){
    for(i in 1:N){
      jahr[N*(t-1)+i] <- 2000+t
    }
  }


  # Kovarianzmarix erstellen

  # OLS
  Cov_Matrix_OLS <- matrix(rep(0, times = N*T*N*T), ncol = N*T)
  for(i in 1:N){
    for(t in 1:T){
      for(j in 1:N){
        for(s in 1:T){
          x <- abs(t-s)
          if(i==j & t==s){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- Var_eta
          }
          else if(x==0){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- D_eta*Var_delta
          }
          else if(D_eta > 1 & x > 0 & x <= (D_eta-1)){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- (D_eta-x)*Var_delta
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
      if(i==j){
        Cov_Matrix_Fama1[i,j] <- Var_eta
      }
      else{
        Cov_Matrix_Fama1[i,j] <- 0
      }
    }
  }
  Cov_Matrix_Fama2 <- vector(mode = "list", length = D_eta)
  Cov_Matrix_Fama2[[1]] <- matrix(rep(0, times = N*N), ncol = N)
  if(D_eta>1){
    for(i in 2:D_eta){
      Cov_Matrix_Fama2[[i]] <- matrix(rep((D_eta-(i-1))*Var_delta, times = N*N), ncol = N)
    }
  }


  # Vorbereitung Simulation

  Ergeb_Matrix <- matrix(rep(0, times = Anz_Sim*8), ncol = 8)
  colnames(Ergeb_Matrix) <- c("OLS", "Fama-MacBeth", "Cluster", "Driscoll-Kraay", "Cluster-Schock",
                              "Newey-West", "wahre Standardabweichung OLS", "wahre Standardabweichung Fama")


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    # Simulation

    # X

    for(r in 1:D_X){
      assign(paste0("zeta", r), stats::rnorm((T/D_X)+1, mean = 0, sd = sqrt(Var_zeta)))
      assign(paste0("zeta", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("zeta", r, "_neu")))
      assign("y", get(paste0("zeta", r)))
      for(t in 1:((T/D_X)-1)){
        for(i in 1:(N*D_X)){
          x[(r-1)*N+(t-1)*(N*D_X)+i] <- y[t]
        }
      }
      if(r == 1){
        for(z in 1:(N*D_X)){
          x[(r-1)*N+((T/D_X)-1)*(N*D_X)+z] <- y[(T/D_X)]
        }
      }else{
        for(l in 1:(N*(r-1))){
          x[l] <- y[(T/D_X)]
        }
        for(f in 1:(N*(D_X-(r-1)))){
          x[(r-1)*N+((T/D_X)-1)*(N*D_X)+f] <- y[(T/D_X)+1]
        }
      }

      assign(paste0("zeta", r, "_neu"), x)
    }


    nu <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_X-D_X*Var_zeta))
    X <- nu
    for(r in 1:D_X){
      assign("X", X + get(paste0("zeta", r, "_neu")))
    }



    # eta

    for(r in 1:D_eta){
      assign(paste0("delta", r), stats::rnorm((T/D_eta)+1, mean = 0, sd = sqrt(Var_delta)))
      assign(paste0("delta", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("delta", r, "_neu")))
      assign("y", get(paste0("delta", r)))
      for(t in 1:((T/D_eta)-1)){
        for(i in 1:(N*D_eta)){
          x[(r-1)*N+(t-1)*(N*D_eta)+i] <- y[t]
        }
      }
      if(r == 1){
        for(z in 1:(N*D_eta)){
          x[(r-1)*N+((T/D_eta)-1)*(N*D_eta)+z] <- y[(T/D_eta)]
        }
      }else{
        for(l in 1:(N*(r-1))){
          x[l] <- y[(T/D_eta)]
        }
        for(f in 1:(N*(D_eta-(r-1)))){
          x[(r-1)*N+((T/D_eta)-1)*(N*D_eta)+f] <- y[(T/D_eta)+1]
        }
      }

      assign(paste0("delta", r, "_neu"), x)
    }


    epsilon <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_eta-D_eta*Var_delta))
    eta <- epsilon
    for(r in 1:D_eta){
      assign("eta", eta + get(paste0("delta", r, "_neu")))
    }



    # y Variable erstellen

    y <- X + eta

    # Paneldaten Regression durchführen

    daten <- data.frame(X, y, unternehmen, jahr)
    colnames(daten) <- c("X", "y", "unternehmen", "jahr")

    paneldaten <- plm::pdata.frame(daten, index = c("unternehmen", "jahr"))

    paneldaten_Fama <- plm::pdata.frame(daten, index = c("jahr", "unternehmen"))

    Reg <- plm::plm(y ~ X, data = paneldaten, model = "pooling")

    # Standardfehler berechnen

    # OLS Standardfehler

    Sd_OLS <- lmtest::coeftest(Reg)[2,2]
    Ergeb_Matrix[D,1] <- Sd_OLS


    # Fama-MacBeth Standardfehler

    Reg_Fama <- plm::pmg(y ~ X, data = paneldaten_Fama, model = "mg")

    Sd_Fama <- lmtest::coeftest(Reg_Fama)[2,2]
    Ergeb_Matrix[D,2] <- Sd_Fama


    # Cluster Standardfehler

    Var_Cluster <- plm::vcovHC(Reg, cluster = "time", type = "sss")
    Sd_Cluster <- sqrt(Var_Cluster[2,2])
    Ergeb_Matrix[D,3] <- Sd_Cluster


    # Discroll-Kraay

    #maxlag <- (D_eta-1)
    #wj <-  function(j, maxlag) 1 - j/(maxlag + 1)
    #Var_Schock <- matrix(rep(0, times = 4), ncol = 2)
    #Var_D_K <- plm::vcovHC(Reg, cluster = "time", type = "sss")
    #if(maxlag > 0){
    #  for(i in seq_len(maxlag)){
    #    Vctl <- plm::vcovG(Reg, type = "sss", cluster = "time",
    #                       l = i, inner = "cluster")
    #    Var_Schock <- Var_Schock + wj(i, maxlag) * (Vctl + t(Vctl))
    #  }
    #}
    #Var_D_K <- Var_D_K + Var_Schock

    Var_D_K <- plm::vcovSCC(Reg, cluster = "time", maxlag = L, type = "sss")

    Sd_D_K <- sqrt(Var_D_K[2,2])
    Ergeb_Matrix[D,4] <- Sd_D_K


    # Cluster und Schock

    w1 <- function(j, maxlag) 1
    Var_Cluster_S <- plm::vcovSCC(Reg, cluster = "time", wj = w1, type = "sss",
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
      X_t[[t]] <- X_Mat_Fama[(N*(t-1)+1):(N*t),]
    }
    X_sum1 <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_sum1[[t]] <- (Var_eta-D_eta*Var_delta)*solve(t(X_t[[t]])%*%X_t[[t]])
    }

    wahre_Var_Fama <- (1/(T^2))*(Reduce("+", X_sum1))
    wahre_Sd_Fama <- sqrt(wahre_Var_Fama[2,2])
    Ergeb_Matrix[D,8] <- wahre_Sd_Fama

  }

  Ergeb_Mat <- colMeans(Ergeb_Matrix, na.rm = TRUE)

  return(Ergeb_Mat)

}

