#' Standardfehler bei einem Unternehmenseffekt und einem Schock
#'
#' Die Funktion berechnet für die angegebene Anzahl an Simulationen die wahren Standardabweichungen des OLS und Fama-MacBeth Schätzers und die geschätzten Standardabweichungen nach OLS, Fama-MacBeth, Cluster, Cluster plus Schock und Newey-West und bildet darüber den Durchschnitt.
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T N bestimmt die Anzahl der Jahre im Paneldatensatz
#' @param Anzahl_Regressoren Legt die Anzahl der Regressoren fest
#' @param Anz_Sim Legt die Anzahl der Simulationen fest. Der Parameter bestimmt wie oft die wahren Standardfehler und die geschätzten für einen Paneldatensatz mit N Unternehmen und T Zeitperioden berechnet werden sollen.
#' @returns Für jede Schätzung eines Betas die wahren Standardfehler für die OLS Regression und die Fama-MacBeth Regression
#' sowie die Schätzung der Standardfehler nach
#' @importFrom plm plm
#' @export
Unternehmenseffekt_Schock <- function(N = 500, T = 10, Anz_Sim = 100,
                                      Anteil_mu = 0.25, Anteil_u = 0.25,
                                      lag_X = 1, lag_eta = 1){

  beta <- 1

  Var_X <- 1
  Var_eta <- 2^2

  Cov_Cluster_X <- Anteil_mu*Var_X
  Cov_Cluster_eta <- Anteil_u*Var_eta

  Var_mu <- Cov_Cluster_X/lag_X
  Var_u <- Cov_Cluster_eta/lag_eta

  # Fehlermeldungen

  if((N %% lag_X) != 0 | (N %% lag_eta) != 0){
    stop(paste(N, "ist nicht durch", lag_eta, "teilbar"))
  }

  if(Var_X < (lag_eta*Var_mu) | Var_eta < (lag_eta*Var_u)){
    stop("Die Varianz der Unternehmenseffekte ist größer als die Varianz des Regressors oder der Fehlervariable")
  }


  # Anteil von individuellem und gemeinsamem Effekt an der gesamten Varianz der Regressoren

  AnteilVarJointEffekt <- 1/1^2


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
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- lag_eta*Var_u
          }
          else if(lag_eta > 1 & x > 0 & x <= (lag_eta-1)){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- (lag_eta-x)*Var_u
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
      else if(lag_eta > 1 & x > 0 & x <= (lag_eta-1)){
        Cov_Matrix_Fama1[i,j] <- (lag_eta-x)*Var_u
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
        Cov_Matrix_Fama2[i,j] <- lag_eta*Var_u
      }
      else if(lag_eta > 1 & x > 0 & x <= (lag_eta-1)){
        Cov_Matrix_Fama1[i,j] <- (lag_eta-x)*Var_u
      }
      else{
        Cov_Matrix_Fama2[i,j] <- 0
      }
    }
  }


  # Vorbereitung Simulation

  Ergeb_Matrix <- matrix(rep(0, times = Anz_Sim*7), ncol = 7)
  colnames(Ergeb_Matrix) <- c("OLS", "Fama-MacBeth", "Cluster", "Cluster, Schock", "Newey-West",
                              "wahre Standardabweichung OLS", "wahre Standardabweichung Fama")


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    # Simulation

    # X

    for(r in 1:lag_X){
      assign(paste0("mu", r), rnorm((N/lag_X)+1, mean = 0, sd = sqrt(Var_mu)))
      assign(paste0("mu", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("mu", r, "_neu")))
      assign("y", get(paste0("mu", r)))
      for(i in 1:((N/lag_X)-1)){
        for(t in 1:(T*lag_X)){
          x[(r-1)*T+(i-1)*(T*lag_X)+t] <- y[i]
        }
      }
      if(r == 1){
        for(z in 1:(T*lag_X)){
          x[(r-1)*T+((N/lag_X)-1)*(T*lag_X)+z] <- y[(N/lag_X)]
        }
      }else{
        for(l in 1:(T*(r-1))){
          x[l] <- y[(N/lag_X)]
        }
        for(f in 1:(T*(lag_X-(r-1)))){
          x[(r-1)*T+((N/lag_X)-1)*(T*lag_X)+f] <- y[(N/lag_X)+1]
        }
      }

      assign(paste0("mu", r, "_neu"), x)
    }


    nu <- rnorm(N*T, mean = 0, sd = sqrt(Var_X-lag_X*Var_mu))
    X <- nu
    for(r in 1:lag_X){
      assign("X", X + get(paste0("mu", r, "_neu")))
    }



    # eta

    for(r in 1:lag_eta){
      assign(paste0("u", r), rnorm((N/lag_eta)+1, mean = 0, sd = sqrt(Var_u)))
      assign(paste0("u", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("u", r, "_neu")))
      assign("y", get(paste0("u", r)))
      for(i in 1:((N/lag_eta)-1)){
        for(t in 1:(T*lag_eta)){
          x[(r-1)*T+(i-1)*(T*lag_eta)+t] <- y[i]
        }
      }
      if(r == 1){
        for(z in 1:(T*lag_eta)){
          x[(r-1)*T+((N/lag_eta)-1)*(T*lag_eta)+z] <- y[(N/lag_eta)]
        }
      }else{
        for(l in 1:(T*(r-1))){
          x[l] <- y[(N/lag_eta)]
        }
        for(f in 1:(T*(lag_eta-(r-1)))){
          x[(r-1)*T+((N/lag_eta)-1)*(T*lag_eta)+f] <- y[(N/lag_eta)+1]
        }
      }

      assign(paste0("u", r, "_neu"), x)
    }


    epsilon <- rnorm(N*T, mean = 0, sd = sqrt(Var_eta-lag_eta*Var_u))
    eta <- epsilon
    for(r in 1:lag_eta){
      assign("eta", eta + get(paste0("u", r, "_neu")))
    }




    # y Variable erstellen

    y <- 0 + X + eta

    # Paneldaten Regression durchführen

    daten <- data.frame(X, y, unternehmen, jahr)
    colnames(daten) <- c("X", "y", "unternehmen", "jahr")

    paneldaten <- plm::pdata.frame(daten, index = c("unternehmen", "jahr"))
    Reg <- plm::plm(y ~ 0 + X, data = paneldaten, model = "pooling")

    paneldaten_Fama <- plm::pdata.frame(daten, index = c("jahr", "unternehmen"))
    Reg2 <- plm::plm(y ~ 0 + X, data = paneldaten_Fama, model = "pooling")

    # Standardfehler berechnen

    # OLS Standardfehler

    Sd_OLS <- lmtest::coeftest(Reg)[1,2]
    Ergeb_Matrix[D,1] <- Sd_OLS


    # Fama-MacBeth Standardfehler

    Reg_Fama <- plm::pmg(y ~ 0 + X, data = paneldaten_Fama, model = "mg")

    Sd_Fama <- lmtest::coeftest(Reg_Fama)[1,2]
    Ergeb_Matrix[D,2] <- Sd_Fama


    # Cluster Standardfehler

    Var_Cluster <- plm::vcovHC(Reg)
    Sd_Cluster <- sqrt(Var_Cluster)
    Ergeb_Matrix[D,3] <- Sd_Cluster


    # Cluster und Schock

    Var_ClusterSchock <- plm::vcovSCC(Reg, cluster = "group",
                                      inner = "cluster", maxlag = (lag_eta-1))
    Sd_ClusterSchock <- sqrt(Var_ClusterSchock)
    Ergeb_Matrix[D,4] <- Sd_ClusterSchock


    # Newey-West Standardfehler

    Var_Newey <- plm::vcovNW(Reg, maxlag = (T-1))
    Sd_Newey <- sqrt(Var_Newey)
    Ergeb_Matrix[D,5] <- Sd_Newey


    # wahre Varianz OLS
    # Var(b|X)=(X'X)^{-1}X'\Omega X(X'X)^{-1}
    X <- as.matrix(X)

    wahre_Varianz_OLS <- solve(t(X)%*%X)%*%t(X)%*%Cov_Matrix_OLS%*%X%*%solve(t(X)%*%X)

    wahre_Sd_OLS <- sqrt(wahre_Varianz_OLS)
    Ergeb_Matrix[D,6] <- wahre_Sd_OLS


    # wahre Varianz Fama

    X_t <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_t[[t]] <- data.matrix(as.data.frame(paneldaten_Fama[(N*(t-1)+1):(N*t),1]))
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
    wahre_Sd_Fama <- sqrt(wahre_Var_Fama)
    Ergeb_Matrix[D,7] <- wahre_Sd_Fama

  }

  Ergeb_Mat <- colMeans(Ergeb_Matrix)

  return(Ergeb_Mat)

}
