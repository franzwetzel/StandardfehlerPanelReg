#' Standardfehler bei einem Zeiteffekt und Schock
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
Zeiteffekt_Schock <- function(N = 500, T = 10, Anz_Sim = 100,
                              Anteil_zeta = 0.25, Anteil_delta = 0.25,
                              lag_X = 1, lag_eta = 1){

  beta <- 1

  Var_X <- 1
  Var_eta <- 2^2

  Cov_Cluster_X <- Anteil_zeta*Var_X
  Cov_Cluster_eta <- Anteil_delta*Var_eta

  Var_zeta <- Cov_Cluster_X/lag_X
  Var_delta <- Cov_Cluster_eta/lag_eta

  # Fehlermeldungen

  if((T %% lag_X) != 0 | (T %% lag_eta) != 0){
    stop(paste(N, "ist nicht durch", lag_eta, "teilbar"))
  }

  if(Var_X < (lag_eta*Var_zeta) | Var_eta < (lag_eta*Var_delta)){
    stop("Die Varianz der Unternehmenseffekte ist größer als die Varianz des Regressors oder der Fehlervariable")
  }


  # Anteil von individuellem und gemeinsamem Effekt an der gesamten Varianz der Regressoren

  AnteilVarJointEffekt <- 1/1^2


  # Unternehmens- und Zeitvariable erstellen

  unternehmen1 <- vector(mode = "double", length = N*T)
  for(i in 1:N){
    for(t in 1:T){
      unternehmen1[T*(i-1)+t] <- i
    }
  }
  jahr1 <- vector(mode = "double", length = N*T)
  for(i in 1:N){
    for(t in 1:T){
      jahr1[T*(i-1)+t] <- 2000+t
    }
  }

  unternehmen2 <- vector(mode = "double", length = N*T)
  for(t in 1:T){
    for(i in 1:N){
      unternehmen2[N*(t-1)+i] <- i
    }
  }
  jahr2 <- vector(mode = "double", length = N*T)
  for(t in 1:T){
    for(i in 1:N){
      jahr2[N*(t-1)+i] <- 2000+t
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
          if(T*(i-1)+t==T*(j-1)+s){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- Var_eta
          }
          else if(t==s){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- lag_eta*Var_delta
          }
          else if(lag_eta > 1 & x > 0 & x <= (lag_eta-1)){
            Cov_Matrix_OLS[T*(i-1)+t,T*(j-1)+s] <- (lag_eta-x)*Var_delta
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
        Cov_Matrix_Fama1[i,j] <- lag_eta*Var_delta
      }
    }
  }
  Cov_Matrix_Fama2 <- vector(mode = "list", length = lag_eta)
  Cov_Matrix_Fama2[[1]] <- matrix(rep(0, times = N*N), ncol = N)
  if(lag_eta>1){
    for(i in 2:lag_eta){
      Cov_Matrix_Fama2[[i]] <- matrix(rep((lag_eta-(i-1))*Var_delta, times = N*N), ncol = N)
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
      assign(paste0("zeta", r), rnorm((T/lag_X)+1, mean = 0, sd = sqrt(Var_zeta)))
      assign(paste0("zeta", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("zeta", r, "_neu")))
      assign("y", get(paste0("zeta", r)))
      for(t in 1:((T/lag_X)-1)){
        for(i in 1:(N*lag_X)){
          x[(r-1)*N+(t-1)*(N*lag_X)+i] <- y[t]
        }
      }
      if(r == 1){
        for(z in 1:(N*lag_X)){
          x[(r-1)*N+((T/lag_X)-1)*(N*lag_X)+z] <- y[(T/lag_X)]
        }
      }else{
        for(l in 1:(N*(r-1))){
          x[l] <- y[(T/lag_X)]
        }
        for(f in 1:(N*(lag_X-(r-1)))){
          x[(r-1)*N+((T/lag_X)-1)*(N*lag_X)+f] <- y[(T/lag_X)+1]
        }
      }

      assign(paste0("zeta", r, "_neu"), x)
    }


    nu <- rnorm(N*T, mean = 0, sd = sqrt(Var_X-lag_X*Var_zeta))
    X <- nu
    for(r in 1:lag_X){
      assign("X", X + get(paste0("zeta", r, "_neu")))
    }



    # eta

    for(r in 1:lag_eta){
      assign(paste0("delta", r), rnorm((T/lag_eta)+1, mean = 0, sd = sqrt(Var_delta)))
      assign(paste0("delta", r, "_neu"), vector(mode = "double", length = N*T))

      assign("x", get(paste0("delta", r, "_neu")))
      assign("y", get(paste0("delta", r)))
      for(t in 1:((T/lag_eta)-1)){
        for(i in 1:(N*lag_eta)){
          x[(r-1)*N+(t-1)*(N*lag_eta)+i] <- y[t]
        }
      }
      if(r == 1){
        for(z in 1:(N*lag_eta)){
          x[(r-1)*N+((T/lag_eta)-1)*(N*lag_eta)+z] <- y[(T/lag_eta)]
        }
      }else{
        for(l in 1:(N*(r-1))){
          x[l] <- y[(T/lag_eta)]
        }
        for(f in 1:(N*(lag_eta-(r-1)))){
          x[(r-1)*N+((T/lag_eta)-1)*(N*lag_eta)+f] <- y[(T/lag_eta)+1]
        }
      }

      assign(paste0("delta", r, "_neu"), x)
    }


    epsilon <- rnorm(N*T, mean = 0, sd = sqrt(Var_eta-lag_eta*Var_delta))
    eta <- epsilon
    for(r in 1:lag_eta){
      assign("eta", eta + get(paste0("delta", r, "_neu")))
    }



    # y Variable erstellen

    y <- 0 + X + eta

    # Paneldaten Regression durchführen

    daten <- data.frame(X, y, unternehmen2, jahr2)
    colnames(daten) <- c("X", "y", "unternehmen", "jahr")

    paneldaten <- plm::pdata.frame(daten, index = c("unternehmen", "jahr"))

    paneldaten_Fama <- plm::pdata.frame(daten, index = c("jahr", "unternehmen"))

    Reg <- plm::plm(y ~ 0 + X, data = paneldaten, model = "pooling")

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

    Var_Cluster <- plm::vcovHC(Reg, cluster = "time")
    Sd_Cluster <- sqrt(Var_Cluster)
    Ergeb_Matrix[D,3] <- Sd_Cluster


    # Cluster und Schock

    Var_ClusterSchock <- plm::vcovSCC(Reg, cluster = "time",
                                      inner = "cluster", maxlag = (lag_eta-1))
    Sd_ClusterSchock <- sqrt(Var_ClusterSchock)
    Ergeb_Matrix[D,4] <- Sd_ClusterSchock


    # Newey-West Standardfehler

    Var_Newey <- plm::vcovNW(Reg, maxlag = (T-1), cluster = "time")
    Sd_Newey <- sqrt(Var_Newey)
    Ergeb_Matrix[D,5] <- Sd_Newey


    # wahre Varianz OLS
    # Var(b|X)=(X'X)^{-1}X'\Omega X(X'X)^{-1}
    X <- as.matrix(as.vector(paneldaten[,1]))

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
        x <- abs(t-s)
        if(x <= (lag_eta-1)){
          X_sum2[[index]] <- solve(t(X_t[[t]])%*%X_t[[t]])%*%t(X_t[[t]])%*%
            Cov_Matrix_Fama2[[x+1]]%*%X_t[[s]]%*%solve(t(X_t[[s]])%*%X_t[[s]])
        }
        else{
          X_sum2[[index]] <- solve(t(X_t[[t]])%*%X_t[[t]])%*%t(X_t[[t]])%*%
            Cov_Matrix_Fama2[[1]]%*%X_t[[s]]%*%solve(t(X_t[[s]])%*%X_t[[s]])
        }
      }
    }

    wahre_Var_Fama <- (1/(T^2))*(Reduce("+", X_sum1))+2*(1/(T^2))*(Reduce("+", X_sum2))
    wahre_Sd_Fama <- sqrt(wahre_Var_Fama)
    Ergeb_Matrix[D,7] <- wahre_Sd_Fama

  }

  Ergeb_Mat <- colMeans(Ergeb_Matrix)

  return(Ergeb_Mat)

}
