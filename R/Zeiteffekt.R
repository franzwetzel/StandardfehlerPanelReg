#' Standardfehler bei einem Zeiteffekt
#'
#' Die Funktion berechnet für die angegebene Anzahl an Simulationen die wahren Standardabweichungen des OLS und Fama-MacBeth Schätzers und die geschätzten Standardabweichungen nach OLS, Fama-MacBeth, Cluster, Cluster plus Schock und Newey-West und bildet darüber den Durchschnitt.
#' @details Da die Anzahl der Regressoren erst in den Argumenten der Funktion bestimmt wird, wird der Nutzer nach dem starten der Funktion aufgefordert, die Korrelationen der Regressoren anzugeben.
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T T bestimmt die Anzahl der Jahre im Paneldatensatz
#' @param Anzahl_Regressoren Legt die Anzahl der Regressoren fest
#' @param Anz_Sim Legt die Anzahl der Simulationen fest. Der Parameter bestimmt wie oft die wahren Standardfehler und die geschätzten für einen Paneldatensatz mit N Unternehmen und T Zeitperioden berechnet werden sollen.
#' @param Anteil_delta Anteil des Zeiteffekts an der gesamten Varianz der Fehlervariable
#' @returns Matrix mit den wahren und geschätzten Standardfehlern
#' @importFrom plm plm
#' @export
Zeiteffekt <- function(N = 500, T = 10, Anzahl_Regressoren = 1, Anz_Sim = 100, Anteil_delta = 0.25){
  beta <- 1

  # Varianz und Korrelation der erklärenden Variablen und der Fehlervariable

  for(i in 1:Anzahl_Regressoren){
    assign(paste0("AnteilVarZeitEffektX",i),
           as.double(readline(prompt = paste("Geben Sie den Anteil der Varianz des Zeiteffekts an der gesamten Varianz von ", paste0("X",i) ,"an ")))
    )
  }
  for(r in 1:Anzahl_Regressoren){
    assign(paste0("Var_X", r), 1^2)
    assign("x", get(paste0("AnteilVarZeitEffektX",r)))
    assign("y", get(paste0("Var_X",r)))
    assign(paste0("Var_zeta", r), x*y)
  }
  Var_eta <- 2^2
  Var_delta <- Anteil_delta*Var_eta


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
  Cov_Matrix <- matrix(rep(0, times = N*T*N*T), ncol = N*T)
  for(i in 1:N){
    for(t in 1:T){
      for(j in 1:N){
        for(s in 1:T){
          if(T*(i-1)+t==T*(j-1)+s){
            Cov_Matrix[T*(i-1)+t,T*(j-1)+s] <- Var_eta
          }
          else if(i != j & t==s){
            Cov_Matrix[T*(i-1)+t,T*(j-1)+s] <- Var_delta
          }
          else{
            Cov_Matrix[T*(i-1)+t,T*(j-1)+s] <- 0
          }
        }
      }
    }
  }

  # Fama
  Cov_Mat_Fama <- matrix(rep(0, times = N*N), ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      if(i==j){
        Cov_Mat_Fama[i,j] <- Var_eta
      }
      else{
        Cov_Mat_Fama[i,j] <- Var_delta
      }
    }
  }

  # Vorbereitung Simulation

  Ergeb_Mat <- vector(mode = "list", length = Anzahl_Regressoren)
  for(i in 1:Anzahl_Regressoren){
    x <- matrix(rep(0, times = Anz_Sim*6), ncol = 6)
    colnames(x) <- c("OLS", "Fama-MacBeth", "Cluster", "Newey-West",
                     "wahre Standardabweichung OLS", "wahre Standardabweichung Fama")
    names(Ergeb_Mat)[[i]] <- paste0("Ergeb_Mat_b", i)
    Ergeb_Mat[[i]] <- x
  }


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    for(r in 1:Anzahl_Regressoren){
      assign(paste0("Varzeta"), get(paste0("Var_zeta",r)))
      zeta <- rnorm(T, mean = 0, sd = sqrt(Varzeta))
      zeta_neu <- vector(mode = "double", length = N*T)
      for(i in 1:N){
        for(t in 1:T){
          zeta_neu[(i-1)*T+t] <- zeta[t]
        }
      }
      assign(paste0("VarX"), get(paste0("Var_X",r)))
      nu <- rnorm(N*T, mean = 0, sd = sqrt(VarX-Varzeta))

      assign(paste0("X", r), zeta_neu+nu)
    }

    delta <- rnorm(T, mean = 0, sd = sqrt(Var_delta))
    delta_neu <- vector(mode = "double", length = N*T)
    for(i in 1:N){
      for(t in 1:T){
        delta_neu[(i-1)*T+t] <- delta[t]
      }
    }
    epsilon <- rnorm(N*T, mean = 0, sd = sqrt(Var_eta-Var_delta))
    eta <- epsilon+delta_neu

    # Regressorenmatrix ertsellen

    X_Mat <- matrix(X1)
    if(Anzahl_Regressoren>1){
      for(i in 2:Anzahl_Regressoren){
        X_Mat <- cbind(X_Mat, get(paste0("X",i)))
      }
    }


    # y Variable erstellen

    y <- eta
    for(i in 1:Anzahl_Regressoren){
      assign(paste0("y"), y+beta*get(paste0("X",i)))
    }

    # Paneldaten Regression durchführen

    multReg_data <- data.frame(X_Mat, unternehmen, jahr, y)
    multReg_data_names <- vector(mode = "character", length = Anzahl_Regressoren + 3)
    multReg_data_names[(Anzahl_Regressoren+1):(Anzahl_Regressoren+3)] <- c("unternehmen", "jahr", "y")
    for(i in 1:Anzahl_Regressoren){
      multReg_data_names[i] <- paste0("X",i)
    }
    colnames(multReg_data) <- multReg_data_names
    multReg_paneldata <- plm::pdata.frame(multReg_data, index = c("unternehmen","jahr"))

    RegGleich <- paste0("y ~ 0 + X1")
    if(Anzahl_Regressoren > 1){
      for(i in 2:Anzahl_Regressoren){
        assign(paste0("RegGleich"), paste0(RegGleich, paste0(" + X",i)))
      }
    }
    multPanelReg <- plm::plm(as.formula(RegGleich), data = multReg_paneldata, model = "pooling")

    mult_Reg_Residuen <- matrix(resid(multPanelReg))

    # Standardfehler berechnen

    # OLS Standardfehler

    multReg_sd_OLS <- lmtest::coeftest(multPanelReg)[,2]
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,1] <- multReg_sd_OLS[i]
    }

    # Fama-MacBeth Standardfehler
    Fama_multReg_paneldata <- plm::pdata.frame(multReg_data, index = c("jahr", "unternehmen"))

    multFamaReg <- plm::pmg(as.formula(RegGleich), data = Fama_multReg_paneldata, model = "mg")

    multReg_sd_Fama <- lmtest::coeftest(multFamaReg)[,2]
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,2] <- multReg_sd_Fama[i]
    }

    # Cluster Standardfehler
    multPanelReg2 <- plm::plm(as.formula(RegGleich), data = Fama_multReg_paneldata, model = "pooling")

    multReg_VarMat_Cluster <- plm::vcovHC(multPanelReg2)
    multReg_Var_Cluster <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      multReg_Var_Cluster[i] <- multReg_VarMat_Cluster[i,i]
    }
    multReg_sd_Cluster <- sqrt(multReg_Var_Cluster)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,3] <- multReg_sd_Cluster[i]
    }

    # Newey-West Standardfehler

    multReg_VarMat_Newey <- plm::vcovNW(multPanelReg, maxlag = (T-1))
    multReg_Var_Newey <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      multReg_Var_Newey[i] <- multReg_VarMat_Newey[i,i]
    }
    multReg_sd_Newey <- sqrt(multReg_Var_Newey)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,4] <- multReg_sd_Newey[i]
    }


    # wahre Varianz
    # Var(b|X)=(X'X)^{-1}X'\Omega X(X'X)^{-1}

    Varianzmatrix <- solve(t(X_Mat)%*%X_Mat)%*%t(X_Mat)%*%Cov_Matrix%*%X_Mat%*%solve(t(X_Mat)%*%X_Mat)
    multReg_wahreVar <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      multReg_wahreVar[i] <- Varianzmatrix[i,i]
    }
    multReg_wahresd <- sqrt(multReg_wahreVar)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,5] <- multReg_wahresd[i]
    }


    # wahre Varianz Fama

    X_Mat_t <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_t[[t]] <- data.matrix(as.data.frame(Fama_multReg_paneldata[(N*(t-1)+1):(N*t),1:Anzahl_Regressoren]))
    }
    X_Mat_sum <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_sum[[t]] <- solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])%*%
        t(X_Mat_t[[t]])%*%Cov_Mat_Fama%*%X_Mat_t[[t]]%*%solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])
    }

    Varianzmatrix_Fama <- (1/(T^2))*(Reduce("+", X_Mat_sum))

    multReg_wahreVar_Fama <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      multReg_wahreVar_Fama[i] <- Varianzmatrix_Fama[i,i]
    }
    multReg_wahresd_Fama <- sqrt(multReg_wahreVar_Fama)

    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,6] <- multReg_wahresd_Fama[i]
    }
  }


  for(i in 1:Anzahl_Regressoren){
    x <- Ergeb_Mat[[i]]
    Ergeb_Mat[[i]] <- colMeans(x)
  }

  return(Ergeb_Mat)

}