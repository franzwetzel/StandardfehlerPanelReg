#' Standardfehler bei einem Zeiteffekt
#'
#' Die Funktion berechnet die geschätzten und wahren Standardfehler bei einem Zeiteffekt.
#' @details Die Anzahl der Regressoren ist ein Argument der Funktion.
#' Nach dem Starten der Funtion wird der Nutzer aufgefordert, für jeden Regressor die
#' Korrelation innerhalb der Zeitcluster anzugeben.
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T T bestimmt die Anzahl der Perioden im Paneldatensatz
#' @param Anzahl_Regressoren Anzahl der Regressoren
#' @param Anz_Sim Anzahl der Simulationen. Anz_Sim bestimmt, wieviele Paneldatensätze
#' mit N Unternehmen und T Perioden simuliert werden sollen.
#' @param Korr_delta Höhe der Korrelation innerhalb der Zeitcluster der Fehlervariable
#' (Anteil der Varianz des Zeiteffekts an der gesamten Varianz der Fehlervariable)
#' @returns Vektor mit Länge 6. Durchschnittliche Schätzung der Standardfehler nach OLS, Fama-MacBeth, Cluster und Newey-West sowie
#' die durchschnittlichen wahren Standardfehler der OLS Regression und der Fama-MacBeth Regression
#' @importFrom plm plm
#' @export
Zeiteffekt <- function(N = 500, T = 10, Anzahl_Regressoren = 1,
                       Anz_Sim = 100, Korr_delta = 0.25){

  beta <- 1

  # Varianz und Korrelation der erklärenden Variablen und der Fehlervariable

  for(i in 1:Anzahl_Regressoren){
    assign(paste0("Korr_zeta",i),
           as.double(readline(prompt = paste("Legen Sie für Regressor", paste0("X",i) , "die Höhe der Korrelation zwischen Variablen der gleichen Periode fest ")))
    )
  }
  for(r in 1:Anzahl_Regressoren){
    assign(paste0("Var_X", r), 1^2)
    assign("x", get(paste0("Korr_zeta",r)))
    assign("y", get(paste0("Var_X",r)))
    assign(paste0("Var_zeta", r), x*y)
  }

  Var_eta <- 2^2
  Var_delta <- Korr_delta*Var_eta


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
        Cov_Mat_Fama[i,j] <- 0
      }
    }
  }

  # Vorbereitung Simulation

  Ergeb_Mat <- vector(mode = "list", length = Anzahl_Regressoren)
  for(i in 1:Anzahl_Regressoren){
    x <- matrix(rep(0, times = Anz_Sim*6), ncol = 6)
    colnames(x) <- c("OLS", "Fama-MacBeth", "Cluster", "Newey-West",
                     "wahrer Stdfehler OLS", "wahrer Stdfehler Fama")
    names(Ergeb_Mat)[[i]] <- paste0("Ergeb_Mat_b", i)
    Ergeb_Mat[[i]] <- x
  }


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    for(r in 1:Anzahl_Regressoren){
      assign(paste0("Varzeta"), get(paste0("Var_zeta",r)))
      zeta <- stats::rnorm(T, mean = 0, sd = sqrt(Varzeta))
      zeta_neu <- vector(mode = "double", length = N*T)
      for(i in 1:N){
        for(t in 1:T){
          zeta_neu[(i-1)*T+t] <- zeta[t]
        }
      }
      assign(paste0("VarX"), get(paste0("Var_X",r)))
      nu <- stats::rnorm(N*T, mean = 0, sd = sqrt(VarX-Varzeta))

      assign(paste0("X", r), zeta_neu+nu)
    }

    delta <- stats::rnorm(T, mean = 0, sd = sqrt(Var_delta))
    delta_neu <- vector(mode = "double", length = N*T)
    for(i in 1:N){
      for(t in 1:T){
        delta_neu[(i-1)*T+t] <- delta[t]
      }
    }
    epsilon <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_eta-Var_delta))
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

    daten <- data.frame(X_Mat, unternehmen, jahr, y)
    daten_colnames <- vector(mode = "character", length = Anzahl_Regressoren + 3)
    daten_colnames[(Anzahl_Regressoren+1):(Anzahl_Regressoren+3)] <- c("unternehmen", "jahr", "y")
    for(i in 1:Anzahl_Regressoren){
      daten_colnames[i] <- paste0("X",i)
    }
    colnames(daten) <- daten_colnames
    paneldaten <- plm::pdata.frame(daten, index = c("unternehmen","jahr"))

    RegGleich <- paste0("y ~ X1")
    if(Anzahl_Regressoren > 1){
      for(i in 2:Anzahl_Regressoren){
        assign(paste0("RegGleich"), paste0(RegGleich, paste0(" + X",i)))
      }
    }

    Reg <- plm::plm(stats::as.formula(RegGleich), data = paneldaten, model = "pooling")

    # Standardfehler berechnen

    # OLS Standardfehler

    sd_OLS <- lmtest::coeftest(Reg)[2:(Anzahl_Regressoren + 1),2]
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,1] <- sd_OLS[i]
    }

    # Fama-MacBeth Standardfehler
    paneldaten_Fama <- plm::pdata.frame(daten, index = c("jahr", "unternehmen"))

    Reg_Fama <- plm::pmg(stats::as.formula(RegGleich), data = paneldaten_Fama, model = "mg")

    sd_Fama <- lmtest::coeftest(Reg_Fama)[2:(Anzahl_Regressoren + 1),2]
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,2] <- sd_Fama[i]
    }

    # Cluster Standardfehler

    VarMat_Cluster <- plm::vcovHC(Reg, cluster = "time", type = "sss")
    Var_Cluster <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 2:(Anzahl_Regressoren + 1)){
      Var_Cluster[i-1] <- VarMat_Cluster[i,i]
    }
    sd_Cluster <- sqrt(Var_Cluster)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,3] <- sd_Cluster[i]
    }

    # Newey-West Standardfehler

    VarMat_Newey <- plm::vcovNW(Reg, maxlag = (T-1))
    Var_Newey <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 2:(Anzahl_Regressoren + 1)){
      Var_Newey[i-1] <- VarMat_Newey[i,i]
    }
    sd_Newey <- sqrt(Var_Newey)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,4] <- sd_Newey[i]
    }


    # wahre Varianz
    # Var(b|X)=(X'X)^{-1}X'\Omega X(X'X)^{-1}
    intercept <- as.matrix(rep(1, times = N*T))
    X_Mat <- cbind(intercept, X_Mat)

    Varianzmatrix <- solve(t(X_Mat)%*%X_Mat)%*%t(X_Mat)%*%Cov_Matrix%*%X_Mat%*%solve(t(X_Mat)%*%X_Mat)
    wahreVar <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 2:(Anzahl_Regressoren + 1)){
      wahreVar[i-1] <- Varianzmatrix[i,i]
    }
    wahresd <- sqrt(wahreVar)
    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,5] <- wahresd[i]
    }


    # wahre Varianz Fama

    X_Mat_Fama <- cbind(intercept, as.matrix(
      as.data.frame(paneldaten_Fama[,(1:Anzahl_Regressoren)])))

    X_Mat_t <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_t[[t]] <- X_Mat_Fama[(N*(t-1)+1):(N*t),]
    }
    X_Mat_sum <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_sum[[t]] <- (Var_eta-Var_delta)*solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])
    }

    Varianzmatrix_Fama <- (1/(T^2))*(Reduce("+", X_Mat_sum))

    wahreVar_Fama <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 2:(Anzahl_Regressoren + 1)){
      wahreVar_Fama[i-1] <- Varianzmatrix_Fama[i,i]
    }
    wahresd_Fama <- sqrt(wahreVar_Fama)

    for(i in 1:Anzahl_Regressoren){
      Ergeb_Mat[[i]][D,6] <- wahresd_Fama[i]
    }
  }


  for(i in 1:Anzahl_Regressoren){
    x <- Ergeb_Mat[[i]]
    Ergeb_Mat[[i]] <- colMeans(x)
  }

  return(Ergeb_Mat)

}
