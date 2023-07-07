#' Standardfehler bei einem Unternehmenseffekt und korrelierten Regressoren
#'
#' Die Funktion berechnet die geschätzten und wahren Standardfehler bei einem Unternehmenseffekt und korrelierten Regressoren.
#' @details Die Anzahl der Regressoren ist ein Argument der Funktion.
#' Nach dem Starten der Funtion wird der Nutzer aufgefordert, für jeden Regressor die
#' Korrelation innerhalb der Unternehmenscluster anzugeben.
#'
#' @param N N bestimmt die Anzahl der Firmen im Paneldatensatz
#' @param T T bestimmt die Anzahl der Perioden im Paneldatensatz
#' @param Anzahl_Regressoren Anzahl der Regressoren
#' @param Anz_Sim Anzahl der Simulationen. Anz_Sim bestimmt, wieviele Paneldatensätze
#' mit N Unternehmen und T Perioden simuliert werden sollen.
#' @param Anteil_u Hohe der Korrelation innerhalb der Unternehmenscluster der Fehlervariable
#' (Anteil der Varianz des Unternehmenseffekts an der gesamten Varianz der Fehlervariable)
#' @param Corr_nu Höhe der Korrelation zwischen den Komponenten der erklärenden Variablen,
#' die sich für jedes Unternehmen und jede Periode neu realisieren
#' @param Corr_mu Hohe der Korrelation zwischen den Unternehmenseffekten der Regressoren
#' @returns Vektor mit Länge 6. Durchschnittliche Schätzung der Standardfehler nach OLS, Fama-MacBeth, Cluster und Newey-West sowie
#' die durchschnittlichen wahren Standardfehler der OLS Regression und der Fama-MacBeth Regression
#' @importFrom plm plm
#' @export
Unternehmeneffekt_korrRegressoren <- function(N = 500, T = 10, Anzahl_Regressoren = 2,
                                              Anz_Sim = 100, Anteil_u = 0.25,
                                              Corr_nu = 0.3, Corr_mu = 0.3){

  beta <- 1

  for(i in 1:Anzahl_Regressoren){
    assign(paste0("Anteil_mu",i),
           as.double(readline(prompt = paste("Legen Sie für Regressor", paste0("X",i) , "die Höhe der Korrelation zwischen Variablen des gleichen Unternehmens fest ")))
    )
  }
  for(r in 1:Anzahl_Regressoren){
    assign(paste0("Var_X", r), 1^2)
    assign("x",get(paste0("Anteil_mu",r)))
    assign("y",get(paste0("Var_X",r)))
    assign(paste0("Var_mu", r), x*y)
  }

  Var_eta <- 2^2
  Var_u <- Anteil_u*Var_eta

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
          else if(i==j){
            Cov_Matrix[T*(i-1)+t,T*(j-1)+s] <- Var_u
          }
          else{
            Cov_Matrix[T*(i-1)+t,T*(j-1)+s] <- 0
          }
        }
      }
    }
  }

  # Fama
  Cov_Mat_Fama1 <- matrix(rep(0, times = N*N), ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      if(i==j){
        Cov_Mat_Fama1[i,j] <- Var_eta
      }
      else{
        Cov_Mat_Fama1[i,j] <- 0
      }
    }
  }
  Cov_Mat_Fama2 <- matrix(rep(0, times = N*N), ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      if(i==j){
        Cov_Mat_Fama2[i,j] <- Var_u
      }
      else{
        Cov_Mat_Fama2[i,j] <- 0
      }
    }
  }

  Corr_mu_Mat <- matrix(rep(0, times = Anzahl_Regressoren^2), ncol = Anzahl_Regressoren)
  for(i in 1:Anzahl_Regressoren){
    for(t in 1:Anzahl_Regressoren){
      if(i==t){
        Corr_mu_Mat[i,t] <- 1
      }
      else{
        Corr_mu_Mat[i,t] <- Corr_mu
      }
    }
  }

  Corr_nu_Mat <- matrix(rep(0, times = Anzahl_Regressoren^2), ncol = Anzahl_Regressoren)
  for(i in 1:Anzahl_Regressoren){
    for(t in 1:Anzahl_Regressoren){
      if(i==t){
        Corr_nu_Mat[i,t] <- 1
      }
      else{
        Corr_nu_Mat[i,t] <- Corr_nu
      }
    }
  }



  # Vorbereitung Simulation

  Ergeb_Mat <- vector(mode = "list", length = Anzahl_Regressoren)
  for(i in 1:Anzahl_Regressoren){
    x <- matrix(rep(0, times = Anz_Sim*6), ncol = 6)
    colnames(x) <- c("OLS", "Fama-MacBeth", "Cluster", "Newey-West",
                     "wahre Standardabweichung", "wahre Standardabweichung Fama")
    names(Ergeb_Mat)[[i]] <- paste0("Ergeb_Mat_b", i)
    Ergeb_Mat[[i]] <- x
  }


  # Simulation

  for(D in 1:Anz_Sim){
    # Simulation der Regressoren und der Fehlervariable

    # Simulation gemeinsame Komponente
    Var_X_Vektor1 <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      assign("x", get(paste0("Var_mu",i)))
      Var_X_Vektor1[i] <- x
    }

    mu_Mat <- as.matrix(faux::rnorm_multi(n = N, vars = Anzahl_Regressoren, mu = rep(0, times = Anzahl_Regressoren),
                                          sd = sqrt(Var_X_Vektor1), r = Corr_mu_Mat))

    mu_Mat_neu <- matrix(rep(0, times = Anzahl_Regressoren*N*T), ncol = Anzahl_Regressoren)
    for(i in 1:N){
      for(j in 1:T){
        mu_Mat_neu[(i-1)*T+j,] <- mu_Mat[i,]
      }
    }

    Var_X <- vector(mode = "double", length = Anzahl_Regressoren)
    for(i in 1:Anzahl_Regressoren){
      assign("x", get(paste0("Var_X",i)))
      Var_X[i] <- x
    }

    Var_X_Vektor2 <- Var_X-Var_X_Vektor1
    nu_Mat <- as.matrix(faux::rnorm_multi(n = N*T, vars = Anzahl_Regressoren, mu = rep(0, times = Anzahl_Regressoren),
                                          sd = sqrt(Var_X_Vektor2), r = Corr_nu_Mat))

    X_Mat <- mu_Mat_neu + nu_Mat


    u <- stats::rnorm(N, mean = 0, sd = sqrt(Var_u))
    u_neu <- vector(mode = "double", length = N*T)
    for(i in 1:N){
      for(t in 1:T){
        u_neu[(i-1)*T+t] <- u[i]
      }
    }
    epsilon <- stats::rnorm(N*T, mean = 0, sd = sqrt(Var_eta-Var_u))
    eta <- epsilon + u_neu


    # y Variable erstellen

    y <- eta
    for(i in 1:Anzahl_Regressoren){
      assign(paste0("y"), y+beta*X_Mat[,i])
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

    VarMat_Cluster <- plm::vcovHC(Reg)
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


    # wahre Varianz OLS
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

    X_Mat_Fama <- cbind(intercept,
                        as.matrix(as.data.frame(paneldaten_Fama[,1:Anzahl_Regressoren])))

    X_Mat_t <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_t[[t]] <- X_Mat_Fama[(N*(t-1)+1):(N*t),]
    }
    X_Mat_sum1 <- vector(mode = "list", length = T)
    for(t in 1:T){
      X_Mat_sum1[[t]] <- solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])%*%t(X_Mat_t[[t]])%*%Cov_Mat_Fama1%*%X_Mat_t[[t]]%*%solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])
    }
    X_Mat_sum2 <- vector(mode = "list", length = ((T-1)*T)/2)
    index <- 0
    for(t in 1:(T-1)){
      for(s in (t+1):T){
        index <- index + 1
        X_Mat_sum2[[index]] <- solve(t(X_Mat_t[[t]])%*%X_Mat_t[[t]])%*%t(X_Mat_t[[t]])%*%Cov_Mat_Fama2%*%X_Mat_t[[s]]%*%solve(t(X_Mat_t[[s]])%*%X_Mat_t[[s]])
      }
    }

    Varianzmatrix_Fama <- (1/(T^2))*(Reduce("+", X_Mat_sum1))+2*(1/(T^2))*(Reduce("+", X_Mat_sum2))

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
