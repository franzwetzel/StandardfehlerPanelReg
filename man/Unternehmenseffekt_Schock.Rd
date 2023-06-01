% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Unternehmenseffekt_Schock.R
\name{Unternehmenseffekt_Schock}
\alias{Unternehmenseffekt_Schock}
\title{Standardfehler bei einem Unternehmenseffekt und einem Schock}
\usage{
Unternehmenseffekt_Schock(
  N = 500,
  T = 10,
  Anz_Sim = 100,
  Anteil_mu = 0.25,
  Anteil_u = 0.25,
  lag_X = 1,
  lag_eta = 1
)
}
\arguments{
\item{N}{N bestimmt die Anzahl der Firmen im Paneldatensatz}

\item{T}{N bestimmt die Anzahl der Jahre im Paneldatensatz}

\item{Anz_Sim}{Legt die Anzahl der Simulationen fest. Der Parameter bestimmt wie oft die wahren Standardfehler und die geschätzten für einen Paneldatensatz mit N Unternehmen und T Zeitperioden berechnet werden sollen.}

\item{Anzahl_Regressoren}{Legt die Anzahl der Regressoren fest}
}
\value{
Für jede Schätzung eines Betas die wahren Standardfehler für die OLS Regression und die Fama-MacBeth Regression
sowie die Schätzung der Standardfehler nach
}
\description{
Die Funktion berechnet für die angegebene Anzahl an Simulationen die wahren Standardabweichungen des OLS und Fama-MacBeth Schätzers und die geschätzten Standardabweichungen nach OLS, Fama-MacBeth, Cluster, Cluster plus Schock und Newey-West und bildet darüber den Durchschnitt.
}