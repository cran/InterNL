% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InterNL.R
\name{InterNL}
\alias{InterNL}
\title{Time Series Intervention Model Using Non-linear Function}
\usage{
InterNL(Data, Time, TSModel, TSOrder = NULL, NLModel, InitialNLM)
}
\arguments{
\item{Data}{Time series data}

\item{Time}{Point of intervention}

\item{TSModel}{Time series model ("arima" or "ann")}

\item{TSOrder}{If model is ANN, then order is lag of the model}

\item{NLModel}{Non-linear models ("gompertz","logistic", "monomolecular", "richard", "hoerl")}

\item{InitialNLM}{Initial value for parameters of non-linear model}
}
\value{
\itemize{
  \item Accuracy: Accuracy metric of the proposed model
  \item PreFitted: Fitted values for the pre intervention series
  \item PostFitted: Prediction for the post intervention series
  \item NLM: Details of fitted non-linear model
}
}
\description{
Time Series Intervention Model Using Non-linear Function
}
\examples{
library("InterNL")
data<- as.ts(rnorm(120,100,50))
Result <- InterNL(Data = data,Time = 90, TSModel = "arima",
TSOrder=NULL, NLModel=NULL, InitialNLM=NULL )
}
\references{
\itemize{
\item Paul, R.K. and Yeasin, M., 2022. COVID-19 and prices of pulses in Major markets of India: Impact of nationwide lockdown. Plos one, 17(8), p.e0272999.
\item Yeasin, M., Paul, R.K., Das, S., Deka, D. and Karak, T., 2023. Change in the air due to the coronavirus outbreak in four major cities of India: What do the statistics say?. Journal of Hazardous Materials Advances, 10, p.100325.
}
}
