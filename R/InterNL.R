#' @title Time Series Intervention Model Using Non-linear Function
#' @param Data Time series data
#' @param Time Point of intervention
#' @param TSModel Time series model ("arima" or "ann")
#' @param TSOrder If model is ANN, then order is lag of the model
#' @param NLModel Non-linear models ("gompertz","logistic", "monomolecular", "richard", "hoerl")
#' @param InitialNLM Initial value for parameters of non-linear model
#' @import stats forecast MLmetrics
#' @return
#' \itemize{
#'   \item Accuracy: Accuracy metric of the proposed model
#'   \item PreFitted: Fitted values for the pre intervention series
#'   \item PostFitted: Prediction for the post intervention series
#'   \item NLM: Details of fitted non-linear model
#' }
#' @export
#'
#' @examples
#' library("InterNL")
#' data<- as.ts(rnorm(120,100,50))
#' Result <- InterNL(Data = data,Time = 90, TSModel = "arima",
#' TSOrder=NULL, NLModel=NULL, InitialNLM=NULL )
#' @references
#' \itemize{
#' \item Paul, R.K. and Yeasin, M., 2022. COVID-19 and prices of pulses in Major markets of India: Impact of nationwide lockdown. Plos one, 17(8), p.e0272999.
#' \item Yeasin, M., Paul, R.K., Das, S., Deka, D. and Karak, T., 2023. Change in the air due to the coronavirus outbreak in four major cities of India: What do the statistics say?. Journal of Hazardous Materials Advances, 10, p.100325.
#' }

InterNL<-function(Data, Time, TSModel, TSOrder=NULL, NLModel, InitialNLM){
  data<-as.ts(Data)
  y1<-data[1:Time]
  y2<-data[-c(1:Time)]
  if(TSModel=="arima"){
    TModel<-auto.arima(y1)

  } else {
    TModel<-nnetar(y1,p=TSOrder)
  }
  Tfit<-na.omit(TModel$fitted)
  TFor<-as.vector(forecast(TModel,length(y2))$mean)

  resid<-y2-TFor
  t<-seq(1:length(resid))
  Gdata<-as.data.frame(cbind(resid,t))
  if(is.null(NLModel)){
    message("NLModel argument is missing. Random results has been generated for NLModel")
    GModel<-NULL
    gfit <-rnorm(length(resid), mean(resid), sd(resid)) ## only for example purpose
  } else if(NLModel=="gompertz"){
    fn<- resid ~ A * exp(-B * exp(-k * t))
    GModel<-nls(fn,start=InitialNLM, data = Gdata)
    gfit<-fitted(GModel)
  } else if(NLModel=="logistic"){
    fn<- resid ~ K / (1 + ((K - N0) / N0) * exp(-r * t))
    GModel<-nls(fn,start=InitialNLM, data = Gdata)
    gfit<-fitted(GModel)
  }else if(NLModel=="monomolecular"){
    fn<- resid ~ A * exp(-k * t)
    GModel<-nls(fn,start=InitialNLM, data = Gdata)
    gfit<-fitted(GModel)
  }else if(NLModel=="richard"){
    fn<- resid ~ A + (K - A) / (1 + exp(-B * (C - t)))^(1/beta)
    GModel<-nls(fn,start=InitialNLM, data = Gdata)
    gfit<-fitted(GModel)
  } else if(NLModel=="hoerl") {
    fn<- resid ~ (a*(b**t)*(t**c))
    GModel<-nls(fn,start=InitialNLM, data = Gdata)
    gfit<-fitted(GModel)
  }

  finalfor<-TFor+gfit
  ACC<-list(RMSE=RMSE(finalfor,y2), MAPE=MAPE(finalfor,y2))

  Result<-list(Accuracy=ACC,PreFitted=Tfit, PostFitted=finalfor, NLM=GModel)
  return(Result)
}
