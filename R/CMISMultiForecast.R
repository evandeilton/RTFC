cmisMultiforecast <- function(tsdata, Control=cmisControl(), fcMethod = NULL, ...) {

  #Controle
  h <- Control$maxHorizon
  level <- Control$level
  onlyfc  <- Control$onlyfc
  cvMethod <- Control$cvMethod
  tsfrequency <- Control$tsfrequency
  OutType <- Control$OutType
  OutlierClean <- Control$OutlierClean
  residlevel <- Control$residlevel
  
  outputFormat <- match.arg(Control$outputFormat, c("forecast","df", "both"))
  
  #if(ncol(as.data.frame(tsdata)) > 2) stop("Check data!\n")

  if(class(tsdata)[1] %in% c("data.frame")) {
    # To format the forecast output
    flag_forecast_format <- TRUE
    # Work with data and date formats
    x <- ConvertDataToTs(Data=tsdata, tsfrequency=tsfrequency, OutType = OutType, OutlierClean=OutlierClean)

    XTS <- ConvertDataToTs(Data=tsdata, tsfrequency=tsfrequency, OutlierClean=OutlierClean, OutType = "xts")

    # Dates in nice format for forecast output
    FCH <- ForecastHorizon(XTS, tsfrequency=tsfrequency, horizon=h)
    ForecastDates <- FCH$FCHorizon
    ForecastDates <- ForecastDates[!is.na(ForecastDates)]
	
	# Fix data.frame para formar a saída com os limites mais o forecast
	
	len <- length(x)+length(ForecastDates)	
	historico <- limInf <-  forecast <- limSup <- rep(NA, len)
	historico[1:length(x)] <- as.numeric(x)
	datas <- c(index(XTS), ForecastDates)
	
	#dd <- data.frame(data = datas, historico = historico, limInf =limInf,  forecast = forecast, limSup = limSup)	
	
  } else if (class(tsdata)[1] %in% c("xts","ts")){
    x <- tsdata
    ## Objeto com a projeção e datas em forma de data.frame para os melhores modelos obtidos
	FCH <- ForecastHorizon(x, tsfrequency=tsfrequency, horizon=h)
	ForecastDates <- FCH$FCHorizon
	ForecastDates <- ForecastDates[!is.na(ForecastDates)]
  } else if (class(tsdata)[1] == "numeric") {
    x <- ts(tsdata, frequency = 1)
    ## Objeto com a projeção e datas em forma de data.frame para os melhores modelos obtidos
	FCH <- ForecastHorizon(x, tsfrequency=tsfrequency, horizon=h)
	ForecastDates <- FCH$FCHorizon
	ForecastDates <- ForecastDates[!is.na(ForecastDates)]
  } else {
    stop("Check data!\n")
  }

  if(is.null(x) | length(as.numeric(x)) < 8) {
    stop(cat("Poucos dados", length(as.numeric(x)), "\n"))
  }

  ## Seleção do método de forecast baseado em tamanho da série, lineraidade, tendência e resíduo.
  
  MFc <- try(MultiForecast(x, fcMethod=fcMethod, Control=Control))

  ## Calcula estatística de bondade e resíduo e ordena os 5 primeiros modelos de melhor para o pior.
  
  CV_names <- try(BestModel(MFc, Control))
  CV_names <- CV_names[[1]]
  
	if(!all(is.na(CV_names)) & class(CV_names) != "try-error") { 
		CV_names <- as.list(CV_names)
		## Seleciona apenas os cinco melhores
		if (length(CV_names) > 5) CV_names <- CV_names[1:6]
	} else {
		# Arima e ETS
		CV_names <- list("auto.arimaForecast","etsForecast")
	}
    ## Objeto com forecast dos Melhores modelos
	best_models <- c()
	for(i in CV_names) {
		best_models[[i]] <- MFc[[i]]
	}
	
  if(outputFormat == "forecast") {
	return(structure(best_models, class = "multiforecast"))
  } else {
	if (outputFormat == "df" & class(tsdata)[1] %in% c("data.frame")) {
			ForecastDf <- plyr::llply(best_models, function(X) {
			forecast[(length(x)+1):len] <- as.numeric(X$mean)
			limSup[(length(x)+1):len] <- as.numeric(X$upper)
			limInf[(length(x)+1):len] <- as.numeric(X$lower)
			data.frame(data = datas, historico = historico, limInf = limInf, forecast = forecast, limSup = limSup)
		})
		return(structure(ForecastDf, class = "multiforecast"))
	} else {
		List <- list(ForecastDf = ForecastDf, Forecasts = best_models)
		return(structure(List, class = "multiforecast"))
	}
  }
}
