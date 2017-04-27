################ FUN??ES AUXILIARES CMIS FORECAST ################

# Require packages
#Rpacks()

# Tratar erros
trycatch_w_e <- function(expr){
	W <- NULL
	w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
	}
	list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
}

# Transforma??o de BoxCox
BoxCox <- function (x, lambda) {
	if (lambda < 0)
		x[x < 0] <- NA
	if (lambda == 0)
		out <- log(x)
	else out <- (sign(x) * abs(x)^lambda - 1)/lambda
	if (!is.null(colnames(x)))
		colnames(out) <- colnames(x)
	return(out)
}

# Transforma??o de BoxCox inversa
InvBoxCox <- function (x, lambda) {
	if (lambda < 0)
		x[x > -1/lambda] <- NA
	if (lambda == 0)
		out <- exp(x)
	else {
		xx <- x * lambda + 1
		out <- sign(xx) * abs(xx)^(1/lambda)
	}
	if (!is.null(colnames(x)))
		colnames(out) <- colnames(x)
	return(out)
}

BoxCox.lambda <- function (x, method = c("guerrero", "loglik"), lower = -1, upper = 2) {
	if (any(x <= 0))
		lower <- 0
	if (length(x) <= 2 * frequency(x))
		return(1)
	method <- match.arg(method)
	if (method == "loglik")
		return(bcloglik(x, lower, upper))
	else return(guerrero(x, lower, upper))
}

#Busca e trata aoutliers: se a m?trica ? de tend?ncia troca outliers por interpola??o linear, se ela for peri?dica troca outliers por valores ajustados por Seasonal Decomposition of Time Series by Loess (STL). LOESS se refere ? interpola??o polinomial para capturar melhor os picos e vales fazendo quebra dos dados em Tend?ncia, Sazonalidade e Restante
tsoutliers <- function (x, iterate = 2, lambda = NULL) {
	missng <- is.na(x)
	n <- length(x)
	freq <- frequency(x)
	if (freq > 1 & n > 2 * freq) {
		xx <- na.interp(x, lambda = lambda)
		if (!is.null(lambda)) {
			xx <- BoxCox(xx, lambda = lambda)
		}
		fit <- stl(xx, s.window = "periodic", robust = TRUE)
		resid <- fit$time.series[, "remainder"]
		resid[missng] <- NA
	}
	else {
		tt <- 1:n
		if (!is.null(lambda)) {
			xx <- BoxCox(x, lambda = lambda)
		}
		else xx <- x
		mod <- loess(xx ~ tt, na.action = na.exclude, family = "symmetric",
			degree = 1, span = min(20/n, 0.75))
		resid <- xx - fitted(mod)
	}
	resid.q <- quantile(resid, prob = c(0.025, 0.975), na.rm = TRUE)
	#resid.q <- quantile(resid, prob = c(0.1, 0.9), na.rm = TRUE)
	iqr <- diff(resid.q)
	limits <- resid.q + 2 * iqr * c(-1, 1)
	outliers <- which((resid < limits[1]) | (resid > limits[2]))
	x[outliers] <- NA
	x <- na.interp(x, lambda = lambda)
	if (iterate > 1 & frequency(x) <= 1) {
		tmp <- tsoutliers(x, iterate = 1, lambda = lambda)
		if (length(tmp$index) > 0) {
			outliers <- sort(c(outliers, tmp$index))
			x[outliers] <- NA
			x <- na.interp(x, lambda = lambda)
		}
	}
	return(list(index = outliers, replacements = x[outliers]))
}

# Utiliza a fun??o tsoutliers para "Limpar" os outliers e retornar o vetor de dados tratados
tsclean <- function (x, replace.missing = TRUE, lambda = NULL) {
	outliers <- tsoutliers(x, lambda = lambda)
	x[outliers$index] <- outliers$replacements
	if (replace.missing)
		x <- na.interp(x, lambda = lambda)
	return(x)
}

# Fun??o para fazer interpola??o pros casos onde h? dados missing na s?rie
na.interp <- function (x, lambda = NULL) {
	missng <- is.na(x)
	if (sum(missng) == 0)
		return(x)
	if (is.null(tsp(x)))
		x <- ts(x)
	freq <- frequency(x)
	n <- length(x)
	tt <- 1:n
	idx <- tt[!missng]
	if (freq <= 1 | n <= 2 * freq) {
		if (!is.null(lambda)) {
			x[idx] <- BoxCox(x[idx], lambda = lambda)
			xx <- as.ts(approx(idx, x[idx], 1:n, rule = 2)$y)
			xx <- InvBoxCox(xx, lambda = lambda)
		}
		else xx <- as.ts(approx(idx, x[idx], 1:n, rule = 2)$y)
		tsp(xx) <- tsp(x)
		return(xx)
	}
	else {
		if (!is.null(lambda)) {
			x <- BoxCox(x, lambda = lambda)
		}
		X <- cbind(fourier(x, 3), poly(tt, degree = 3))
		fit <- lm(x ~ X, na.action = na.exclude)
		pred <- predict(fit, newdata = data.frame(X))
		x[missng] <- pred[missng]
		if (!is.null(dim(x)))
			stop("The time series is not univariate.")
		fit <- stl(x, s.window = 11, robust = TRUE)
		sa <- seasadj(fit)
		sa <- approx(idx, sa[idx], 1:n, rule = 2)$y
		x[missng] <- sa[missng] + fit$time.series[missng, "seasonal"]
		if (!is.null(lambda)) {
			x <- InvBoxCox(x, lambda = lambda)
		}
		return(x)
	}
}


## Função genérica para converter dados em séries temporais
#1) dados com data em forma factor ou character (data.frame) data(diario)
#2) dados em formato zoo e xts (data(MSFT))
#3) dados em formato ts e mts (data(MSFT))
#4) dados em forma numeric (simulações)

# Estratégia: 
# - Isolar a data e transformar em character
# - Isolar valores e transformar em numeric ou data.frame
# - Juntar e retornar em formato TimeSeries (dates, minute, hour, day, week, month, year e value)

# Exemplos
#require(RCapaTFC)
#require(forecast)
#require(xts)
#require(zoo)
#data(mensal, package="RCapaTFC")
#data(MSFT, package="timeSeries")
#daddf <- mensal[,1:2]
#daxts <- as.xts(MSFT[,1:3])
#dazoo <- as.zoo(daxts)
#damts <- ConvertDataToTs(daddf, "month")
#dasim <- rnorm(100)

#head(ConvertData(dasim, tsfrequency = "sec", dateformat='%d/%m/%Y %H:%M:%S', outType="df"))

ConvertData <- function(Data, tsfrequency = "day", dateformat='%d/%m/%Y %H:%M:%S', outType = "ts", tz = Sys.getenv("TZ"), ...) {

	## Check data, dates and frequencies
	if(is.null(Data)) stop("Empty dataset!\n")
	if(is.null(dateformat) && any(class(Data)=="data.frame")) stop("data.frame needs first column as date format. ex: '%Y/%m/%d', '%Y-%m-%d', '%Y/%m/%d %H:%M:%S', '%Y-%m-%d %H:%M:%S', etc.!\n")
	#check data frequancies
	tsfrequency <- match.arg(tsfrequency, c("year","month","day","hour","min","sec"))
	
	outType <- match.arg(outType, c("ts","xts","df"))
	
	# check ts frequency
	if (tsfrequency %in% c("year","month")){
	freq <- 12
	} else if (tsfrequency=="day") {
	freq <- 7
	} else if (tsfrequency=="hour") {
	freq <- 24
	} else if (tsfrequency=="min"){
	freq <- 60
	} else {
	freq <- 1
	}
	## Rules for converting process
	  if (any(class(Data)=="data.frame")){
		date  <- as.character(Data[,1])
		value <- data.frame(Data[,-1, drop=FALSE])
		out <- TimeSeries(date, dateformat, value)
	} else if(any(class(Data) %in% c("zoo","xts"))) {
		stopifnot(inherits(time(Data), "Date") || inherits(time(Data), "POSIXt"))
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		if (length(dim(Data)) < 2) {
			date <- time(Data)
			value <- data.frame(zoo::coredata(Data))
			names(value) <- deparse(substitute(Data))
			out <- TimeSeries(date, fmt, value)
		} else {
			date <- time(Data)+0.01
			value <- as.data.frame(zoo::coredata(Data))	
			if (length(dim(Data)) < 2) names(value) <- deparse(substitute(Data))
			out <- TimeSeries(date, fmt, value)
		}
	} else if(any(class(Data) %in% c("ts","mts"))) {
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		date  <- as.POSIXlt(as.Date(index(Data)))+0.01
		value <- data.frame(zoo::coredata(Data))
		if (length(dim(Data)) < 2) names(value) <- deparse(substitute(Data))
		out <- TimeSeries(date, fmt, value)
	} else if(any(class(Data) %in% c("numeric","integer"))){
		Data <- ts(Data)
		fmt <- if (inherits(time(Data), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
		date  <- as.POSIXlt(as.Date(index(Data)))+0.01
		value <- data.frame(value = zoo::coredata(Data))
		out <- TimeSeries(date, fmt, value)
	} else {
		stop("Check your dataset, dates and/or class!\n")
	}	
	if(outType=='ts') {
		O <- ts(value, start = Start(tsfrequency, out), frequency=freq)
	} else if (outType == "xts"){
		O <- xts::as.xts(zoo::zoo(x=out[-seq(7)], out$dates), tzone = tz)
	} else {
		O <- out
	}
	return(O)
}

## Tratar dados para a série temporal com datas no formato "%d/%m/%Y %H:%M:%S"
ts.dados <- function(dados_hist, n_passos_frente=20, freq=7, normalize = FALSE, ...){

	names(dados_hist) <- tolower(names(dados_hist))	
	tsfrequency <- if(freq == 7) "day" else if (freq == 12) "month" else if (freq == 24) "hour" else "min"
	
	#if(!class(dados_hist) %in% c("data.frame","xts")) {
	#	cat("Exigido dados no formato data-valor 'data.frame' ou 'xts'!\n")
	#}
	
	# Converte dados para séries
	dados_xts <- Try_error(ConvertData(dados_hist, tsfrequency, outType = "xts"))
	
	if (class(dados_xts)[1] !="try-error") {	
		# Normaliza os dados, se preciso
		if(normalize) dados_xts <- scale(dados_xts)
		
		# Calacula sequancia de pontos para forecast
		Horizonte <- ForecastHorizon(dados_xts, tsfrequency, horizon=n_passos_frente)
		
		dados_historicos <- data.frame(data=index(dados_xts), as.data.frame(dados_xts), stringsAsFactors = FALSE)	
		names(dados_historicos) <- c("data","realizado")
		rownames(dados_historicos) <- NULL
	} else {
		stop("Check data!")
	}
	out <- c()
	out$historico_novo_ts <- Horizonte[[1]]
	out$dados_historicos <- dados_historicos	
	out$n_dados_historicos <- nrow(dados_historicos)	
	out$data_previsao <- Horizonte[[2]]
	out$n_passos_frente <- n_passos_frente
	out$frequencia <- freq	
	return(out)
}


## Tratar dados para a s?rie temporal com datas no formato "%d/%m/%Y %H:%M:%S"
ts.dados.old <- function(dados_hist, n_passos_frente=20, freq=7, normalize = FALSE, ...){
	out <- c()

	names(dados_hist) <- tolower(names(dados_hist))
	dados <- as.numeric(dados_hist[,2])

	dados_ts <- TimeSeries(dados_hist[,1],  "%d/%m/%Y %H:%M:%S", dados_hist[,2])

	if (normalize == TRUE) {
		temp_da <- try(TimeSeries(dados_hist[,1],  "%d/%m/%Y %H:%M:%S", as.numeric(scale(dados_hist[,2]))))
		if (class(temp_da) != "try-error") {
			dados_ts <- temp_da
		} else {
			temp_da <- interpNA(dados_hist[,2], method = "linear")
			temp_da <- try(TimeSeries(dados_hist[,1],  "%d/%m/%Y %H:%M:%S", temp_da))
		}
		dados_ts <- temp_da
	}

	#definindo vari?veis globais
	assign("ano_histo", dados_ts$year[1], envir = .GlobalEnv)
	assign("mes_histo", dados_ts$month[1], envir = .GlobalEnv)
	assign("semana_histo", dados_ts$week[1], envir = .GlobalEnv)
	assign("dia_histo", dados_ts$day[1], envir = .GlobalEnv)
	assign("hora_histo", dados_ts$hour[1], envir = .GlobalEnv)
	assign("dia_da_semana_histo", as.numeric(format(as.Date(dados_ts$dates[1]),'%w'))+1, envir = .GlobalEnv)

	#definindo vari?veis internas para a fun??o
	ultimo_ano_histo <- dados_ts$year[length(dados)]
	ultimo_mes_histo <- dados_ts$month[length(dados)]
	ultimo_semana_histo <- dados_ts$week[length(dados)]
	ultimo_dia_histo <- dados_ts$day[length(dados)]
	ultimo_hora_histo <- dados_ts$hour[length(dados)]
	ultimo_dia_da_semana_histo <- as.numeric(format(as.Date(dados_ts$dates[length(dados)]),'%w'))+1

	#criando as datas das previsoes
	if (freq==12) {
	#datas de inicio e fim de previsao para freq = 12
		ultima_data_historica <- ISOdate(ultimo_ano_histo, ultimo_mes_histo, ultimo_dia_histo)
		ultima_data_historica <- as.Date(ultima_data_historica)

		primeiro_mes_previsao <- ultima_data_historica + 1
		ultimo_mes_previsao   <- ultima_data_historica + n_passos_frente		
		
		pri_mes_previsao      <- format(primeiro_mes_previsao,'%m')
		pri_ano_previsao      <- format(primeiro_mes_previsao,'%Y')
		primeiro_mes_previsao <- ISOdate(pri_ano_previsao,pri_mes_previsao,1,0)
		ult_mes_previsao      <- format(ultimo_mes_previsao,'%m')
		ult_ano_previsao      <- format(ultimo_mes_previsao,'%Y')
		ultimo_mes_previsao   <- ISOdate(ult_ano_previsao,ult_mes_previsao,1,0)
		data_previsao         <- seq(from = primeiro_mes_previsao, to=ultimo_mes_previsao,by = "month")
	} else if (freq==7) {
	#datas de inicio e fim de previsao para freq = 7
		ultima_data_historica <- ISOdate(ultimo_ano_histo,ultimo_mes_histo,ultimo_dia_histo)
		ultima_data_historica <- as.Date(ultima_data_historica)
		
		primeiro_dia_previsao <- ultima_data_historica + 1
		ultimo_dia_previsao   <- ultima_data_historica + n_passos_frente
		
		pri_dia_previsao      <- format(primeiro_dia_previsao,'%d')
		pri_mes_previsao      <- format(primeiro_dia_previsao,'%m')
		pri_ano_previsao      <- format(primeiro_dia_previsao,'%Y')
		primeiro_dia_previsao <- ISOdate(pri_ano_previsao,pri_mes_previsao,pri_dia_previsao,0)
		ult_dia_previsao      <- format(ultimo_dia_previsao,'%d')
		ult_mes_previsao      <- format(ultimo_dia_previsao,'%m')
		ult_ano_previsao      <- format(ultimo_dia_previsao,'%Y')
		ultimo_dia_previsao   <- ISOdate(ult_ano_previsao,ult_mes_previsao,ult_dia_previsao,0)
		data_previsao         <- seq(from = primeiro_dia_previsao, to=ultimo_dia_previsao,by = "day")
	} else if (freq==24){
	#datas de inicio e fim de previsao para freq = 24
		ultima_data_historica  <- ISOdate(ultimo_ano_histo, ultimo_mes_histo, ultimo_dia_histo, ultimo_hora_histo)
		primeira_hora_previsao <- ultima_data_historica+3600
		ultima_hora_previsao   <- ultima_data_historica+(n_passos_frente*3600)
		data_previsao          <- seq(from=primeira_hora_previsao, to=ultima_hora_previsao,by = "hour")
	} else {
		ultima_data_historica  <- ISOdate(ultimo_ano_histo, ultimo_mes_histo, ultimo_dia_histo, ultimo_hora_histo)
		primeira_hora_previsao <- ultima_data_historica
		ultima_hora_previsao   <- ultima_data_historica+n_passos_frente
		data_previsao          <- seq(from=primeira_hora_previsao, to=ultima_hora_previsao, by = 1)
	}


	# confere se a primeira e a ?ltima observa??o dos dados s?o NA
	if(paste(dados[1], collapse = "") == "NA")	{
		dados[1] <- mean(dados, na.rm=TRUE)
	}

	if(paste(dados[length(dados)], collapse = "") == "NA")	{
		dados[length(dados)] <- mean(dados, na.rm=TRUE)
	}

	#interpola linearmente valores aos dados missing.
	dados <- interpNA(dados, method = "linear")

	#transforma os dados historicos em series temporais
	if (freq==12) {
		# transformar dados em series temporais - meses
		dados <- (ts(dados, frequency=12,start=c(ano_histo, mes_histo)))
	} else if (freq==7) {
		# transformar dados em series temporais ## start = c(semana, dia da semana)## - dias
		dados  <- (ts(dados, frequency = 7, start = c(semana_histo, dia_da_semana_histo)))
	} else if (freq==24) {
		# transformar dados em series temporais ( dias com 24 horas)- horas
		dados <- (ts(dados, frequency = 24, start = c(1)))
	} else {
		# Qualquer situa??o diferente de 7, 12 e 24
		dados <- (ts(dados, frequency = 1, start = c(1)))
	}

	qtd_dados_historicos <- length(dados)

	historico_novo <- dados
	out$historico_novo_ts <- historico_novo[,1]
	dados_historicos  <- data.frame(dados_ts$dates, historico_novo, stringsAsFactors = FALSE)
	names(dados_historicos) <- c("data","realizado")
	out$dados_historicos <- dados_historicos
	out$n_dados_historicos <- qtd_dados_historicos
	out$data_previsao <- data_previsao
	out$n_passos_frente <- n_passos_frente
	out$frequencia <- freq
	return(out)
}

## Tratar os outliers para a s?rie temporal utilizando STL para m?tricas sazonais e lm para aquelas com tend?ncia apenas.
trata.outliers.ts <- function(dados.ts, ...) {
	obj <- dados.ts
	out <- c()
	historico_novo_ts <- obj$historico_novo_ts
	dados_historicos  <- obj$dados_historicos

	## exportar  os dados tratados para o banco de dados
	temp_outliers <-  trycatch_w_e(tsoutliers(historico_novo_ts))$value
	outliers <- dados_historicos[temp_outliers$index,]
	outliers$realizado <- temp_outliers$replacements

	if (length(temp_outliers$index) > 0) {
		dados.clean <- trycatch_w_e(tsclean(historico_novo_ts))$value
	} else {
		dados.clean  <- historico_novo_ts
	}

	out$historico_novo_ts  <- dados.clean
	out$outliers <- outliers
	out$dados_historicos <- dados_historicos
  #out$historico_novo_ts <- dados.clean
	out$data_previsao <- obj$data_previsao
	out$n_passos_frente <- obj$n_passos_frente
	out$frequencia <- obj$frequencia
	return(out)
}

analise_residuos <- function(residuos, historico_novo, modelo) {
out <- c()

	x <- as.numeric(historico_novo)
	r <- as.numeric(residuos)

	di <- abs(length(x)-length(r))
	if (di > 0) {
		x <- x[-c(1:di)]
	}

#### testes para independencia dos residuos
	independencia <- trycatch_w_e(Box.test(r, lag=10, type = "Ljung-Box"))$value
	# Teste para ver se a media tende a zero
	media_zero <- trycatch_w_e(t.test(r, alternative='two.sided', mu=0.0, conf.level=.95))$value
	# Teste para ver se os residuos sao ruido branco
	ruido_branco <- trycatch_w_e(LB.test(modelo, no.error=TRUE))$value

	# Teste para normalidade dos res?duos jarque-bera
	normalidade <- trycatch_w_e(jarque.bera.test(r))$value
	# Teste de heterocedasticidade dos res?duos p-valor >0,05 indica homocedasticidade
	homocedasticidade <- trycatch_w_e(bptest(r ~ x))$value

	# Teste de durbin-watson para autocorrelacao dos res?duos se dw~2 ? independente
	autocorrelacao <- trycatch_w_e(dwtest(r ~ x))$value

	## Quando o teste for siginificativo a 95% retorna valores diferentes de zero ##
	if (class(independencia$p.value) == "numeric")		{p0 <- as.numeric(independencia$p.value)} else {p0 <- NA}
	if (class(media_zero$p.value) == "numeric") 		{p1 <- as.numeric(media_zero$p.value)} else {p1 <- NA}
	if (class(ruido_branco$p.value) == "numeric") 		{p2 <- as.numeric(ruido_branco$p.value)} else {p2 <- NA}
	if (class(normalidade$p.value) == "numeric") 		{p3 <- as.numeric(normalidade$p.value)} else {p3 <- NA}
	if (class(homocedasticidade$p.value) == "numeric") 	{p4 <- as.numeric(homocedasticidade$p.value)} else {p4 <- NA}
	if (class(autocorrelacao$p.value) == "numeric") 	{p5 <- as.numeric(autocorrelacao$p.value)} else {p5 <- NA}

	df.peso <-c(ifelse(p0 > 0.05, 2.0, 0),
				ifelse(p1 > 0.05, 1.0, 0),
				ifelse(p2 > 0.05, 1.0, 0),
				ifelse(p3 > 0.05, 1.0, 0),
				ifelse(p4 > 0.05, 2.0, 0),
				ifelse(p5 > 0.05, 3.0, 0))

	df.pvalor <- c(p0, p1, p2, p3, p4, p5)

	analise_residual <- data.frame(df.peso, df.pvalor = round(df.pvalor, 5))
	rownames(analise_residual) <- c("independencia","media_zero","ruido_branco","normalidade","homocedasticidade","autocorrelacao")

	out$analise_residual <- analise_residual
	out$soma_peso <- colSums(analise_residual)[1]
	out$soma_peso_inv <- 1/(colSums(analise_residual)[1])
	return(out)
}

decisao_arima_ets_v2 <- function(previsao_teste) {
	if (tolower(class(previsao_teste))[1] == "forecast"){
		x <- previsao_teste$x
		f <- previsao_teste$fitted
		r <- previsao_teste$residuals
		m <- previsao_teste$model
		res  <- as.numeric(analise_residuos(r, x, m)$soma_peso_inv)

		gof <-  c(accuracy(m))
		out <- round(c(gof, res), 8)
		names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
	} else {
		#objeto vazio para saida
		out <- c(rep(NA, 8))
		names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
	}
	return(out)
}


decisao_holtwinters_v2 <- function(modelo_teste) {
	if (tolower(class(modelo_teste))[1] == "forecast"){
		modelo <- modelo_teste$model
		x  <- modelo$x
		f0 <- modelo$fitted
		m  <- modelo
		f1 <- f0[,1]
		r <- x-f1
		res  <- as.numeric(analise_residuos(r, x, m)$soma_peso_inv)

		gof <-  c(accuracy(modelo_teste))
		out <- round(c(gof, res), 8)
		names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
	} else {
		#objeto vazio para saida
		out <- c(rep(NA, 8))
		names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
	}
	return(out)
}

decisao_rneurais <- function(previsao_teste) {
		if (tolower(class(previsao_teste))[1] == "forecast"){
			x <- previsao_teste$x
			f <- previsao_teste$fitted
			r <- previsao_teste$residuals
			m <- previsao_teste
			res <- NA
			res  <- as.numeric(analise_residuos(r, x, m)$soma_peso_inv)

			gof <-  c(accuracy(m))
			out <- round(c(gof, res), 8)
			names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
		} else {
			#objeto vazio para saida
			out <- c(rep(NA, 8))
			names(out) <- c("me", "rmse", "mae", "mpe", "mape", "mase", "acf1", "residuo_testes")
		}
		return(out)
}

### Ajusta 5 tipos de modelos iniciais
ajusta_modelos_iniciais <- function(obj.dados, ...){
	out <- out1 <- out2 <- out3 <- out4 <- out5 <- out6 <- c()
	obj <- obj.dados
	dados_historicos <- obj$dados_historicos
	historico_novo   <- obj$historico_novo
	n_passos_frente  <- obj$n_passos_frente
	freq <- obj$frequencia
	###tentativa com modelos sem transforma??es nos dados

	##############################################################################
	# Modelo ETS
	##############################################################################
	temp_ets <- trycatch_w_e(ets(historico_novo, damped = TRUE, additive.only= TRUE))$value
	if (class(temp_ets)[1] == "ets") {
		modelo_teste_ets   <- temp_ets
		previsao_teste_ets <- trycatch_w_e(forecast.ets(modelo_teste_ets, h = n_passos_frente))$value
	} else { # Preenche o objeto com vazio
		modelo_teste_ets   <- "NA"
		previsao_teste_ets <- "NA"
	}	# Bondade do ajuste por ETS
		gof_ets <- trycatch_w_e(decisao_arima_ets_v2(previsao_teste_ets))$value
		out1$modelo_teste_ets   <- modelo_teste_ets
		out1$previsao_teste_ets <- previsao_teste_ets
		out1$testes_ets <- gof_ets

	##############################################################################
	# Modelo ARIMA(p, d, q)X(P, D, Q)
	##############################################################################
	temp   <- trycatch_w_e(auto.arima(historico_novo))$value
	if (tolower(class(temp))[1] == "arima"){
		modelo_teste_arima   <- temp
		previsao_teste_arima <- trycatch_w_e(forecast(modelo_teste_arima, h = n_passos_frente))$value
	} else {
	# Tenta transforma??o de BoxCox aos dados. O resultado volta na escala original
		lambda <- trycatch_w_e(BoxCox.lambda(historico_novo))$value
		temp0  <- trycatch_w_e(auto.arima(historico_novo, lambda = lambda))$value
		temp1  <- trycatch_w_e(forecast(temp0, h = n_passos_frente))$value

		if (tolower(class(temp0))[1] == "arima" & tolower(class(temp1))[1] == "forecast") {
			modelo_teste_arima   <- temp0
			previsao_teste_arima <- temp1
		} else {
			modelo_teste_arima   <- "NA"
			previsao_teste_arima <- "NA"
		}
	}	# Bondade do ajuste por ARIMA(p, d, q)X(P, D, Q)
		gof_arima  <- trycatch_w_e(decisao_arima_ets_v2(previsao_teste_arima))$value
		out2$modelo_teste_arima   <- modelo_teste_arima
		out2$previsao_teste_arima <- previsao_teste_arima
		out2$gof_arima <- gof_arima

	##############################################################################
	# Modelo Holt Winters sazonal
	##############################################################################
	temp_hw_sazonal <- trycatch_w_e(HoltWinters(historico_novo))$value
	if (tolower(class(temp_hw_sazonal))[1] == "holtwinters"){
		modelo_teste_HoltWinters_sazonalidade   <- temp_hw_sazonal
		previsao_teste_HoltWinters_sazonalidade <- trycatch_w_e(forecast(modelo_teste_HoltWinters_sazonalidade, h = n_passos_frente,prediction.interval = TRUE))$value
	} else {
		modelo_teste_HoltWinters_sazonalidade <- "NA"
		previsao_teste_HoltWinters_sazonalidade <- "NA"
	}
		gof_hw_sazonal  <- trycatch_w_e(decisao_holtwinters_v2(previsao_teste_HoltWinters_sazonalidade))$value
		out3$modelo_teste_HoltWinters_sazonalidade   <- modelo_teste_HoltWinters_sazonalidade
		out3$previsao_teste_HoltWinters_sazonalidade <- previsao_teste_HoltWinters_sazonalidade
		out3$gof_hw_sazonal <- gof_hw_sazonal

	##############################################################################
	# Modelo Holt Winters n?o sazonalidade
	##############################################################################
	tem_hw_n_sazonal <- trycatch_w_e(HoltWinters(historico_novo, gamma = FALSE))$value
	if (tolower(class(tem_hw_n_sazonal))[1] == "holtwinters"){
		modelo_teste_HoltWinters_sem_sazonalidade <- tem_hw_n_sazonal
		previsao_teste_HoltWinters_sem_sazonalidade <- trycatch_w_e(forecast(modelo_teste_HoltWinters_sem_sazonalidade, h = n_passos_frente,prediction.interval = TRUE))$value
	} else {
		modelo_teste_HoltWinters_sem_sazonalidade <- "NA"
		previsao_teste_HoltWinters_sem_sazonalidade <- "NA"
	}
	# Bondade do ajuste por Holt Winters n?o sazonalidade
	gof_hw_n_sazonal  <- trycatch_w_e(decisao_holtwinters_v2(previsao_teste_HoltWinters_sem_sazonalidade))$value
	out4$modelo_teste_HoltWinters_sem_sazonalidade   <- modelo_teste_HoltWinters_sem_sazonalidade
	out4$previsao_teste_HoltWinters_sem_sazonalidade <- previsao_teste_HoltWinters_sem_sazonalidade
	out4$gof_hw_n_sazonal <- gof_hw_n_sazonal

	##############################################################################
	# Modelo Holt Winters com alisamento exponencial
	##############################################################################
	temp_hw_es <- trycatch_w_e(HoltWinters(historico_novo, gamma = FALSE, beta = FALSE))$value
	if (tolower(class(temp_hw_es))[1] == "holtwinters"){
		modelo_teste_HoltWinters_Exponential_Smoothing <- temp_hw_es
		previsao_teste_HoltWinters_Exponential_Smoothing <- trycatch_w_e(forecast(modelo_teste_HoltWinters_Exponential_Smoothing, h=n_passos_frente,prediction.interval = TRUE))$value
	} else {
		modelo_teste_HoltWinters_Exponential_Smoothing <- "NA"
		previsao_teste_HoltWinters_Exponential_Smoothing <- "NA"
	}

	# Bondade do ajuste por Holt Winters com alisamento exponencial
	gof_hw_es  <- trycatch_w_e(decisao_holtwinters_v2(previsao_teste_HoltWinters_Exponential_Smoothing))$value
	out5$modelo_teste_HoltWinters_Exponential_Smoothing   <- modelo_teste_HoltWinters_Exponential_Smoothing
	out5$previsao_teste_HoltWinters_Exponential_Smoothing <- previsao_teste_HoltWinters_Exponential_Smoothing
	out5$gof_hw_es <- gof_hw_es

	## Se os dados n?o possuirem sazonalidade ou se a frequ?ncia for menor que o comprimento dos dados, aplicar alisamento exponencial nos dados do objeto HoltWinters Sazonal, pois n?o h? ajuste sazonais quando o modelo possui menos de dois per?odos (freq < 2)

	if (length(historico_novo) < 2*freq) {
		gof_hw_sazonal <- gof_hw_es
		out3$modelo_teste_HoltWinters_sazonalidade   <- modelo_teste_HoltWinters_Exponential_Smoothing
		out3$previsao_teste_HoltWinters_sazonalidade <- previsao_teste_HoltWinters_Exponential_Smoothing
		out3$gof_hw_sazonal <- gof_hw_sazonal
	}

	##############################################################################
	# Modelo Redes Neurais (N?o tem limites Superior e Inferior, apenas proje??o
	##############################################################################
	temp_rneurais <- "NA"

	#trycatch_w_e(nnetar(historico_novo, lambda=0))$value
	if (class(temp_rneurais)[1] == "nnetar") {
		modelo_teste_rneurais   <- temp_rneurais
		previsao_teste_rneurais <- trycatch_w_e(forecast.nnetar(modelo_teste_rneurais, h = n_passos_frente))$value
	} else { # Preenche o objeto com vazio
		modelo_teste_rneurais   <- "NA"
		previsao_teste_rneurais <- "NA"
	}	# Bondade do ajuste por ETS
		gof_rneurais <- trycatch_w_e(decisao_rneurais(previsao_teste_rneurais))$value
		out6$modelo_teste_rneurais   <- modelo_teste_rneurais
		out6$previsao_teste_rneurais <- previsao_teste_rneurais
		out6$testes_rneurais <- gof_rneurais

	gof_modelo <- cbind(c("arima","ets","hw_sazonal","hw_n_sazonal","hw_es","rneural"))
	gof <- data.frame(gof_modelo, rbind(gof_arima, gof_ets, gof_hw_sazonal,  gof_hw_n_sazonal, gof_hw_es, gof_rneurais))

	# OBS.: o indicador residuo_testes cont?m o resultado dos rankings de 8 tipos de testes para os res?duos. Como para este indicador maior ? melhor, levamos em conta o inverso, ou seja, 1/residuo_testes, assim ele fica adequado para o ranking com mape, mad e rmse que quanto menores melhores s?o.

	gof$rk_mape <- rank(gof[, 6, drop = FALSE], na.last = TRUE, ties.method = "first") # mape
	gof$rk_rmse <- rank(gof[, 3, drop = FALSE], na.last = TRUE, ties.method = "first") # rmse
	gof$rk_resi <- rank(gof[, 9, drop = FALSE], na.last = TRUE, ties.method = "first") # residuo_testes
	gof$rk_mae  <- rank(gof[, 4, drop = FALSE], na.last = TRUE, ties.method = "first") # mae
	gof$rk_mase <- rank(gof[, 5, drop = FALSE], na.last = TRUE, ties.method = "first") # mase

	# Soma de rankings
	gof$rk_soma <- rowSums(gof[, c("rk_mape", "rk_rmse", "rk_resi", "rk_mae", "rk_mase")], na.rm = TRUE)
	#gofout <- gof[, c("gof_modelo", "rk_mape", "rk_rmse", "rk_resi", "rk_mae", "rk_mase", "rk_soma"), TRUE]

	# Ordena pelos menores rankings
	gof <- gof[order(gof$rk_soma, gof$rk_mape, decreasing = FALSE, na.last=NA), ]

	escolha <- gof[1,1]

	if (escolha[1] == "ets") {
		modelo <- out1
	} else if  (escolha[1] == "arima") { modelo <- out2
	} else if  (escolha[1] == "hw_sazonal") { modelo <- out3
	} else if  (escolha[1] == "hw_n_sazonal") {	modelo <- out4
	} else if  (escolha[1] == "hw_es") {modelo <- out5
	} else if  (escolha[1] == "rneural") {modelo <- out6
	} else {
		stop("Ajuste impossivel com estes tipos de modelos, tentar outra opcao!")
	}

	out$modelo_escolhido <- modelo
	out$gof <- gof
	out$modelos_testados <- list(ets = out1, arima = out2, hw_sazonal = out3, hw_n_sazonal = out4, hw_es = out5, gof = gof, historico_novo = historico_novo, n_passos_frente = n_passos_frente, frequencia = obj$frequencia)
	out$dados_historicos <- dados_historicos
	out$data_previsao <- obj$data_previsao
	out$outliers <- obj$outliers
	out$escolha <- escolha
	return(out)
}

## Salva resultados do modelo escolhido
salva.resultados_v2 <- function(obj.ajuste, id_metrica_frequencia=0, primeiro_modelo=0, parametro_1=NA, parametro_2=NA, parametro_3=NA, parametro_4=NA, parametro_5=NA, parametro_6=NA, parametro_7=NA, modelo=0, foreplot=TRUE, ...) {

	out <- tipo_modelo_ativo_id <- c()

	obj <- obj.ajuste
	obj_model <- obj$modelo_escolhido
	model     <- obj_model[[1]]
	previsao  <- obj_model[[2]]
	data_previsao <- obj$data_previsao
	n_passos_frente <- obj$modelos_testados$n_passos_frente
	id_modelo_teste <- as.numeric(modelo)

	funplot <- function(obj){
		mtest <- obj$modelos_testados
		nmesc <- as.character(obj$escolha)
		nm1 <- names(mtest)[1:5]
		nmplo <- nm1[nm1 != nmesc]

		m0 <- model
		m1 <- mtest[nmplo[1]][[1]][[1]]
		m2 <- mtest[nmplo[2]][[1]][[1]]
		m3 <- mtest[nmplo[3]][[1]][[1]]
		m4 <- mtest[nmplo[4]][[1]][[1]]

		dpar <- par(no.readonly = TRUE) #Padr?o de tela gr?fica do R
		nf <- layout(matrix(c(1,1,2,3,4,5), ncol = 2, byrow = TRUE))
		par(nf, mar = c(3,3,2,1), las=1)

		if(class(m0)[1] =="character") plot(rep(1, n_passos_frente), type="l", main="Sem projecao!") else plot(forecast(m0, h = n_passos_frente))
		if(class(m1)[1] =="character") plot(rep(1, n_passos_frente), type="l", main="Sem projecao!") else plot(forecast(m1, h = n_passos_frente))
		if(class(m2)[1] =="character") plot(rep(1, n_passos_frente), type="l", main="Sem projecao!") else plot(forecast(m2, h = n_passos_frente))
		if(class(m3)[1] =="character") plot(rep(1, n_passos_frente), type="l", main="Sem projecao!") else plot(forecast(m3, h = n_passos_frente))
		if(class(m4)[1] =="character") plot(rep(1, n_passos_frente), type="l", main="Sem projecao!") else plot(forecast(m4, h = n_passos_frente))

		par(dpar)
	}

	if (foreplot == TRUE) {
		funplot(obj)
	}

	tipo_modelo_ativo <- tolower(class(model))
	siz <- length(data_previsao)

	# Grava data frame com proje??o mais limites de superior e inferior de 95% de confian?a

	if (tipo_modelo_ativo[1] == "ets") {
		resultados <- data.frame(
		data  = as.character(data_previsao[1:length(previsao$mean)]),
		valor = as.numeric(previsao$mean),
		upp   = previsao$upper[,c(2)],
		low   = previsao$lower[,c(2)])
	} else if (tipo_modelo_ativo[1] == "arima") {
		resultados <- data.frame(
		data  = as.character(data_previsao[1:length(previsao$mean)]),
		valor = as.numeric(previsao$mean),
		upp   = previsao$upper[,c(2)],
		low   = previsao$lower[,c(2)])
	} else if (tipo_modelo_ativo[1]=="holtwinters") {
		resultados <- data.frame(
		data  = as.character(data_previsao[1:length(previsao$mean)]),
		valor = as.numeric(previsao$mean),
		upp   = previsao$upper[,c(2)],
		low   = previsao$lower[,c(2)])
	} else if (tipo_modelo_ativo[1]=="nnetar"){
		resultados <- data.frame(
		data  = as.character(data_previsao),
		valor = as.numeric(previsao$mean),
		upp   = NA,
		low   = NA)
	} else {
		resultados <- data.frame(
		data  = as.character(NA),
		valor =  NA,
		upp   = NA,
		low   = NA)
	}

	if (primeiro_modelo == 0) {

		if(tipo_modelo_ativo[1] == "ets") tipo_modelo_ativo_id <- 1
		if(tipo_modelo_ativo[1] == "arima") tipo_modelo_ativo_id <- 2
		if(tipo_modelo_ativo[1] == "holtwinters") tipo_modelo_ativo_id <- 3
		if(tipo_modelo_ativo[1] == "nnetar") tipo_modelo_ativo_id <- 4

		teste_modelo_ativo <- FALSE

		if(tipo_modelo_ativo_id[1] == id_modelo_teste) {

			if (tipo_modelo_ativo[1] == "ets") {
				teste_modelo_ativo = parametro_1 == paste(model$components[c(1,2,3)], collapse="")
			}

			if (tipo_modelo_ativo[1] == "arima") {

				arima_tendencia_p    <- as.numeric(parametro_1)
				arima_tendencia_d    <- as.numeric(parametro_6)
				arima_tendencia_q    <- as.numeric(parametro_2)
				arima_sazonalidade_p <- as.numeric(parametro_3)
				arima_sazonalidade_d <- as.numeric(parametro_7)
				arima_sazonalidade_q <- as.numeric(parametro_4)
				arima_freq <- as.numeric(parametro_5)

				teste_modelo_ativo1 <- (arima_tendencia_p == model$arma[1])
				teste_modelo_ativo2 <- (arima_tendencia_d == model$arma[6])
				teste_modelo_ativo3 <- (arima_tendencia_q == model$arma[2])
				teste_modelo_ativo4 <- (arima_sazonalidade_p == model$arma[3])
				teste_modelo_ativo5 <- (arima_sazonalidade_d == model$arma[7])
				teste_modelo_ativo6 <- (arima_sazonalidade_q == model$arma[4])
				teste_modelo_ativo7 <- (arima_freq == model$arma[5])

				if((teste_modelo_ativo1 & teste_modelo_ativo2 & teste_modelo_ativo3 & teste_modelo_ativo4 & teste_modelo_ativo5 & teste_modelo_ativo6 & teste_modelo_ativo7) == TRUE) {
					teste_modelo_ativo = TRUE
				} else {
					teste_modelo_ativo=FALSE
				}
			}

			if (tipo_modelo_ativo[1] == "holtwinters") {
				HoltWinters_gamma <- parametro_1
				HoltWinters_beta  <- parametro_2
				HoltWinters_seazonal  <- paste(parametro_3, sep="")

				teste_modelo_ativo1 <- (HoltWinters_gamma == paste(model$gamma, sep=" "))
				teste_modelo_ativo2 <- (HoltWinters_beta == paste(model$beta, sep=" "))
				teste_modelo_ativo3 <- (HoltWinters_seazonal == paste(model$seasonal, sep=" "))

				if((teste_modelo_ativo1 & teste_modelo_ativo2 & teste_modelo_ativo3) == TRUE) {
					teste_modelo_ativo=TRUE
				} else {
					teste_modelo_ativo=FALSE
				}
			}
		}

		##parametros do modelo ETS
		if (tipo_modelo_ativo[1] == "ets") {
			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse = ""), componentes = paste(model$components[c(1,2,3)], collapse="")," "," "," "," "," "," ")
		}

		if (tipo_modelo_ativo[1] == "arima") {
		##parametros do modelo arima-
			parametro_arima  <- paste(model$arma, sep = " ")
			parametro_arima_1 <- parametro_arima[c(1)]
			parametro_arima_2 <- parametro_arima[c(2)]
			parametro_arima_3 <- parametro_arima[c(3)]
			parametro_arima_4 <- parametro_arima[c(4)]
			parametro_arima_5 <- parametro_arima[c(5)]
			parametro_arima_6 <- parametro_arima[c(6)]
			parametro_arima_7 <- parametro_arima[c(7)]

			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse=""), parametro_arima_1, parametro_arima_2, parametro_arima_3, parametro_arima_4, parametro_arima_5, parametro_arima_6, parametro_arima_7)

		}
		if (tipo_modelo_ativo[1] == "holtwinters") {
		##parametros do modelo HoltWinters-
			parametro_gamma   <- paste(model$gamma, sep=" ")
			parametro_beta    <- paste(model$beta, sep=" ")
			parametro_sazonal <- paste(model$seasonal, sep=" ")
			##parametros_HoltWinters=data.frame(gamma=parametro_gamma,beta=parametro_beta,seasonal=parametro_sazonal)

			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse=""), parametro_gamma, parametro_beta,parametro_sazonal," "," "," "," ")

		}
	#ore.save(par_modelo, name=paste("ds_metrica_mod_", id_metrica_frequencia, sep=''), overwrite = TRUE)
	} else {

	## salvar o modelo ativo na tabela. somente quando mudar o modelo
		if (tipo_modelo_ativo[1] == "ets") {
		##parametros do modelo ETS
			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse=""), componentes = paste(model$components[c(1,2,3)], collapse="")," "," "," "," "," "," ")
		}

		if (tipo_modelo_ativo[1] == "arima") {
		##parametros do modelo arima-
			parametro_arima   <- paste(model$arma, sep=" ")
			parametro_arima_1 <- parametro_arima[c(1)]
			parametro_arima_2 <- parametro_arima[c(2)]
			parametro_arima_3 <- parametro_arima[c(3)]
			parametro_arima_4 <- parametro_arima[c(4)]
			parametro_arima_5 <- parametro_arima[c(5)]
			parametro_arima_6 <- parametro_arima[c(6)]
			parametro_arima_7 <- parametro_arima[c(7)]

			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse = ""), parametro_arima_1, parametro_arima_2, parametro_arima_3, parametro_arima_4, parametro_arima_5, parametro_arima_6, parametro_arima_7)
		}

		if (tipo_modelo_ativo[1]=="holtwinters") {
		##parametros do modelo HoltWinters-
			parametro_gamma <- paste(model$gamma, sep=" ")
			parametro_beta  <- paste(model$beta, sep=" ")
			parametro_sazonal <- paste(model$seasonal, sep=" ")

			par_modelo <- data.frame(modelo = paste(model$call[c(1)], collapse=""), parametro_gamma, parametro_beta,parametro_sazonal," "," "," "," ")
		}
	}

	out$modelo     <- model
	out$par_modelo <- par_modelo
	out$bondade    <- obj$gof
	out$previsao   <- previsao
	out$projecao   <- resultados
	out$dados_historicos <- obj$dados_historicos
	out$outliers   <- obj$outliers
	out$modelos_testados <- obj$modelos_testados
	return(out)
}

cmis_tendencia_npar <-
function(x, metodo = "stl", plot = TRUE, ...) {
## Faz extração dq tendência utilizando decomposição em STL ou Médias Móveis e depois Calcula a direção desta tendência utilizando O Slope de Sen e O teste de Mann-Kendal que são testes não paramétricos.

	## Verifica se a métrica é constante
	if (identico(x) < 1) {
		#Caso haja valores missing: Interpolação linear!
		if (any(is.na(x))) {
			x <- interpNA(x)
		}
		
		if(class(x) != "ts") {
			cat("Variavel nao possui classe 'ts'. Utilizando cycle = 2 e frequency = 1!\n")
			x <- ts(x, frequency=2)
		}

		#lo <- lowess(time(x), x, iter = 5); try(plot(lo))
		#try(lines(lo, lwd=3, col=3))
		
		## Extrai a tendência para aplicar os testes
		## Utilizando decomposição da tendência e sazonalidade por médias móveis
		#da <- decompose(x)$trend
		
		if (metodo == "stl") {
			tsmod <- try(suppressMessages(stlm(x)))
			if (class(tsmod) !="try-error") {
				#cat("Extracao da tendencia por STL!\n")
				da <- tsmod$stl$time.series[,2]
			}
			if (plot) plot(tsmod$stl)
		} else if (metodo == "mm"){
			#cat("Problema com STL, extracao da temdencia por Medias Moveis!\n")
			tsmod <- try(suppressMessages(decompose(x)))
			if (class(tsmod) !="try-error") {
				#cat("Extracao da tendencia por Medias Moveis!\n")
				da <- tsmod$trend
			}
			if (plot) plot(tsmod)
		} else {
			stop(cat("O Metodo", metodo, "nao foi reconhecido. Utilize 'stl' ou 'mm'!\n"))
		}			
		
		#mm <- MannKendall(da);mm
		rk <- rkt(time(da), da)

		# Slope diferente de zero, testar os valroes de \tau
		if (round(rk$B, 6) != 0) {
			if (rk$tau > 0.10) {    
				if (rk$tau <= 0.95) sinal <- 1L
				else sinal <- 2L
			} else if (rk$tau < -0.10){  
				if (rk$tau > -0.95) sinal <- -1L
				else sinal <- -2L
			} else {
				sinal <- 0L
			}
		} else {
			sinal <- 0L
		}
		out <- c(slope = rk$B, tau = rk$tau, score = rk$S, p.valor = rk$sl, sinal=sinal)
		round(out, 4)
	}
}

cmis_tendencia_utilizacao <-
function (metrica, thx = 80, forma = "tcm", plotm = FALSE) {
	
	options(digits = 8)
	temp <- tc_m <- out <- c()

	#Trata tamanho dos dados
	if (length(metrica) < 6) {stop(paste("Pelo menos 6 observacoes sao necessarias! \n"))}

	# Trata outliers considerando sazonalidade e tendência (pacote forecast)
	temp <- try(tsclean(metrica))
	if (class(temp)!="try-error") metrica <- temp

	# Verifica presença de NA's
	if (any(is.na(metrica))) {
	  cat("Localizados ", sum(is.na(metrica)), " valores nulos na metrica. Feita remocao! \n")
	  metrica <- na.omit(metrica)
	}
	

	## Indicador de Crescimento usando médias 
	
	X <- metrica
	n <- length(X)
	p1<-X[1]
	pn<-X[n]
	tc <- tu1 <- tu2 <- c()
	
	# tcme: Indicador taxa de crescimento médio na forma exponencial corrigida
	if (forma == "tcme") {
		tc_m <- exp(log(pn/p1)/n)-1
		tcprod <- NA
		tc <- NA
	}

	# Indicador taxa de crescimento médio usando a taxa de crescimento a cada dois pontos dos dados
	# tc = (Xt-Xt-1)/Xt-1
	# tcm = raiz(t-1){(1+tc1)*(1+tc2)*...*(1+tct)}-1
	if (forma == "tcm"){

		for(i in 2:(length(X))){
		tc[i] <- ((X[i]-X[i-1])/X[i-1] + 0.000009) # Correçao de erro
		tu1[i] <- ((X[i]-thx)/thx + 0.000009) 
		}
		avgtutil1 <- (1 - mean(abs(tu1), na.rm = TRUE))
		#Tratar -Inf e +Inf
		is.na(tc) <- do.call(cbind, lapply(tc, is.infinite)) #Rápida para dados grandes.
		tcprod <- prod(1+tc, na.rm=TRUE)
		tc_m   <- sign(tcprod) * (abs(tcprod)^(1/(length(tc)-1)))-1  
	}

	if (tc_m >= 0){
		if (tc_m <= 0.001) {sinal = 0}
		else if (tc_m <= 0.70){sinal = 1}
		else {sinal = 2}
	}
	
	if (tc_m < 0) {
		if (tc_m >= -0.001) {
			sinal = 0
		} else if (tc_m >= -0.70){
			sinal = -1
		} else {
			sinal = -2
		}
	}

	# Cálculo do indicador de utilização: Tomar os 50% de pontos mais recentes da média móvel exponencial e tirar destes as taxas de utilização em relação ao treshold
	
	Y <- suppressWarnings(na.omit(WMA(X, n = 5))) ## Garante maiores pesos para valores mais recentes
	le <- length(Y)
	if (le < 24) {
		ce <- ceiling(le/2) #Os 50% mais recentes
	} else {
		ce <- ceiling(le*(3/4)) # Os 25% mais recentes
	}
	
	for(i in ce:le){
		tu2[i] <- ((Y[i]-thx)/thx + 0.000009) # Taxa de utilização em #relação ao theshold
	}
	
	avgtutil2 <- (1 - mean(abs(tu2), na.rm = TRUE))

	if(plotm) {
		dev.new(width=18, height = 6)
		par(mar = c(4,4,2,1), bg = colors()[350])
		plot(X, type="h", col=3, lwd=2, ylab=names(metrica), xlab='Tempo', main = 'Metrica vs Media Movel Ponderada (Maiores pesos obs mais recentes)');
		lines(suppressWarnings(WMA(na.omit(X), n=5)), type="l", col=4, lwd=2)
		legend('topright', legend=c('Métrica','WMA'), col=3:4,lwd=2)
		par(las = 0)
	}
	
	out <- data.frame(
		indicador_tendencia = sinal,
		tx_cresc_medio = tc_m,
		#tx_utilizacao_total_simples = avgtutil1,
		#vl_utilizacao_med_simples = mean(X, na.rm = TRUE),
		tx_utilizacao_total_ponderada = avgtutil2,
		vl_utilizacao_med_ponderada = mean(Y[ce:le], na.rm = TRUE))
	return(out)
}


