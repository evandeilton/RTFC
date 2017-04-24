cmisControl <- function (maxHorizon = 20, level=95, cvMethod = "MAPE",tsfrequency="day", OutType="ts", outputFormat = "forecast", OutlierClean=TRUE, onlyfc=FALSE, residlevel = 0.05,
userJDBC = "CMIS_OWNER", passJDBC = "CMIS_OWNER", driverJDBC = "ojdbc6.jar"

) {
  list(maxHorizon = maxHorizon, level = level, onlyfc = onlyfc, residlevel=residlevel, cvMethod=cvMethod, tsfrequency=tsfrequency, OutType=OutType, OutlierClean=OutlierClean, outputFormat=outputFormat, userJDBC=userJDBC, passJDBC=passJDBC, driverJDBC=driverJDBC)
}

forecastMethod <- function(x, fcMethod=NULL) {
  if(!all(is.null(fcMethod))){
  	metodo <- list('auto.arimaForecast','naiveForecast','etsForecast','rwForecast','lmForecast','stsForecast','HWsForecast','snaiveForecast','thetaForecast','HWnsForecast','HWesForecast','holtForecast')	
	  if(!class(fcMethod)[1] %in% c("character","list")) {
		stop("Indique os tipos de modelos em forma de lista.\nDica: list('stsForecast','auto.arimaForecast')")
	  } else {
		return(as.list(fcMethod))
	  }
	} else {
		if(is.null(x) | length(as.numeric(x)) < 8) {
			stop(cat("Poucos dados", length(as.numeric(x)), "\n"))
		}
		if (class(x)[1]!='ts' & class(x)[1] == "numeric"){
			x <- ts(x, frequency=1)
		}
		# série curta
		if(length(as.numeric(x)) < 2*max(cycle(x)) | max(cycle(x))==1) {
			short <- TRUE
		} else {
			short <- FALSE
		}
		if (!short) {# Só inicia testes se a série for longa
			# Objetos vazios para os resultados
			trend <- check_trend <- nlTests <- c()
		
			# Teste de tendência não paramétrico
			trend <- nparTrend(x)
			check_trend <- unname(trend["trend_sign"] != 0)
		
			# Teste de linearidade para a série
			nlTests <- list("terasvirta","white", "keenan", "tsay","tarTest")
			linear  <- Try_error(na.omit(sapply(nlTests,  function(n) linearityTest(x, n)$p.value)))
		
			if (class(linear)[1] == "numeric") {
				linear <- ifelse(sum(linear > 0.01) < 5, TRUE, FALSE)
			} else {
				linear <- TRUE
			}
		}
		## Decisão
		if (short) {
			## Modelos para Séries curtas ou sem frequencia definida
			metodo <- list("auto.arimaForecast","etsForecast","lmForecast", "holtForecast")
		} else if(linear & !check_trend) {
			## Modelos para Séries lineares mas sem tendência
			metodo <- list("naiveForecast","rwForecast","stsForecast","thetaForecast","HWnsForecast","HWesForecast", "HWsForecast")
		} else if(linear & check_trend) {
			## Modelos para Séries lineares e com tendência
			metodo <- list("auto.arimaForecast","etsForecast","HWsForecast","snaiveForecast", "holtForecast")
		} else {
			## Modelos para Séries não lineares com ou sem tendência
			metodo <- list("snaiveForecast","lmForecast","HWsForecast", "HWnsForecast")
		}
	}
  return(metodo)
}

MultiForecast <- function(x, fcMethod = NULL, Control = cmisControl()) {
  maxHorizon <- Control$maxHorizon
  level <- Control$level
  onlyfc <- Control$onlyfc

  metodo <- forecastMethod(x, fcMethod=fcMethod)

  Forecasts <- plyr::llply(metodo, function(X, ...) {
    temp <- Try_error(switch.cvforecast(x, nmodelo = X, h=maxHorizon, level=level, onlyfc=onlyfc))
    if (class(temp)[1]!="try-error") temp
    else NULL
  })

  if(length(Forecasts) == 0 | all(is.null(Forecasts))) {
    metodo <- list("auto.arimaForecast","etsForecast")
    Forecasts <- Try_error(switch.cvforecast(x, nmodelo = metodo, h=maxHorizon, level=level, onlyfc=onlyfc))
  }

  names(Forecasts) <- as.character(metodo)

  # Remove possíveis valores missing
  Forecasts <- Forecasts[!is.na(Forecasts)]

  # Remove possíveis valores missings
  Forecasts <- Filter(Negate(function(X) is.null(unlist(X))), Forecasts)
  return(structure(Forecasts, class = "multiforecast"))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Função para plot e estatísticas estilo ggplot2
plotFc <- function(df_dados, linesize=1, v_titgraf=NULL, dtfreq = 7, ...) { 
  stopifnot(require(ggplot2, quietly = TRUE))
  stopifnot(require(scales, quietly = TRUE))
  df_plotr  <- na.omit(df_dados[,c("data","historico")])
  df_plotp  <- na.omit(df_dados[,c("data","forecast","limInf","limSup")])
  
  tdformat <- if(nrow(df_plotr) == 90) "%m/%d/%y" else if(nrow(df_plotr) == 120) "%m/%d %H:%M:%S" else "%b/%Y"
  
  v_titgraf <- if(is.null(v_titgraf)) "Projecao" else v_titgraf  
  ggpl_grafico <- ggplot()
  ggpl_grafico + geom_line(data=df_plotr,aes(x=data,y = historico),color="blue",size=linesize) +
  geom_line(data=df_plotp,aes(x=data,y=forecast),color="darkgreen",size=linesize) +
  geom_line(data=df_plotp,aes(x=data,y=limInf),color="grey", size=linesize) +
  geom_line(data=df_plotp,aes(x=data,y=limSup),color="grey", size=linesize) + 
  scale_x_datetime(labels = date_format(tdformat))+
  xlab("Data") + 
  ylab("Forecast") + ggtitle(v_titgraf)+
  theme(axis.title = element_text(size = 9, angle = 0), plot.title = element_text(size = 11))  
}


plot.multiforecast <- function(obj, linesize = 1) {
  if(class(obj) != "multiforecast") stop("Objeto deve ter classe 'multiforecast'!\n")
  
  	linhas <- function(fit) {
		lines(fitted(fit), col='green')
		legend("topleft", c("observado","predito"), pch = '_', col = c("black","green"))
	}
  
  ## Sequencia com impar para gera??o dos gr?ficos
  ## Plot forecasts
  if (length(obj) == 6) {
    lfor <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 6)
  } else if (length(obj) == 5) {
    lfor <- c(1, 1, 2, 2, 3, 3, 4, 5)
  } else if (length(obj) == 4) {
    lfor <- c(1, 1, 2, 2, 3, 3, 4, 4)
  } else if (length(obj) == 3){
    lfor <- c(1, 1, 2, 2, 3, 3)
  } else if (length(obj) == 2){
    lfor <- c(1, 1, 2, 2)
  } else if (length(obj) == 1){
    lfor <- c(1, 1)
  }
  nf <- matrix(c(lfor), byrow=TRUE, ncol=2)   

  if(all(sapply(obj, class) == "data.frame")) {
	ggPlot <- lapply(1:length(obj), function(i) {
		Try_error(plotFc(obj[[i]], linesize=linesize, v_titgraf=names(obj)[i]))
	})
	
	multiplot(plotlist = ggPlot, layout = nf)

	} else {
	  par(layout(nf), mar = c(2.5, 5, 1.5, 1.5))
	  l_ply(obj, function(X) {
		if (class(X)[1] != "try-error" | !is.null(X))
		  try({
			plot(X, las=1)
			linhas(X)
			})
		}
	)
  }
}

summary.multiforecast <- function(obj, digits=5) {
	if(class(obj) != "multiforecast") stop("Objeto deve ter classe 'multiforecast'!\n")
	gof <- c()
	if(all(sapply(obj, class) == "data.frame")) {
		cat("Data frame contendo forecasts!\n")
		cat("Estatisticas descritivas: ===================================\n")
		descr <- llply(obj, function(X) descritiva(X[, -1, drop=FALSE], dig = digits))
		print(descr)
		cat("=============================================================\n")
		plot(obj)
		return(invisible(list(descr)))
	} else {
		cat("Lista contendo forecasts!\n")
		cat("\nEstatisticas de bondade do forecast: ===================================\n")
	  ## monta tabela com as estat
	  gof <- ldply(obj, function(X) {
		if (class(X)[1] != "try-error") {
		  st <- round(tsSummary(X), digits)
		  st[st %in% c(Inf, -Inf, NA, NULL)] <- 999999999
		  st
		} else 999999999 # Trata casos de erros com 999.999.999
	  })
	  print(gof)
	  cat("Testes sob os residuos: ===================================\n")
	  TR <- ldply(obj, function(X) Mresid(X))
	  print(TR)

	  cat("\nEstatisticas descritivas do forecast: ====================================\n")
	  descr <- llply(obj, function(X) descritiva(X, dig = digits))
	  print(descr)
	  cat("=============================================================\n")
	  plot(obj)
	  return(invisible(list(gof, TR, descr)))
	}
}


Mresid <- function(forecast) {
	out <- c()
	x <- as.numeric(forecast$x)
	r <- as.numeric(forecast$residuals)

	di <- abs(length(x)-length(r))
	if (di > 0) {
		x <- x[-c(1:di)]
	}

	# testes para independencia dos residuos
	independencia <- trycatch_w_e(Box.test(r, lag=10, type = "Ljung-Box"))$value
	# Teste para ver se a media tende a zero
	media_zero <- trycatch_w_e(t.test(r, alternative='two.sided', mu=0.0, conf.level=.95))$value
	# Teste para ver se os residuos sao ruido branco
	ruido_branco <- trycatch_w_e(LB.test(forecast, no.error=TRUE))$value

	# Teste para normalidade dos resíduos jarque-bera
	normalidade <- trycatch_w_e(jarque.bera.test(r))$value
	# Teste de heterocedasticidade dos resíduos p-valor >0,05 indica homocedasticidade
	homocedasticidade <- trycatch_w_e(bptest(r ~ x))$value

	# Teste de durbin-watson para autocorrelacao dos resíduos se dw~2 é independente
	autocorrelacao <- trycatch_w_e(dwtest(r ~ x))$value

	if (class(independencia$p.value) == "numeric")		{p0 <- as.numeric(independencia$p.value)} else {p0 <- NA}
	if (class(media_zero$p.value) == "numeric") 		{p1 <- as.numeric(media_zero$p.value)} else {p1 <- NA}
	if (class(ruido_branco$p.value) == "numeric") 		{p2 <- as.numeric(ruido_branco$p.value)} else {p2 <- NA}
	if (class(normalidade$p.value) == "numeric") 		{p3 <- as.numeric(normalidade$p.value)} else {p3 <- NA}
	if (class(homocedasticidade$p.value) == "numeric") 	{p4 <- as.numeric(homocedasticidade$p.value)} else {p4 <- NA}
	if (class(autocorrelacao$p.value) == "numeric") 	{p5 <- as.numeric(autocorrelacao$p.value)} else {p5 <- NA}

	df.pvalor <- round(c(p0, p1, p2, p3, p4, p5), 4)
	names(df.pvalor) <- c("independencia","media_zero","ruido_branco","normalidade","homocedasticidade","autocorrelacao")
	return(df.pvalor)
}

cmis_JDBC_Connect <- function(local_jdbcdriver, user, pass, JavaRAM = "-XmX2048m", setenv=FALSE, JAVA_HOME = NULL, autocommit = FALSE, ...) {

#local_jdbcdriver <- "C://Users//G0047743//Documents//sqldeveloper//jdbc//lib/ojdbc6.jar"
#user <- "CMIS_OWNER"
#pass <- "CMIS_OWNER"
#connect <- cmis_JDBC_Connect(local_jdbcdriver, user, pass)

	if(!require(RJDBC)) install.packages("RJDBC")
	
	if (setenv | !is.null(JAVA_HOME)) {
		cat("Definir variavel de ambiente do Java!\n")
		Sys.setenv(JAVA_HOME = JAVA_HOME)
	}
	if (autocommit) {
		cat("ATENCAO: autocomit ativado. Dados podem ser removidos ou salvos no banco automaticamante!\n")
		#set autocommit false - nao comita automaticamente no banco
		.jcall(jdbcConnection@jc,"V","setAutoCommit", autocommit)
	}

	options(parameters = JavaRAM)
	.jinit(force.init=TRUE)

	# Create connection driver and open connection
	jdbcDriver <- JDBC(driverClass="oracle.jdbc.OracleDriver", classPath=local_jdbcdriver)
	jdbcConnection <- try(dbConnect(jdbcDriver, "jdbc:oracle:thin:@//svuxhcap1:1521/HRCAP", user, pass))

	if(class(jdbcConnection) != "try-error") {
		return(jdbcConnection)
	} else {
		stop(cat("Problema na tentativa de conxao, reveja os argumentos!\n"))
	}
}


CMISForecastData <- function(dados=NULL, freq, idmini, idmfin, connect = FALSE, Control = cmisControl(), ...) {

	# 1-24	Hora
	# 2-12	Mensal
	# 3- 7	Diário

	if (!is.null(dados)) connect <- FALSE

	cat("-------------------------------------------------------------------\n")
	cat("Log: ETAPA 1 - Coleta e processamento dos dados.\n")
	cat("-------------------------------------------------------------------\n")
	#cat("Testa argumentos de entrada, conexao e coleta de dados!\n")

	if (connect & is.null(dados)) {
		Freq <- if(freq == 7) "'3'" else if (freq == 12) "'2'" else if (freq == 24) "'1'" else stop(paste("Exigido, '7=diario', '12=mensal', '24=horario'\n"))
		if(!is.numeric(idmini) | !is.numeric(idmfin)) stop("Numerico >= 0 exigido!\n")

		cat("Conexao com banco de dados exigida!\n")

		userJDBC   <- Control$userJDBC
		passJDBC   <- Control$passJDBC
		driverJDBC <- Control$driverJDBC

		cat("Tipo conexao  : RJDBC\n")
		cat("Driver conexao:", driverJDBC,"\n")
		cat("Nome banco    : SVUXPCAP2\n")
		cat("Nome usuario  :", userJDBC, "\n")
		cat("Nome senha    : *****\n")
		cat("Tentativa de conexao\n")

		connect <- c()
		connect <- Try_error(cmis_JDBC_Connect(driverJDBC, userJDBC, passJDBC))

		if(tolower(class(connect)[1]) == "jdbcconnection") {
		cat("Conexao estabelecida com sucesso!\n")
		} else stop(cat("Falha de conexao! Verifique os parametros de entrada!\n"))

		#cat("Lancando query SQL para coleta de dados!\n")
		Query <- 	paste("SELECT to_char(ger.DATA,'dd/mm/yyyy hh24:mi:ss') as data,
		  dat.realizado,
		  ger.id_metrica_frequencia
		FROM (
		  SELECT *
		  FROM (
			SELECT id_metrica_frequencia
			FROM cmis_owner.dsh_metrica_frequencia m
			WHERE id_metrica_frequencia BETWEEN ", idmini, " AND ", idmfin, " AND m.frequencia = ", Freq, "
			GROUP BY id_metrica_frequencia
		  ) md, (
			SELECT CASE ", Freq, "
			  WHEN '3' THEN trunc(SYSDATE-1) - LEVEL+1
			  WHEN '2' THEN add_months(trunc(add_months(SYSDATE, -1), 'MM'), - LEVEL+1)
			  WHEN '1' THEN trunc(SYSDATE) - ((1/24) * LEVEL)
			END AS DATA
			FROM dual
			CONNECT BY LEVEL <=
			  CASE ", Freq, "
				WHEN '3' THEN 90
				WHEN '2' THEN 24
				WHEN '1' THEN 120
			  END
		  )
		) ger LEFT JOIN (
		  SELECT DATA,
			realizado,
			metrica_freq,
			nova_observacao
		  FROM cmis_owner.dsh_dados_historicos h
		  INNER JOIN cmis_owner.dsh_metrica_frequencia m ON h.metrica_freq = m.id_metrica_frequencia
		  WHERE metrica_freq BETWEEN ", idmini, " AND ", idmfin, " AND DATA < SYSDATE AND m.frequencia = ", Freq, "
		) dat ON dat.metrica_freq = ger.id_metrica_frequencia AND dat.DATA = ger.data where exists ( select null
								from cmis_owner.dsh_dados_historicos tmp   
								where tmp.metrica_freq = dat.metrica_freq
								and nova_observacao = 1
								)		
		ORDER BY ger.id_metrica_frequencia, to_date(ger.DATA,'dd/mm/yyyy hh24:mi:ss')", sep="")

		cat("Coleta iniciada...\n")
		dados <- Try_error(DBI::dbGetQuery(connect, Query))

		if(class(dados) != "try-error") {
			cat("")
		} else {
			stop(cat("Coleta mal sucedida, verifique os parametros de coleta!\n"))
		}

		## Split dos dados por nome da métrica e repassagen para lista.
		dados <- split(dados, dados$ID_METRICA_FREQUENCIA, drop = TRUE)

		cat("Coleta bem sucedida...\nProcessando dados!\n\n")
		D <- Try_error(llply(dados, function(X) {
			if (all(is.na(X$REALIZADO))) NULL else X}
			, .progress =  "time"))
		if(class(D) !="try-error") {
			cat("\nDados processados com sucesso!\n")
		} else {
			stop(cat("Coleta mal sucedida, verifique os parametros de coleta!\n"))
		}

		Nm <- table(sapply(D, class))
		Ndf <- Nm[names(Nm)=="data.frame"]
		Nnu <- Nm[!names(Nm)=="data.frame"]

		cat("Total de metricas no range  :", sum(Nm), "\n")
		cat("Total de metricas COM dados :", Ndf, "\n")
		cat("Total de metricas SEM dados :", Nnu, "\n")

		#dados <- D[-(which(sapply(D, is.null), arr.ind=TRUE))]
		dados <- Filter(Negate(function(X) is.null(unlist(X))), D)
		
		DBI::dbDisconnect(connect)
		cat("Conexao encerrada!\n")
		cat("-------------------------------------------------------------------\n\n")
	}
	return(invisible(dados))
}

#' Test if an object exists
#' @export
testObject <- function(object){
  exists(as.character(substitute(object)))
}

#' Default summary function
#' @export
tsSummary <- function(P,A) {
data.frame((as.data.frame(accuracy(P,A))))
}


#' Mean forecast wrapper
#' @export
meanForecast <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::meanf(x, h, ..., level = level)
  if (onlyfc) fc$mean
  else fc
}

#' Naive forecast wrapper
#' @export
naiveForecast <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::naive(x, h, ..., level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Seasonal naive forecast wrapper
#' @export
snaiveForecast <- function(x, h, level = 95, onlyfc=TRUE, ...) {
  fc <- forecast::snaive(x, h, ..., level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Random walk forecast wrapper
#' @export
rwForecast <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fc <- forecast::rwf(x, h, ..., level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Theta forecast wrapper
#' @export
thetaForecast <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fc <- forecast::thetaf(x, h, ..., level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Linear model forecast wrapper
#' @export
lmForecast <- function(x,h,level=95, onlyfc=TRUE, xreg=NULL,newxreg=NULL,...) {
  x <- data.frame(x)
  colnames(x) <- 'x'
  if (is.null(xreg) & is.null(newxreg)) {
    fit <- tslm(x ~ trend + season, data=x, ...)
	fc <- forecast(fit, h=h, level=level)
	if (onlyfc) {
		fc$mean
	}  else {
		fc
	}

  } else if ((!is.null(xreg)) & !(is.null(newxreg))) {
    newnames <- c('x',colnames(xreg))
    x <- cbind(x,xreg)
    colnames(x) <- newnames
    fmla <- as.formula(paste("x ~ trend + season +", paste(colnames(xreg), collapse= "+")))
    fit <- tslm(fmla, data=x, ...)

	fc <- forecast(fit, h=h, level=level, newdata=newxreg)
	fc_mean <- fc$mean
	if (onlyfc) {
		fc$mean
	} else {
		fc
	}
  } else {
    stop('xreg and newxreg must both be NULL or both be provided')
  }
}

#' Structural time series forecast wrapper
#' @export
stsForecast <- function(x,h,level=95,  onlyfc=TRUE, ...) {
  fit <- StructTS(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Stl forecast wrapper
#' @export
stl.Forecast <- function(x, h, level=95, method='ets',  onlyfc=TRUE, ...) {
  fc <- forecast::stlf(x, h=h, method=method, level=level, ...)
  if (onlyfc) fc$mean
  else fc
}

#' Arima forecast wrapper
#' @export
arimaForecast <- function(x,h,level=95, onlyfc=TRUE, xreg=NULL,newxreg=NULL,...) {
  fit <- forecast::Arima(x, xreg=xreg, ...)
  fc <- forecast::forecast(fit, h=h, level=level, xreg=newxreg)
  if (onlyfc) fc$mean
  else fc
}

#' auto.arima forecast wrapper
#' @export
auto.arimaForecast <- function(x,h,level=95, onlyfc=TRUE, xreg=NULL,newxreg=NULL,...) {
  fit <- forecast::auto.arima(x, xreg=xreg, ...)
  fc <- forecast::forecast(fit, h=h, level=level, xreg=newxreg)
  if (onlyfc) fc$mean
  else fc
}

#' Ets forecast wrapper
#' @export
etsForecast <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fit <- forecast::ets(x, allow.multiplicative.trend = FALSE, additive.only = TRUE, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' BATS forecast wrapper
#' @export
batsForecast <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fit <- forecast::bats(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' TBATS forecast wrapper
#' @export
tbatsForecast <- function(x,h,level=95, onlyfc=TRUE, ...) {
  fit <- forecast::tbats(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' NNetar forecast wrapper
#' @export
nnetarForecast <- function(x, h, level=95,  onlyfc=TRUE, nn_p=1, ...) {
  fit <- forecast::nnetar(x, p=nn_p, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}


#' ses forecast wrapper
#' @export
sesForecast <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- forecast::ses(x, h=h, level=level, ...)
  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Holt forecast wrapper
#' @export
holtForecast <- function(x, h, level=95,  onlyfc=TRUE, ...) {
  fit <- forecast::holt(x, h=h, level=level, ...)
  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters forecast wrapper
#' @export
hwForecast <- function(x, h, level=95,  onlyfc=TRUE, ...) {
	## Check Seasonality, if YES, use Holt, else use HW
	if (sum(cycle(x))==length(x)) {
		fit <- forecast::holt(x, ...)
	} else {
		fit <- forecast::hw(x, ...)
	}

  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' Meanf forecast wrapper
#' @export
meanForecast <- function(x, h, level=95,  onlyfc=TRUE, ...) {
  fit <- forecast::meanf(x, ...)
  fc  <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}


#' HoltWinters Sazonal forecast wrapper by stats
#' @export
HWsForecast <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters Non Seazonal forecast wrapper
#' @export
HWnsForecast <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, gamma = FALSE, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else fc
}

#' HoltWinters Exponential Smoothing forecast wrapper
#' @export
HWesForecast <- function(x, h, level=95, onlyfc=TRUE, ...) {
  fit <- HoltWinters(x, gamma = FALSE, beta = FALSE, ...)
  fc <- forecast::forecast(fit, h=h, level=level)
  if (onlyfc) fc$mean
  else (fc)
}


switch.cvforecast <- function(x, nmodelo, h, level=95, onlyfc=FALSE) {
  switch(nmodelo,
		stsForecast = stsForecast(x, h=h, level=level, onlyfc=onlyfc),
		hwForecast = hwForecast(x, h=h, level=level, onlyfc=onlyfc),
		tbatsForecast = tbatsForecast(x, h=h, level=level),
		auto.arimaForecast = auto.arimaForecast(x, h=h, level=level, onlyfc=onlyfc),
		#stl.Forecast = stl.Forecast(x, h=h, level=level),
		sesForecast = sesForecast(x, h=h, level=level, onlyfc=onlyfc),
		meanForecast = meanForecast(x, h=h, level=level, onlyfc=onlyfc),
		holtForecast = holtForecast(x, h=h, level=level, onlyfc=onlyfc),
		batsForecast = batsForecast(x, h=h, level=level, onlyfc=onlyfc),
		etsForecast = etsForecast(x, h=h, level=level, onlyfc=onlyfc),
		arimaForecast = arimaForecast(x, h=h, level=level, onlyfc=onlyfc),
		lmForecast = lmForecast(x, h=h, level=level, onlyfc=onlyfc),
		thetaForecast = thetaForecast(x, h=h, level=level, onlyfc=onlyfc),
		rwForecast = rwForecast(x, h=h, level=level, onlyfc=onlyfc),
		snaiveForecast = snaiveForecast(x, h=h, level=level, onlyfc=onlyfc),
		naiveForecast = naiveForecast(x, h=h, level=level, onlyfc=onlyfc),
		meanForecast = meanForecast(x, h=h, level=level, onlyfc=onlyfc),
		nnetarForecast = nnetarForecast(x, h=h, level=level, onlyfc=onlyfc),
		HWsForecast = HWsForecast(x, h=h, level=level, onlyfc=onlyfc),
		HWnsForecast = HWnsForecast(x, h=h, level=level, onlyfc=onlyfc),
		HWesForecast = HWesForecast(x, h=h, level=level, onlyfc=onlyfc)
	)
}


## Use Mann-Kendall test (MK) and the Seasonal and the Regional Kendall Tests for trend (SKT and RKT) and Theil-Sen's slope estimator for checking trend
nparTrend <- function (x, npoints = 30, ...) {
#cat("Check trend in the last", npoints, "data points!\n")
if(length(x) > 12 & npoints < length(x)) x <- x[(length(x)-(npoints-1)):length(x)] 

    rk <- rkt::rkt(time(x), x)
    if (round(rk$B, 6) != 0) {
        if (rk$tau > 0.1) {
            if (rk$tau <= 0.95) 
                sinal <- 1L
            else sinal <- 2L
        }
        else if (rk$tau < -0.1) {
            if (rk$tau > -0.95) 
                sinal <- -1L
            else sinal <- -2L
        }
        else {
            sinal <- 0L
        }
    }
    else {
        sinal <- 0L
    }
    out <- c(slope = rk$B, tau = rk$tau, score = rk$S, p.value = rk$sl, 
        trend_sign = sinal)
    round(out, 4)
}


## Error handler improved
Try_error <- function(code, silent = TRUE) {
	W <- NULL
	w.handler <- function(w){
		W <<- w
		invokeRestart("muffleWarning")
	}
	withCallingHandlers(tryCatch(code, error = function(c) {
    msg <- conditionMessage(c)
    if (!silent) message(c)
    invisible(structure(msg, class = "try-error"))
  }), warning = w.handler)
}



linearityTest <- function(x, Test) {
	if (class(x)[1] != "ts") x <- ts(x)
	if (missing(Test)) Test <- "keenan"
	else {
		Test <- match.arg(Test, c("terasvirta","white", "keenan", "mcleodLi", "tsay","tarTest"))
	}
	
	test <- switch(Test,
			terasvirta = tseries::terasvirta.test(x = x, type = "Chisq"),
			white = tseries::white.test(x),
			keenan = TSA::Keenan.test(x),
			mcleodLi = TSA::McLeod.Li.test(y = x, plot = FALSE),
			tsay = TSA::Tsay.test(x),
			tarTest = TSA::tlrt(x)
	)

	if (Test == "terasvirta") {
		out <- data.frame(statistic = test$statistic, p.value = test$p.value)
	} else if (Test == "white") {
		out <- data.frame(statistic = test$statistic, p.value = test$p.value)
	} else if (Test == "keenan") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else if (Test == "tsay") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else if (Test == "tarTest") {
		out <- data.frame(statistic = test$test.stat, p.value = test$p.value)
	} else if (Test == "mcleodLi") {
		out <- data.frame(statistic = NA, p.value = max(unlist(test$p.values)))
	} else {
		cat("Nenhum dos testes se aplica!\n\n")
	}

	out <- round(out, 4)
	rownames(out) <- Test
	return(out)
}


ConvertDataToTs <- function(Data, tsfrequency="month", OutType = "ts", OutlierClean = TRUE, ...) {
#require(timeSeries)
require(lubridate)
  #check dataset and type
  OutType <- match.arg(OutType, c("ts","xts"))

  #check data frequancies
  tsfrequency <- match.arg(tsfrequency, c("year","month","day","hour","min","sec"))

  if(is.null(Data)) stop("Empty dataset!\n")
  
  if(class(Data)[1]  %in% c("ts","mts")) {
	cat("Data already 'time series'!\n")
	return(Data)
	}

  # check ts frequency
  if (tsfrequency %in% c("year","month")){
    freq <- 12
  } else if (tsfrequency=="day") {
    freq <- 7
  } else if (tsfrequency=="hour") {
    freq <- 24
  } else if (tsfrequency=="minute"){
    freq <- 60
  } else {
    freq <- 1
  }

  # check if data is a simple vector
  if(ncol(as.data.frame(Data)) == 1 & class(Data)[1] == "numeric") {
	DateSim <- seq(as.POSIXct(Sys.Date()), by = tsfrequency, l=length(Data))
	Datats <- ts(Data, frequency=freq)
    if (OutType == "ts") {
		return(Datats)
	} else {
		return(xts::xts(Datats, DateSim, frequency=freq))
	}
  }

  ## dates
  date  <- Data[, 1]
  value <- as.matrix(Data[,-1, drop=FALSE])

  #check as.character date pattern
  if (class(date) %in% c("factor","character")) {
    if (any(grepl("/",date))) {
		tdata <- TimeSeries(date, "%d/%m/%Y %H:%M:%S", value)
		date  <- tdata[, 1]
    } else if (any(grepl("-",date))) {
		tdata <- TimeSeries(date, "%d-%m-%Y %H:%M:%S", value)
		date  <- tdata[, 1]
    } else (stop("Incorrect date format!\n"))
  } else if (inherits(time(Data), "Date") || inherits(time(Data), "POSIXt")){
  		tdata <- zoo2TimeSeries(Data)
		date  <- tdata[, 1]
  } else (stop("Incorrect date format!\n"))

  outXTS <- xts::xts(value, date, frequency = freq)
  outTS  <- ts(value, start = Start(tsfrequency, tdata), frequency=freq)

  ## Remove outliers from the data after transformation
  if(OutlierClean) {
	value <- sapply(outTS, function(X) {
		## first of last missing? Insert avg + sd
		if(is.na(X[1])) X[1] <- c(mean(X, na.rm=TRUE) + sd(X, na.rm=TRUE))
		if(is.na(X[length(X)])) X[length(X)] <- c(mean(X, na.rm=TRUE) + sd(X, na.rm=TRUE))
		## Otherwise uses adequate method of treating missing and outliers
		temp <- try(tsclean(X))
		if (class(temp)!="try-error") temp else X
		}
	)
	outXTS <- xts::xts(value, date, frequency = freq)
    outTS  <- ts(value, start = Start(tsfrequency, tdata), frequency=freq)
	#outZOO <- zoo::as.zoo(outXTS, frequency=freq)
	}
	if (OutType == "ts") {
		return(outTS)
	} else if (OutType == "xts"){
		return(outXTS)
	} else {
		stop("Check data!\n")
	}
}


ForecastHorizon <- function(XtsData, tsfrequency, horizon) {
	if(!class(XtsData)[1] %in% c("zoo","xts","ts","numeric")) stop("Data must be of 'zoo', 'xts', 'ts' or 'numeric' classes!\n")

	#check date arg
	tsfrequency <- match.arg(tsfrequency, c("year","month","day","hour","min","sec"))

	# check ts frequency
	if (tsfrequency %in% c("year","month")){freq <- 12
	} else if (tsfrequency=="day") {freq <- 7
	} else if (tsfrequency=="hour") {freq <- 24
	} else {freq <- 1}

	if (class(XtsData)[1] %in% c("ts","numeric")) {
		SE <- seq(as.POSIXct(Sys.Date(), tz = "UTC"), by = tsfrequency, l=horizon)

		if (class(XtsData)[1] == "numeric") {
			x <- ts(XtsData, frequency=freq)
		} else x <- XtsData

		return(list(x, FCHorizon = SE))
	} else {
		IndexXTS <- index(XtsData)
		IndexXTS <- IndexXTS[!is.na(IndexXTS)] # Remove missing cado haja
		SE <- seq(from=max(IndexXTS), by=tsfrequency, length.out = horizon+1)[-1]
		value  <- as.numeric(XtsData)

		tdata <- TimeSeries(as.character(IndexXTS), "%Y-%m-%d", value)

		if(class(XtsData)[1] == "xts") {
			x <- ts(value, start = Start(tsfrequency, tdata), frequency = freq)
		} else {
			x <- XtsData
		}
	}
	return(list(x, FCHorizon = SE))
}


find.freq <- function(x) {
    n <- length(x)
    spec <- spec.ar(c(x),plot=FALSE)
    if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
    {
        period <- round(1/spec$freq[which.max(spec$spec)])
        if(period==Inf) # Find next local maximum
        {
            j <- which(diff(spec$spec)>0)
            if(length(j)>0)
            {
                nextmax <- j[1] + which.max(spec$spec[j[1]:500])
                period <- round(1/spec$freq[nextmax])
            }
            else
                period <- 1
        }
    }
    else
        period <- 1
    return(period)
}

Start <- function(freq, tdata) {
  if (!freq %in% c("year","month","day","hour","min")) {
	return(c(1))
  } else {
	  switch(freq,
			year  = c(tdata$year[1], tdata$month[1]),
			month = c(tdata$year[1], tdata$month[1]),
			day   = c(tdata$week[1], tdata$day[1]),
			hour  = c(tdata$day[1],  tdata$hour[1]),
			min   = c(tdata$hour[1], tdata$minute[1])
			)
	   }
}

descritiva <-
function (x, basic = TRUE, desc = TRUE, norm = FALSE, p = 0.95, dig = 6) {
options(digits = dig)
    stat.desc.vec <- function(x, basic, desc, norm, p) {
        x <- unlist(x)
        if (!is.numeric(x)) {
            Nbrval <- NA
            Nbrnull <- NA
            Nbrna <- NA
            Median <- NA
            Mean <- NA
            StdDev <- NA
            if (basic == TRUE) {
                Res1 <- list(nbr.val = NA, nbr.null = NA, nbr.na = NA,
                  min = NA, max = NA, range = NA, sum = NA)
            }
            else Res1 <- NULL
            if (desc == TRUE) {
                CIMean <- NA
                names(CIMean) <- p
                Res2 <- list(median = NA, mean = NA, SE.mean = NA,
                  CI.mean = NA, var = NA, std.dev = NA, coef.var = NA)
            }
            else Res2 <- NULL
            if (norm == TRUE) {
                Res3 <- list(skewness = NA, skew.2SE = NA, kurtosis = NA,
                  kurt.2SE = NA, normtest.W = NA, normtest.p = NA)
            }
            else Res3 <- NULL
        }
        else {
            Nbrna <- sum(as.numeric(is.na(x)))
            x <- x[!is.na(x)]
            Nbrval <- length(x)
            Nbrnull <- sum(as.numeric(x == 0))
            if (basic == TRUE) {
                Min <- min(x)
                Max <- max(x)
                Range <- Max - Min
                Sum <- sum(x)
                Res1 <- list(nbr.val = Nbrval, nbr.null = Nbrnull,
                  nbr.na = Nbrna, min = Min, max = Max, range = Range,
                  sum = Sum)
            }
            else Res1 <- NULL
            Median <- median(x)
            names(Median) <- NULL
            Mean <- mean(x)
            Var <- var(x)
            StdDev <- sqrt(Var)
            SEMean <- StdDev/sqrt(Nbrval)
            if (desc == TRUE) {
                CIMean <- qt((0.5 + p/2), (Nbrval - 1)) * SEMean
                names(CIMean) <- p
                CoefVar <- StdDev/Mean
                Res2 <- list(median = Median, mean = Mean, SE.mean = SEMean,
                  CI.mean = CIMean, var = Var, std.dev = StdDev,
                  coef.var = CoefVar)
            }
            else Res2 <- NULL
            if (norm == TRUE) {
                Skew <- sum((x - mean(x))^3)/(length(x) * sqrt(var(x))^3)
                Kurt <- sum((x - mean(x))^4)/(length(x) * var(x)^2) -
                  3
                SE <- sqrt(6 * Nbrval * (Nbrval - 1)/(Nbrval -
                  2)/(Nbrval + 1)/(Nbrval + 3))
                Skew.2SE <- Skew/(2 * SE)
                SE <- sqrt(24 * Nbrval * ((Nbrval - 1)^2)/(Nbrval -
                  3)/(Nbrval - 2)/(Nbrval + 3)/(Nbrval + 5))
                Kurt.2SE <- Kurt/(2 * SE)
                Ntest <- shapiro.test(x)
                Ntest.W <- Ntest$statistic
                names(Ntest.W) <- NULL
                Ntest.p <- Ntest$p.value
                Res3 <- list(skewness = Skew, skew.2SE = Skew.2SE,
                  kurtosis = Kurt, kurt.2SE = Kurt.2SE, normtest.W = Ntest.W,
                  normtest.p = Ntest.p)
            }
            else Res3 <- NULL
        }
        Res <- unlist(c(Res1, Res2, Res3))
        if (length(Res) == 0)
            Res <- unlist(list(nbr.val = Nbrval, nbr.null = Nbrnull,
                nbr.na = Nbrna, median = Median, mean = Mean,
                std.dev = StdDev))
        Res
    }
    Basic <- basic
    Desc <- desc
    Norm <- norm
    P <- p
    if (is.vector(x))
        stat.desc.vec(x, Basic, Desc, Norm, P)
    else {
        x <- as.data.frame(x)
        NamesV <- names(x)
        StatM <- NULL
        for (i in 1:ncol(x)) {
            StatV <- stat.desc.vec(x[i], Basic, Desc, Norm, P)
            if (is.null(StatM) == TRUE)
                StatM <- data.frame(StatV)
            else StatM <- cbind(StatM, StatV)
        }
        names(StatM) <- NamesV
        t(StatM)
    }
}

## Mudança: Adicionado LSD, SMAPE e LAR como estatisticas de bondade e mudada a decisão para accuracia ao invés de RANK SOMA
##x <- ConvertDataToTs(dd[[10]][,1:2], tsfrequency = "day",OutType = "ts")
##MFc <- try(MultiForecast(x, fcMethod=fcMethod, Control=Control))

BestModel <- function(objMforecast, Control=cmisControl(), ...) {

	residlevel <- Control$residlevel
	cvMethod <- Control$cvMethod

	if(class(objMforecast) != "multiforecast") stop("Objeto deve ser de classe 'multiforecast'")

	## Análise de residuos dos modelos
    Resid <- try(plyr::ldply(objMforecast, function(X) {
    if (class(X)[1] != "try-error") {
      sm <- sum(Mresid(X) > residlevel)
	  sm[sm %in% c(Inf, -Inf, NA, NULL)] <- 99
	  sm
    } else NULL # Trata casos de erros com -9
	}))
	
	if (class(Resid) != "try-error") {
		# Remove entradas missing (se houver)
		Resid <- Filter(Negate(function(X) is.null(unlist(X))), Resid)
		names(Resid) <- c(".id","RESID.PVALUE")
	  
		# Cria rank de resíduo
		Resid$RES.RANK <- rank(Resid$RESID.PVALUE, na.last = TRUE, ties.method = "first")
		Resid <- Resid[,-2]
	} else {
		Resid <- data.frame(.id = NA, RES.RANK = NA)
	}
	
	## Estatisticas de bondade
	STATS <- Try_error(plyr::ldply(objMforecast, function(X) {
		if (class(X)[1] != "try-error") {
			st <- round(Accuracy(X), 4)
			st[st %in% c(Inf, -Inf, NA, NULL)] <- 999999999
			st
		} else 999999999 # Trata casos de erros com 999.999.999
	}))	
	
	if(class(STATS) != "try-error") {
		STATS <- Filter(Negate(function(X) is.null(unlist(X))), STATS)
	    ## Cria rank de estatística de bondade
		STATS$RESID.GOF <- rank(STATS[,cvMethod], na.last = TRUE, ties.method = "first")
	} else {
		STATS <- data.frame(.id=NA, ME=NA, RMSE=NA, MAE=NA, MPE=NA, MAPE=NA, MASE=NA, SMAPE=NA, LSD = NA, LAR = NA, RESID.GOF=NA)
	} 
  
	if (!all(is.na(STATS)) & !all(is.na(Resid))) {  
		ESTFINAL <- merge(STATS, Resid, by.x = ".id")
		ESTFINAL$SOMA.RK <- rowSums(ESTFINAL[,c("RESID.GOF", "RES.RANK")], na.rm = FALSE)	
		ESTFINAL <- ESTFINAL[order(ESTFINAL[,cvMethod]),]	
		CV_names <- ESTFINAL[,".id"]
	} else {
		CV_names <- NA
	}
	
	return(list(as.character(CV_names), GOF=ESTFINAL))
}

#' Estatisticas de acuracia
#'
#' Calcula as estatísticas "ME", "RMSE", "MAE", "MPE", "MAPE", "SMAPE", "LSD" e "LAR"
#'
#' @param fit objeto da classe forecast
#' @param dig número de digitos da saida
#' 
#' @export
Accuracy <- function(fit, dig = 4) {
  atuals    <- getResponse(fit)
  forecasts <- fitted(fit)
  residuos  <- atuals-forecasts  
  perc.err  <- 100*(residuos/forecasts)
  perc.err  <- na.omit(perc.err[!is.infinite(perc.err)])  
  N <- length(atuals)  
  Q <- abs(as.numeric(forecasts/atuals))
  Q <- na.omit(Q[!is.infinite(Q)])
  Q1 <- log(Q)^2
  L <- (0.5*var(Q, na.rm = TRUE, use = "complete.obs") - log(Q))^2
  L <- na.omit(L[!is.infinite(L)])
  
  me  <- mean(residuos, na.rm = TRUE)
  mse <- mean(residuos^2, na.rm = TRUE)
  mae <- mean(abs(residuos), na.rm = TRUE)
  mpe <- mean(perc.err, na.rm = TRUE)
  mape  <- mean(abs(perc.err), na.rm = TRUE)
  
  smape <- mean((abs(residuos)/((abs(atuals)+abs(forecasts))/2)), na.rm = TRUE)
  lsd <- sqrt(sum(L, na.rm=TRUE)/(length(L)-1))
  lar <- sum(na.omit(Q1[!is.infinite(Q1)]), na.rm=TRUE)  
  
  out <- c(me, sqrt(mse), mae, mpe, mape, 100*smape, 100*lsd, lar)
  names(out) <-  c("ME", "RMSE", "MAE", "MPE", "MAPE", "SMAPE", "LSD", "LAR")
  
  return(round(out, dig))
}


