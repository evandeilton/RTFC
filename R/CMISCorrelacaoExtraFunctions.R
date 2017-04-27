################ FUN??ES AUXILIARES CMIS FORECAST ################

# Require packages
#Rpacks()

## Testa os dados antes de iniciar os c?lculos
dacheck <- function(dados) {
	# Confere numero de linhas e colunas

	if (is.null(dim(dados)) | sum(dim(dados)) < 1 | class(dados) != "data.frame" | sum(nrow(dados)) < 6) {
		cat("ERRO: Verifique a tabela de dados, deve ser data.frame e possuir variaveis explicativas e explicadas,\n pelo menos 6 linhas e 3 colunas inluindo dada!\n")

		cat("ERRO: Os dados informados possuem:", sum(nrow(dados)), "linhas", sum(ncol(dados)), "colunas", "e classe", class(dados), "!\n")
		stop()
	}
	# Verifica variaveis rexplicativas e explicadas
	tb_prefixo <- c("bco","apl","wbs","tbs","hdw","neg")
	nm <- tolower(names(dados))
	if (substr(nm[1],1,3) == "dat") {
		sb <- table(substr(nm[-1], 1, 3))
	} else {
		sb <- table(substr(nm, 1, 3))
	}
	#vars <- tb_prefixo[which(tb_prefixo == names(sb))]
	#vars <- intersect(names(sb),tb_prefixo)
	vars <- tb_prefixo[match(names(sb),tb_prefixo)]

	if (dim(sb) < 2 | length(vars) < 2 | any(is.na(vars))) {
		stop(cat("ERRO: Verifique o conjunto de dados. Apenas os prefixos ", names(sb), " foram encontrados!\n Nao foram encontradas variaveis resposta e/ou explicativas! \n\n"))
	}
	if (dim(sb) > 1 | length(vars) > 1) {
		cat("OK: Encontrados", sb[vars[1]], "tipo", vars[1], "e", sb[vars[2]], "tipo", vars[2], "!\n")
		return(dados)
	}
}


## Estat?stica de for?a da predi??o do modelo (PRESS - PREdicted Residual Sum of Squares)
## Fonte: https://github.com/WinVector/Examples/blob/master/PRESS/fastPRESS.R
## 		  http://www.r-bloggers.com/estimating-generalization-error-with-the-press-statistic/
my_lm_press <- function(obj, wts=c()) {
	## Amostras de valida??o para Y
	hold1outmeans <- function(y, wts=c()) {
	  # get per-datum hold-1 out grand means
	  if(is.null(wts)) {
		wts <- rep(1.0,length(y))
	  }
	  sumY <- sum(y*wts)
	  sumW <- sum(wts)
	  meanP <- (sumY - y*wts)/(sumW - wts)
	  meanP[is.na(meanP)] <- 0.5
	  meanP
	}

	## Amostras de valida??o para X
	hold1outlmpreds <- function(formula, data, wts=c()) {
	  formula <- as.formula(formula)
	  nRows <- dim(data)[[1]]
	  if(is.null(wts)) {
		wts <- rep(1.0,nRows)
	  }
	  x <- model.matrix(formula,data=data)
	  nVars <- dim(x)[[2]]
	  terms <- terms(formula)
	  yvarName <- as.character(as.list((attr(terms,'variables')))[[attr(terms,'response')+1]])
	  y <- data[,yvarName]
	  xTx <- t(x) %*% ( wts * x ) + 1.0e-5*diag(nVars)
	  xTy <- t(x) %*% ( wts * y )
	  pi <- function(i) {
		xTxi <- xTx - wts[i] * (x[i,] %o% x[i,])
		xTyi <- xTy - wts[i] * (x[i,] * y[i])
		betai <- solve(xTxi,xTyi)
		predi <- as.numeric(x[i,] %*% betai)
	  }
	  preds <- sapply(1:nRows,pi)
	  preds
	}

	# Chamada da fun??o

	data <- obj$model
	formul <- as.formula(formula(obj))

	ind <- as.character(formul[[2]])
	y <- data[, ind]
	n = length(y)
	hopreds = hold1outlmpreds(formul, data)
	homeans = hold1outmeans(y)
	devs = y-hopreds
	PRESS = sum(devs^2)
	rmPRESS = sqrt(mean(devs^2))
	dely = y-homeans
	PRESS.r2= abs(1 - (PRESS/sum(dely^2)))
	return(round(data.frame(PRESS=PRESS, rmPRESS=rmPRESS, PRESS.r.squared=PRESS.r2), 5))
}

#' Estat?sticas de qualidade e predi??o de modelos lm
#' Retorna as estat?sticas R-quadrado, R-quadrado parcial e PRESS de modelos lineares
#' @return Data frame com as estat?sticas do modelo
#' @param objeto das classe \code{lm()} model.
lm_press_stat <- function(obj) {
	if(class(obj)!="lm") stop("Apenas modelos lm sao validos!\n")
	#' @title PRESS
	#' @author Thomas Hopper
	#' @description Returns the PRESS statistic (predictive residual sum of squares).
	#'              Useful for evaluating predictive power of regression models.
	#' @param obj A linear regression model (class 'lm'). Required.

	PRESS <- function(obj) {
	  #' calculate the predictive residuals
	  pr <- residuals(obj)/(1-lm.influence(obj)$hat)
	  #' calculate the PRESS
	  PRESS <- sum(pr^2)
	  return(PRESS)
	}

	#' @title Predictive R-squared
	#' @author Thomas Hopper
	#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
	#'              the PRESS statistic.
	#' @param obj A linear regression model (class 'lm'). Required.

	pred_r_squared <- function(obj) {
	  #' Use anova() to get the sum of squares for the linear model
	  lm.anova <- anova(obj)
	  #' Calculate the total sum of squares
	  tss <- sum(lm.anova$'Sum Sq')
	  # Calculate the predictive R^2
	  pred.r.squared <- 1-PRESS(obj)/(tss)

	  return(pred.r.squared)
	}

	r.sqr <- summary(obj)$r.squared
	adj.r.sqr <- summary(obj)$adj.r.squared
	pre.r.sqr <- pred_r_squared(obj)
	PRESS <- PRESS(obj)

	return.df <- data.frame(press.r.squared = pre.r.sqr, press = PRESS, r.squared = r.sqr, adj.r.squared = adj.r.sqr)
	return(return.df)
}


identico <- function(x, interp=FALSE) {
	if (interp) {
	#Caso haja valores missing: Interpola??o linear!
		if (any(is.na(x))) {

			# M?dia nos extremos
			if (is.na(x[1])) x[1] <- mean(x, na.rm=TRUE)
			if (is.na(x[length(x)])) x[length(x)] <- mean(x, na.rm=TRUE)

			# Interpola??o linear nos outros pontos.
			x <- na.interp(x)
		}
	}
	# Comprimento dos dados
	ntot <- length(x)
	ndup <- length(x[duplicated(x)])
	sdto <- sd(x, na.rm=TRUE)

	## Remover m?tricas com mais de 90% repetidos ou com desvio padr?o muito baixo!
	if(ndup/ntot > 0.90 | sdto < 0.00001 | is.logical(x)) out <- 1 else out <- 0
	# 0 = OK, 1 = Deleta
	out
}

# Tratar outliers em data.frames
dataframe.outlier <- function(dados, n_passos_frente=20, freq=7, normalize = FALSE, ...) {
options (digits = 10)
#Rpacks()
	bb <- aa <- temp_a <- out <- c()

	## Remove colunas com todos elementos missing
	dados_temp0 <- dados[sapply(dados, function(x) length(unique(na.omit(x)))) > 1];

	## Remove colunas com valores constantes. Desvio padr?o muito pr?ximo de ZERO

	dt1 <- dados_temp0[, -1, FALSE]
	condicao <- sapply(dt1, identico)
	dt1 <- dt1[,names(condicao[condicao==0])]

	dados_temp1 <- data.frame(dados_temp0[, 1, FALSE], dt1)

	## Arrumando o cabe?alho
	nomes <- names(dados_temp1)
	nomes[1] <- "data"
	names(dados_temp1) <- nomes
	datas <- dados_temp1[,1]

	for(i in 2:length(names(dados_temp1))){
	temp_a <- try(ts.dados(dados_temp1[,c(1, i)], n_passos_frente, freq, normalize = normalize)$historico_novo_ts)

		if (class(temp_a) != "try-error") {
			#temp_a <- temp_a[,1]
			assign(paste("a_", names(dados_temp1)[i], sep=''),  temp_a)
			temp_b <- try(tsclean(get(paste("a_", names(dados_temp1)[i], sep=""))))

			if (class(temp_b) != "try-error") {
				assign(paste(names(dados_temp1)[i], sep=''), temp_b)
				bb[i] <- paste(names(dados_temp1)[i], sep='')
			}
		}
	}
	aa <-as.character(na.omit(bb))
	dados_ts <- do.call(ts.union, mget(aa))

	if (normalize == TRUE) {
		dados_ts <- scale(dados_ts)
	}
	out$dados_ts <- dados_ts
	out$dados <- data.frame(datas, dados_ts)

	return(out)
}

# Recebe um objeto lm ou glm e calcula o VIF (Variance Infraction Factor)
vif<-function (obj, digits = 5) {

	objout <- function(obj, ...){
		if(length(coef(obj)) < 2) {
			ou <- NA
		} else {
			Qr <- obj$qr
			if (is.null(obj$terms) || is.null(Qr))
				stop("invalid 'lm' object:  no terms or qr component")
			tt <- terms(obj)
			hasintercept <- attr(tt, "intercept") > 0
			p <- Qr$rank
			if (hasintercept)
				p1 <- 2:p
			else p1 <- 1:p
			R <- Qr$qr[p1, p1, drop = FALSE]
			if (length(p1) > 1)
				R[row(R) > col(R)] <- 0
			Rinv <- qr.solve(R)
			vv <- apply(Rinv, 1, function(x) sum(x^2))
			ss <- apply(R, 2, function(x) sum(x^2))
			vif <- ss * vv
			ou <- signif(vif, digits)
		}
		return(ou)
	}
	temp <- try(objout(obj))
	if (class(temp) != "try-error") {
		out <- temp
	} else {
		out <- 9999
	}
	return(out)
}

# Faz an?lise de VIF (Variance Infraction Factor) nas m?tricas dos dados de modelagem, eliminando possibilidade de multicolinearidade que pode gerar estimativas falsas.
my_glm_vif <- function(dados, vresposta, explicativas, nvif=50, regsimples=FALSE, ...) {

	metricas <- c(vresposta, explicativas)
	ds_modelo <- dados[, metricas, drop = FALSE]

	## Utiliza OLS para escolher as vari?veis mais fortes com base no VIF
	mreg <- paste(vresposta, " ~ 1")
	fun <- as.formula(paste(c(mreg, explicativas), collapse=" + "))
	templm <- lm(fun, data = ds_modelo)
	if(class(templm) != "try-error") {
		modelocompleto <- templm

	## Dispensa da an?lise de VIF modelos lineares simples
		if (regsimples == TRUE){
			nvif <- Inf
		} else {
			nvif <- nvif
		}
	} else {
		modelocompleto <- "Nao foi possivel fazer a regressao inicial!"
	}

	# Remove as vari?veis com VIF maiores que um valor dado
	vifs <- vif(modelocompleto)
	nmvif <- names(vifs[vifs < nvif])
	nomesvif <- c(vresposta, nmvif)
	ds_refit <- ds_modelo[, nomesvif]

	funfinal <- as.formula(paste(c(mreg, nmvif), collapse=" + "))

	out <- c()
	out$dados <- ds_refit
	out$formula <- funfinal
	out$vif <- vifs[vifs < nvif]
	out$vresposta <- vresposta
	out$vexplicativas  <- nmvif
	return(out)
}

my_glm_cv <- function(dados, vresposta, explicativas, nfolds = 5, ...) {

	y <- dados[, vresposta]
	x <- as.matrix(dados[, explicativas])
	nfolds <- nfolds
	## Fixar semente aleat?ria para possibilidade de replica de codes
	set.seed(47743)

	cvfit <- try(cv.glmnet(x, y, alpha=1, type.measure = "deviance", nfolds = nfolds))

	if (class(cvfit) !="try-error"){
		co <- as.matrix(coef(cvfit, s = "lambda.1se"))
		co <- data.frame(coefnames = rownames(co), coef = as.numeric(co))
		co <- co[-1, ] # remove a primeira linha (do intercepto)
		rownames(co) <- NULL
		dd <- co[co$coef > 0, ]


		# Caso haja alguma vari?vel explicativa com coeficiente significativo (ou seja, que foi estimado via CrossValidation)
		if (length(dd$coef)>0) {
		#cat("Log: A variavel ", vresposta, "possui", length(dd$coef),"covariavel(eis)!\n")
			## Monta formula depois da crossValidation
			nmcv <- as.character(dd$coefnames)
			mreg <- paste(vresposta, " ~ 1")

			ds_refit <- dados[, c(vresposta, as.character(dd$coefnames))]

			funfinal <- as.formula(paste(c(mreg, nmcv), collapse=" + "))

			out <- c()
			out$dados <- ds_refit
			out$formula <- funfinal
			out$lambdacv <- cvfit$lambda.1se
			out$vresposta <- vresposta
			out$vexplicativas  <- nmcv
			return(out)
		}
	}
}



# Recebe um objeto lm ou glm e faz an?lise de residuos
analise_residuos_glm <- function(glm_obj) {
	out <- c()

		x <- glm_obj$model
		r <- as.numeric(residuals(glm_obj))

		di <- abs(nrow(x)-length(r))
		if (di > 0) {
			x <- x[-c(1:di),]
		}

	#### testes para independencia dos residuos
		independencia <- trycatch_w_e(Box.test(r, lag=10))$value
		# Teste para ver se a media tende a zero
		media_zero <- trycatch_w_e(t.test(r, alternative='two.sided', mu=0.0, conf.level=.95))$value

		# Teste para ver se os residuos sao ruido branco
		ruido_branco <- trycatch_w_e(LB.test(glm_obj, no.error=TRUE))$value

		# Teste para normalidade dos residuos jarque-bera
		normalidade <- trycatch_w_e(jarque.bera.test(na.omit(r)))$value

		# Teste de heterocedasticidade dos residuos p-valor >0,05 indica homocedasticidade
		homocedasticidade <- trycatch_w_e(bptest(r ~ x[,1]))$value

		# Teste de durbin-watson para autocorrelacao dos residuos se dw~2 ? independente
		autocorrelacao <- trycatch_w_e(dwtest(r ~ x[,1]))$value

		if (class(independencia$p.value) == "numeric")		{p0 <- as.numeric(independencia$p.value)} else {p0 <- NA}
		if (class(media_zero$p.value) == "numeric") 		{p1 <- as.numeric(media_zero$p.value)} else {p1 <- NA}
		if (class(ruido_branco$p.value) == "numeric") 		{p2 <- as.numeric(ruido_branco$p.value)} else {p2 <- NA}
		if (class(normalidade$p.value) == "numeric") 		{p3 <- as.numeric(normalidade$p.value)} else {p3 <- NA}
		if (class(homocedasticidade$p.value) == "numeric") 	{p4 <- as.numeric(homocedasticidade$p.value)} else {p4 <- NA}
		if (class(autocorrelacao$p.value) == "numeric") 	{p5 <- as.numeric(autocorrelacao$p.value)} else {p5 <- NA}


		df.peso <-c(ifelse(p0 > 0.05, 1, 0),
					ifelse(p1 > 0.05, 1, 0),
					ifelse(p2 > 0.05, 3, 0),
					ifelse(p3 > 0.05, 3, 0),
					ifelse(p4 > 0.05, 1, 0),
					ifelse(p5 > 0.05, 1, 0))

		df.pvalor <- c(p0, p1, p2, p3, p4, p5)

		analise_residual <- data.frame(df.peso, round(df.pvalor, 5))
		rownames(analise_residual) <- c("independencia","media_zero","ruido_branco","normalidade","homocedasticidade","autocorrelacao")

		# Indicador de qualidade da an?lise de residuos: menor = 0.10 = MELHOR, maior = PIOR
		out <- try(1/(colSums(analise_residual, na.rm=TRUE)[1]))
		if (class(out) != "try-erros") {
			out <- as.numeric(out)
		} else {
			out <- 0.00001
		}
		#out$analise_residual <- analise_residual
		#out$soma_peso <- colSums(analise_residual)[1]
		#out$soma_peso_inv <- 1/(colSums(analise_residual)[1])
		return(out)
}


# Tratar dados para entrada no algoritimo das regress?es
my_glm_data <- function(dados, n_passos_frente=20, freq=7, prefixo.expl, prefixo.resp, ... ) {
## Etapa 1: entrada e tratamento dos dados
	out <- dd <- c()

	#cat("\n========================================================================\n")
	#cat("01) Tramento dos dados! \n")
	#cat("========================================================================\n")

	if (!is.null(dados)) {
		names(dados) <- tolower(names(dados))

		# Selciona dados brutos, trata as datas e constr?i as s?ries de acordo com a frequencia
		# ds_temp  <- subset.dados(dados, aplicativo)

		# Faz a remo??o dos outliers
		ds_temp  <- dataframe.outlier(dados, n_passos_frente, freq)$dados
		#ds_modelo <- ds_temp[, c(vresposta, as.character(regressoras))]


		nomes <- names(ds_temp)
		regressoras <- respostas <- c()
		# Coloca nos nomes em objetos
		for (n in nomes) {
			if (substr(n, 1, 3) == prefixo.expl) {
				regressoras[n] <- as.character(n)
			}
			if (substr(n, 1, 3) == prefixo.resp) {
				respostas[n] <- as.character(n)
			}
		}

		dataapp <- nomes[1:2]
		respostas <- as.character(respostas)
		regressoras  <- as.character(regressoras)
		metricas_modelo <- c(dataapp, respostas, regressoras)

		if (length(respostas) < 1 | length(regressoras) < 1) {
			stop(cat("Para modelagem eh preciso pelo menos uma metrica explicativa e uma resposta!\n"),
			cat("Metricas resposta =", length(respostas), "\n"),
			cat("Metricas explicativas =", length(regressoras), "\n"))
		} else {
			ds_modelagem <- ds_temp[, c(respostas, regressoras)]
		}
	} else stop(cat("Conjunto de dados ", dados, " vazio!" ))


	out <- c()
	out$explicativas <- regressoras
	out$respostas    <- respostas
	out$ds_modelagem <- ds_modelagem
	return(out)
}


glm_extract_list <- function(fit, r2=0.50, residteste=0.50, ...) {
	## Coeficientes
	dfcoefnum   <-  as.data.frame(t(coef(summary(fit))[, "Estimate"]))

	## VIF Variance Inflaction Factor (Multicolinearidade das regressoras
	vifvalues   <- as.data.frame(t(vif(fit)))

	## N?mero de coeficientes
	#ncoef <- unlist(apply(regmat, 1, sum))

	## R^2 (R- quadrado, ou coeficiente de correla??o m?ltiplo do modelo
	R2       <- unlist(summary(fit)$r.squared)

	## R^2 (R- quadrado, ou coeficiente de correla??o m?ltiplo do modelo ajustado
	adjR2    <- unlist(summary(fit)$adj.r.squared)

	## Root Mean Squared Error (Raiz quadrada do erro m?dio quadr?tico)
	RMSE     <- unlist(summary(fit)$sigma)

	## AIC e BIC do modelo
	AIC     <- unlist(AIC(fit))

	## p-Valor do teste F dos modelos
	calcpval <- function(x){
		fstat <- summary(x)
		if(class(x)[1] == "glm") {
			fstat <- fstat$df
		} else {
			fstat <- fstat$fstatistic
		}
		pval <- trycatch_w_e(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE))$value

		if (class(pval) != "numeric")
		pval <- 99999
		return(pval)
	}

	fstats   <- unlist(calcpval(fit))

	## An?lise de residuos utilizando o indicador de qualidade do ajuste.
	residtest <- unlist(analise_residuos_glm(fit))

	## Estat?stica PRESS
	press <- unlist(my_lm_press(fit))

	formulas = deparse(formula(fit))

	coefs <- unlist(as.data.frame(t(coef(fit))))
	gridmod <- expand.grid(names(fit$model)[1], names(coefs))
	gridmod$coe    <- coefs
	gridmod$AIC <- round(AIC,5)
	gridmod$AIC <- round(AIC,5)
	gridmod$AIC = round(AIC,5)
	gridmod$R2 = round(R2, 5)
	gridmod$RMSE = round(RMSE,5)
	gridmod$residtest = round(residtest,5)
	gridmod$pF = round(fstats, 5)
	gridmod$PRESS = round(press[1], 5)
	gridmod$PRESSr2 =  round(press[3], 5)
	#gridmod$formula = formulas

	names(gridmod) <- c("v.resposta","explicativa", "coeficiente","aic","r2","rmse","resid.test","f.test", "press", "press.r2")

	## Seleciona apenas modelos com R2 > um valor definido na entrada
	out <- gridmod[(gridmod$r2 > r2 & gridmod$resid.test < residteste), ]
	out <- out[order(out$r2, decreasing = TRUE), ]
	return(out)
}

glm_extract_list_ci <- function(fit, r2=0.50, residteste=0.50, ...) {
	## Coeficientes
	dfcoefnum   <-  as.data.frame(t(coef(summary(fit))[, "Estimate"]))

	## VIF Variance Inflaction Factor (Multicolinearidade das regressoras
	vifvalues   <- as.data.frame(t(vif(fit)))

	## Intervalos de confian?a
	ci <-  as.data.frame(confint(fit))

	## p-valor
	pv <- round(summary(fit)$coefficients[, "Pr(>|t|)"], 5)

	## N?mero de coeficientes
	#ncoef <- unlist(apply(regmat, 1, sum))

	## R^2 (R- quadrado, ou coeficiente de correla??o m?ltiplo do modelo
	R2       <- unlist(summary(fit)$r.squared)

	## R^2 (R- quadrado, ou coeficiente de correla??o m?ltiplo do modelo ajustado
	adjR2    <- unlist(summary(fit)$adj.r.squared)

	## R^2 Parcial para cada regressora (Pacote relaimpo)



	# Se MRLM, extrair o R^2 parcial
	if (length(names(coef(fit))) > 2) {
		parc <- suppressWarnings(try(calc.relimp(fit)))
		if (class(parc) != "try-error") {
			r2parcial <- unlist(c(0, parc@lmg))
		} else {
			r2parcial <- rep(0, length(names(coef(fit))))
		}
	} else {
		# MRLS, utilizar o R^2 simples
		r2parcial <- unlist(c(0, R2))
	}

	## Root Mean Squared Error (Raiz quadrada do erro m?dio quadr?tico)
	RMSE     <- unlist(summary(fit)$sigma)

	## AIC e BIC do modelo
	AIC     <- unlist(AIC(fit))

	## p-Valor do teste F dos modelos
	calcpval <- function(x){
		fstat <- summary(x)
		if(class(x)[1] == "glm") {
			fstat <- fstat$df
		} else {
			fstat <- fstat$fstatistic
		}
		pval <- trycatch_w_e(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE))$value

		if (class(pval) != "numeric")
		pval <- 99999
		return(pval)
	}

	fstats   <- unlist(calcpval(fit))

	## An?lise de residuos utilizando o indicador de qualidade do ajuste.
	residtest <- unlist(analise_residuos_glm(fit))

	## Estat?stica PRESS
	## press <- unlist(my_lm_press(fit))
	press <- unlist(lm_press_stat(fit))

	formulas = deparse(formula(fit))

	coefs <- unlist(as.data.frame(t(coef(fit))))
	gridmod <- expand.grid(names(fit$model)[1], names(coefs))
	gridmod$coe <- coefs
	gridmod$LI  <- round(ci[, 1], 5)
	gridmod$LS  <- round(ci[, 2], 5)
	gridmod$p.valor <- pv
	gridmod$AIC <- round(AIC,5)
	gridmod$AIC <- round(AIC,5)
	gridmod$AIC = round(AIC,5)
	gridmod$R2 = round(R2, 5)
	gridmod$R2Parcial <- round(r2parcial, 5)
	gridmod$RMSE = round(RMSE,5)
	gridmod$residtest = round(residtest,5)
	gridmod$pF = round(fstats, 5)
	gridmod$PRESS = round(press[2], 5)
	gridmod$PRESSr2 =  round(press[1], 5)
	#gridmod$formula = formulas

	names(gridmod) <- c("v.resposta","explicativa", "coeficiente", "ic_2.5", "ic_97.5", "p.valor", "aic","r2","parcial.r2","rmse","resid.test","f.test", "press", "press.r2")

	## Seleciona apenas modelos com R2 > um valor definido na entrada
	out <- gridmod[(gridmod$r2 > r2 & gridmod$resid.test < residteste), ]
	out <- out[order(out$r2, decreasing = TRUE), ]
	return(out)
}


#####################################
# Automated model selection
# Author      : Joris Meys
# version     : 0.2
# date        : 12/01/09
#####################################
#CHANGE LOG
# 0.2   : check for empty scopevar vector
#####################################

# Function has.interaction checks whether x is part of a term in terms
# terms is a vector with names of terms from a model
has.interaction <- function(x, terms, ...){
    out <- sapply(terms,function(i){
        sum(1-(strsplit(x,":")[[1]] %in% strsplit(i,":")[[1]]))==0
    })
    return(sum(out)>0)
}



# Function Model.select
# model is the lm object of the full model
# keep is a list of model terms to keep in the model at all times
# sig gives the significance for removal of a variable. Can be 0.1 too (see SPSS)
# verbose=T gives the F-tests, dropped var and resulting model after
model.select <- function(model, keep, sig = 0.10, verbose=FALSE, ...){
	## Fixar semente aleat?ria para possibilidade de replica de codes
	  set.seed(47743)
      counter=1
      # check input
      if(!is(model,"lm")) stop(paste(deparse(substitute(model)),"is not an lm object\n"))
      # calculate scope for drop1 function
      terms <- attr(model$terms,"term.labels")
      if(missing(keep)){ # set scopevars to all terms
          scopevars <- terms
      } else{            # select the scopevars if keep is used
          index <- match(keep,terms)
          # check if all is specified correctly
          if(sum(is.na(index))>0){
              novar <- keep[is.na(index)]
              warning(paste(
                  c(novar,"nao foi encontrada no modelo",
                  "\nEstes termos foram ignorados na selecao de modelos."),
                  collapse=" "))
              index <- as.vector(na.omit(index))
          }
          scopevars <- terms[-index]
      }

      # Backward model selection :

      while(T){
          # extract the test statistics from drop.
          test <- drop1(model, scope=scopevars,test="F") # test="Chisq"

          if(verbose){
              cat("------------- STEP ",counter,"-------------\n",
              "\nEstatistica de exclusao : \n")
              print(test)
          }

          pval <- test[,dim(test)[2]]

          names(pval) <- rownames(test)
          pval <- sort(pval,decreasing=T)

          if(sum(is.na(pval))>0) stop(paste("Modelo",
              deparse(substitute(model)),"invalido. Verifique se todos os coeficientes foram estimados."))

          # check if all significant
          if(pval[1]<sig) break # stops the loop if all remaining vars are sign.

          # select var to drop
          i=1
          while(T){
              dropvar <- names(pval)[i]
              check.terms <- terms[-match(dropvar,terms)]
              x <- trycatch_w_e(has.interaction(dropvar,check.terms))$value
			  #if(class(x)[1] != "numeric") x <- 0
              if(x){i=i+1;next} else {break}
          } # end while(T) drop var

          if(pval[i]<sig) break # stops the loop if var to remove is significant

          if(verbose){
             cat("\n--------\nTermos removidos no STEP",counter,":",dropvar,"\n--------\n\n")
          }

          #update terms, scopevars and model
          scopevars <- scopevars[-match(dropvar,scopevars)]
          terms <- terms[-match(dropvar,terms)]

          formul <- as.formula(paste(".~.-",dropvar))
          model  <- update(model, formul)

          if(length(scopevars)==0) {
              warning("Todas as variaveis foram removidas do modelo.\n",
              "Nenhum modelo especificado.")
              return()
          }
          counter=counter+1
      } # end while(T) main loop
      return(model)
}


	## Fun??o geral para o c?lculo dos n melhores modelos entre as explicativas e as preditoras baseado em WALD e AIC

exaustive_lm <- function(vresposta, explicativas, dados, r2=0.50, nmod=5, residteste=0.5, sig=0.10, verbose=FALSE, regsimples = FALSE, ...){
	## Cria lista de modelos
	lista_de_modelos <- lapply(seq_along((explicativas)), function(n) {
		lado_esquerdo_da_funcao  <- vresposta
		lado_direiro_da_funcao <- apply(X = combn(explicativas, n), MARGIN = 2, paste, collapse = " + ")
		paste(lado_esquerdo_da_funcao, lado_direiro_da_funcao, sep = "  ~  ")
	})

	if (regsimples == TRUE) {
	## Converte a lista para vetor
		vetor_modelos <- as.character(unlist(lista_de_modelos[1]))
	} else {
		vetor_modelos <- as.character(unlist(lista_de_modelos))
	}

	## Ajusta modelos lineares para todos os elementos do vetor de modelos
	# Faz o otimiza??o dos modelos lms da etapa anterior, removendo todas as vari?veis com WALD menor que 5% (ajustes n?o significativso)
	lista_de_ajustes <- lapply(vetor_modelos,
		function(x)
		{
			formula <- as.formula(x)
			lmtemp  <- lm(formula, data=dados)
			tmp     <- trycatch_w_e(model.select(lmtemp, sig=sig, verbose=verbose))$value
			if (class(tmp)[1] == "lm") {
				return(tmp)
			}
		}
	)

	# Remove poss?veis valores missings
	lista_de_ajustes <- Filter(Negate(function(x) is.null(unlist(x))), lista_de_ajustes)

	# Remove modelos duplicados
	objects <- lista_de_ajustes[!duplicated(vapply(lista_de_ajustes, function(m) deparse(formula(m)), FUN.VALUE = "a"))]

	# Realiza a regra de sele??o por R2 e por Qualidade do res?duo
	objects <- lapply(objects,
		fun <- function(x, ...) {
			rsquared   <- summary(x)$r.squared
			resid.test <- analise_residuos_glm(x)
			if(rsquared > r2 & resid.test < residteste) {
			return(x)
			} else NULL
		}
	)
	# Remove poss?veis valores missing ap?s filtro
	objects <- Filter(Negate(function(x) is.null(unlist(x))), objects)

	formulas  <- sapply(objects, function(m) deparse(formula(m)))
	vresposta <- sapply(objects, function(m) names(m$model)[1])

	modelos   <- try(ldply(objects, glm_extract_list_ci, r2=r2, residteste=residteste))
  if (class(modelos) !="try-error") {
	modelos   <- na.omit(modelos)
	modelos   <- modelos[order(modelos$r2, decreasing = TRUE), ]
	modelos$id <- NULL
	# Realiza split do data.frame para resgatar os 5 ou menos primeiros modelos pelo R^2

	dfsplit <- split(modelos, modelos$r2)

	# Trata n?mero de modelos a deixar para cada vari?vel resposta
	if (length(dfsplit) < nmod) {
		nmod <- length(dfsplit)
	} else {
		dfsplit <- dfsplit[length(dfsplit):(length(dfsplit)-(nmod-1))]
	}

	# Agrupa novamente os dados
	modelos <- do.call("rbind", dfsplit)
	rownames(modelos) <- c()

	# Corrige ID
	modelos$id <-  cumsum(!duplicated(modelos[, c("v.resposta","r2")]))

	out <- c()
	out$objects  <- objects
	out$formulas <- formulas
	out$vresposta<- vresposta
	out$modelos  <- modelos

	if (length(objects) > 0) {
		cat("Log: Total de modelos avaliados na iteracao:", length(objects), "\n")
	} else {
		cat("Log: Total de modelos avaliados na iteracao:", 0, "\n")
	}
} else {
    cat("Log: Nao foram encontrados modelos!\n")
}
	return(out)
}

cmis_corr <- function(tabela_y, tabela_x=NULL, method = "pearson", digits = 4, verbose = FALSE, corsimples = FALSE, ...) {
	if(is.null(tabela_x)) {
		tab_cor <- round(cor(tabela_y, use = "complete.obs", method=method, ...), digits)
		dados <- tabela_y
	}
	else if(is.null(tabela_y)){
		tab_cor <- round(cor(tabela_x, use = "complete.obs", method=method, ...), digits)
		dados <- tabela_x
	} else {
		tab_cor <- round(cor(tabela_y, tabela_x, use = "complete.obs", method=method, ...), digits)
		dados <- cbind(tabela_y, tabela_x)
	}
	if (corsimples) { # Tabela de correlação simples nxn
	tab_cor <- round(cor(dados, use = "complete.obs", method=method, ...), digits)
	return(tab_cor)	
	} else {
	tab_cor[lower.tri(tab_cor, diag=TRUE)] <- NA 	# Prepara para dropar correlações duplicadas
	tab_cor <- as.data.frame(as.table(tab_cor))  	# Transforma em uma tabela de três colunas
	tab_cor <- na.omit(tab_cor)  					# Elimina os valores missings
	names(tab_cor) <- c("Y","X", "Cor")
	tab_cor <- tab_cor[order(tab_cor$Y, -abs(tab_cor$Cor)),]  # Ordena pelas maiores correlações
	rownames(tab_cor) <- 1:nrow(tab_cor)	
	## Estimar e anexar o R-Quadrado à relação
	r2 <- nm <- c()
	for (i in 1:nrow(tab_cor)) {
		formu <- as.formula(paste(as.character(unlist(tab_cor[i,1:2])), collapse="~"))
		mod <- lm(formu, data=dados)
		r2[[i]] <- round(lm_press_stat(mod), digits)
		nm[i]   <- deparse(formu)
		if(verbose) print(summary(mod))
	}	
	dcal <- do.call("rbind", r2)
	out  <- data.frame(tab_cor, R2=dcal$r.squared, R2Ajustado=dcal$adj.r.squared, R2PRESS = dcal$press.r.squared)
	rownames(out) <- as.character(nm)
	return(out)
	}
}

descritiva <- function (x, basic = TRUE, desc = TRUE, norm = FALSE, p = 0.95, dig = 6) {
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

identico <- function(x, interp=FALSE) {

	if (interp) {
	#Caso haja valores missing: Interpolação linear!
		if (any(is.na(x))) {
		
			# Média nos extremos
			if (is.na(x[1])) x[1] <- mean(x, na.rm=TRUE)
			if (is.na(x[length(x)])) x[length(x)] <- mean(x, na.rm=TRUE)
			
			# Interpolação linear nos outros pontos.
			x <- na.interp(x)
		}
	}
	# Comprimento dos dados
	ntot <- length(x)	
	ndup <- length(x[duplicated(x)])
	sdto <- sd(x)
	
	## Remover métricas com mais de 90% repetidos ou com desvio padrão muito baixo!
	ifelse (ndup/ntot > 0.90 | sd(x) < 0.00001 | is.logical(x), 1, 0) 
	# 0 = OK, 1 = Deleta
}

ordem_diferencas <-  function(x, alpha = 0.05, plot = TRUE, ...) {
    if (class(x) == "numeric") {
    	cat("Serie numerica, convertendo em ts com frequency = 2\n")
    	x <- ts(x, frequency = 2)}
    ns <- nsdiffs(x)
    if(ns > 0) {
      xstar <- diff(x, lag=frequency(x), differences=ns)
    } else {
      xstar <- x
    }
    nd <- ndiffs(xstar, alpha=alpha)
    if(nd > 0) {
      xstar <- diff(xstar,differences=nd)
    }
    if (plot & nd > 0){
      par(mfrow=c(2,1))
      plot(x, col=2, ylab="Serie antes", xlab="Tempo", main="")
      plot(xstar, col=3, ylab="Serie depois", xlab="Tempo", main="")
      par(mfrow=c(1,1))
    }
    return(list(antes = x, depois = xstar, ndifs = nd))
  }

## Calcula estatisticas descritivas de um vetor, tabela ou matriz de dados
Desc <- function(dados, nivel=0.95, tipoci = "basic", nsimu = 500, dig=2) {
 
 tipoci <- match.arg(tipoci, c("norm","basic", "stud", "perc", "bca"))
 nivel <- nivel
 
  aux <- function(x,...){  
    x <- na.omit(x)
    minimo  <- media <- mediana <- maximo <- desviop <- meddesv <- na <- null <- B <- cv <- NA 
    minimo  <- min(x, na.rm=TRUE)
    media   <- mean(x, na.rm=TRUE)
    mediana <- median(x, na.rm=TRUE)
    maximo  <- max(x, na.rm=TRUE)
    desviop <- sd(x, na.rm=TRUE)
    meddesv <- media+desviop
    na      <- x[is.na(x)]
    null    <- x[is.null(x)]
    cv <- desviop/media
    
	out <- as.data.frame(t(round(c(Mediana = mediana, DesvP = desviop, "Media+DesvP" = meddesv, Nulos = na + null, Min = minimo, Max = maximo, CV = cv), dig)))

	# Funcao para intervalode confiança da média
    cimean<- function(x, i, ...) {
		m <- mean(x[i])
		n <- length(i)
		v <- (n-1)*var(x[i])/n^2
		c(m, v)
    } 
	da.boot <- boot(x, cimean, R = nsimu)
	
    B <- boot.ci(da.boot, conf = nivel, type = tipoci)
	if (!is.null(B)){
		if (tipoci == "norm") {Li <- B$normal[2];Ls <- B$normal[3]
		} else if(tipoci == "basic") {Li <- B$basic[4]; Ls <- B$basic[5]		
		} else if(tipoci == "stud") {Li <- B$student[4]; Ls <- B$student[5]		
		} else if(tipoci == "perc") {Li <- B$percent[4]; Ls <- B$percent[5]		
		} else {Li <- B$bca[4]; Ls <- B$bca[5]
		}
	} else {Li <- Ls <- media}
	
    media <- paste(round(media,dig), "(", round(Li,dig), ";", round(Ls,dig),")", sep="")
	
	nm <- c(paste("Media(IC", 100*nivel, "%)", sep=""), names(out))
	
	out <- as.data.frame(cbind(media, out))
	colnames(out) <- nm	
    return(out)
  }  
  dados <- as.data.frame(dados)
  out <- do.call("rbind", 
	lapply(dados, function(X) {if(!is.numeric(X)) NULL else aux(X)}))
  return(out)
}

aggreg <- function(daxts, FUN, freq = "daily", dig = 4, ...) {
	stopifnot(require(xts))
	freq <- match.arg(freq, c("hourly","daily","weekly","monthly","quarterly","yearly"))

	apply.hourly <- function(x, FUN, ...) {
		ends <- endpoints(x, 'hours', 1)
		period.apply(x, ends, FUN=FUN)
	}
	
	fun <- function(freq, ...) {
	  switch(freq,
		daily     = apply.daily,
		weekly    = apply.weekly,
		monthly   = apply.monthly,
		quarterly = apply.quarterly,
		yearly    = apply.yearly,
		hourly    = apply.hourly)
	}

	Fun <- fun(freq)
	La <- lapply(daxts, function(X) Fun(X, FUN=FUN))
	Dc <- do.call("cbind.xts", La)
	return(round(Dc, dig))
}

###################################################
## Funções do pacote lmtest
###################################################

## Teste de Breusch-Pagan
bptest <- function (formula, varformula = NULL, studentize = TRUE, data = list()) 
{
    dname <- paste(deparse(substitute(formula)))
    if (!inherits(formula, "formula")) {
        X <- if (is.matrix(formula$x)) 
            formula$x
        else model.matrix(terms(formula), model.frame(formula))
        y <- if (is.vector(formula$y)) 
            formula$y
        else model.response(model.frame(formula))
        Z <- if (is.null(varformula)) 
            X
        else model.matrix(varformula, data = data)
    }
    else {
        mf <- model.frame(formula, data = data)
        y <- model.response(mf)
        X <- model.matrix(formula, data = data)
        Z <- if (is.null(varformula)) 
            X
        else model.matrix(varformula, data = data)
    }
    if (!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% 
        row.names(X))))) {
        allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
        X <- X[allnames, ]
        Z <- Z[allnames, ]
        y <- y[allnames]
    }
    k <- ncol(X)
    n <- nrow(X)
    resi <- lm.fit(X, y)$residuals
    sigma2 <- sum(resi^2)/n
    if (studentize) {
        w <- resi^2 - sigma2
        fv <- lm.fit(Z, w)$fitted
        bp <- n * sum(fv^2)/sum(w^2)
        method <- "studentized Breusch-Pagan test"
    }
    else {
        f <- resi^2/sigma2 - 1
        fv <- lm.fit(Z, f)$fitted
        bp <- 0.5 * sum(fv^2)
        method <- "Breusch-Pagan test"
    }
    names(bp) <- "BP"
    df <- ncol(Z) - 1
    names(df) <- "df"
    RVAL <- list(statistic = bp, parameter = df, method = method, 
        p.value = pchisq(bp, df, lower.tail = FALSE), data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

## Teste de Durbin-Watson
dwtest <- function (formula, order.by = NULL, alternative = c("greater", 
    "two.sided", "less"), iterations = 15, exact = NULL, tol = 1e-10, 
    data = list()) 
{
    dname <- paste(deparse(substitute(formula)))
    alternative <- match.arg(alternative)
    if (!inherits(formula, "formula")) {
        if (!is.null(w <- weights(formula))) {
            if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w))))) 
                stop("weighted regressions are not supported")
        }
        X <- if (is.matrix(formula$x)) 
            formula$x
        else model.matrix(terms(formula), model.frame(formula))
        y <- if (is.vector(formula$y)) 
            formula$y
        else model.response(model.frame(formula))
    }
    else {
        mf <- model.frame(formula, data = data)
        y <- model.response(mf)
        X <- model.matrix(formula, data = data)
    }
    if (!is.null(order.by)) {
        if (inherits(order.by, "formula")) {
            z <- model.matrix(order.by, data = data)
            z <- as.vector(z[, ncol(z)])
        }
        else {
            z <- order.by
        }
        X <- as.matrix(X[order(z), ])
        y <- y[order(z)]
    }
    n <- nrow(X)
    if (is.null(exact)) 
        exact <- (n < 100)
    k <- ncol(X)
    res <- lm.fit(X, y)$residuals
    dw <- sum(diff(res)^2)/sum(res^2)
    Q1 <- chol2inv(qr.R(qr(X)))
    if (n < 3) {
        warning("not enough observations for computing a p value, set to 1")
        pval <- 1
    }
    else {
        if (exact) {
            A <- diag(c(1, rep(2, n - 2), 1))
            A[abs(row(A) - col(A)) == 1] <- -1
            MA <- diag(rep(1, n)) - X %*% Q1 %*% t(X)
            MA <- MA %*% A
            ev <- eigen(MA)$values[1:(n - k)]
            if (any(Im(ev) > tol)) 
                warning("imaginary parts of eigenvalues discarded")
            ev <- Re(ev)
            ev <- ev[ev > tol]
            pdw <- function(dw) .Fortran("pan", as.double(c(dw, 
                ev)), as.integer(length(ev)), as.double(0), as.integer(iterations), 
                x = double(1), PACKAGE = "lmtest")$x
            pval <- switch(alternative, two.sided = (2 * min(pdw(dw), 
                1 - pdw(dw))), less = (1 - pdw(dw)), greater = pdw(dw))
            if (is.na(pval) || ((pval > 1) | (pval < 0))) {
                warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
                exact <- FALSE
            }
        }
        if (!exact) {
            if (n < max(5, k)) {
                warning("not enough observations for computing an approximate p value, set to 1")
                pval <- 1
            }
            else {
                AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), 
                  ncol = k)
                AX[1, ] <- X[1, ] - X[2, ]
                AX[n, ] <- X[n, ] - X[(n - 1), ]
                XAXQ <- t(X) %*% AX %*% Q1
                P <- 2 * (n - 1) - sum(diag(XAXQ))
                Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% 
                  Q1)) + sum(diag(XAXQ %*% XAXQ))
                dmean <- P/(n - k)
                dvar <- 2/((n - k) * (n - k + 2)) * (Q - P * 
                  dmean)
                pval <- switch(alternative, two.sided = (2 * 
                  pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE)), 
                  less = pnorm(dw, mean = dmean, sd = sqrt(dvar), 
                    lower.tail = FALSE), greater = pnorm(dw, 
                    mean = dmean, sd = sqrt(dvar)))
            }
        }
    }
    alternative <- switch(alternative, two.sided = "true autocorrelation is not 0", 
        less = "true autocorrelation is less than 0", greater = "true autocorrelation is greater than 0")
    names(dw) <- "DW"
    RVAL <- list(statistic = dw, method = "Durbin-Watson test", 
        alternative = alternative, p.value = pval, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

###################################################
## Funcões do pacote TSA
###################################################

## Teste de Portmanteau Tests for Fitted ARIMA

LB.test <- function (model, lag = 12, type = c("Ljung-Box", "Box-Pierce"), 
    no.error = FALSE, omit.initial = TRUE) 
{
    x = residuals(model)
    d1 = length(model$mod$Delta)
    if (omit.initial) 
        x = window(x, start = time(x)[d1 + 1])
    narma = sum(eval(model$arma)[1:4])
    if (is.null(model$call$fixed)) 
        nparm = narma
    else nparm = sum(is.na(eval(model$call$fixed)[1:narma]))
    if (lag <= nparm) {
        if (no.error) 
            return(list(p.value = NA))
        else stop("number of lags is less than the number of parameters: increase the lag")
    }
    if (NCOL(x) > 1) 
        stop("x is not a vector or univariate time series")
    DNAME <- paste("residuals from ", deparse(substitute(model)))
    type <- match.arg(type)
    cor <- acf(x, lag.max = lag, plot = FALSE, na.action = na.pass, 
        drop.lag.0 = TRUE)
    n <- sum(!is.na(x))
    PARAMETER <- lag - nparm
    obs <- cor$acf[1:lag]
    if (type == "Box-Pierce") {
        METHOD <- "Box-Pierce test"
        STATISTIC <- n * sum(obs^2)
        PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
    }
    else {
        METHOD <- "Box-Ljung test"
        STATISTIC <- n * (n + 2) * sum(1/seq(n - 1, n - lag) * 
            obs^2)
        PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
    }
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    names(lag) <- "lag"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME, lag = lag), 
        class = "htest")
}

###################################################
## Funções do pacote tseries
###################################################

## Teste de JarqueBera
jarque.bera.test <- function (x) 
{
    if (NCOL(x) > 1) 
        stop("x is not a vector or univariate time series")
    if (any(is.na(x))) 
        stop("NAs in x")
    DNAME <- deparse(substitute(x))
    n <- length(x)
    m1 <- sum(x)/n
    m2 <- sum((x - m1)^2)/n
    m3 <- sum((x - m1)^3)/n
    m4 <- sum((x - m1)^4)/n
    b1 <- (m3/m2^(3/2))^2
    b2 <- (m4/m2^2)
    STATISTIC <- n * b1/6 + n * (b2 - 3)^2/24
    names(STATISTIC) <- "X-squared"
    PARAMETER <- 2
    names(PARAMETER) <- "df"
    PVAL <- 1 - pchisq(STATISTIC, df = 2)
    METHOD <- "Jarque Bera Test"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME), 
        class = "htest")
}
