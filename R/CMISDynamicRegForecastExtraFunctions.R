# Require packages
#Rpacks()

step_reduction <- function (dados, vary, varx, sig = 0.10, verbose = FALSE, ..){

	varx <- as.character(varx)
	vary <- as.character(vary)

	cat("Log: Reducao de variaveis por STEPWISE com teste T ANOVA tipo II!\n")
	cat("--------------------------------------------------------\n")
	cat("Variaveis de entrada:\n", paste(varx, collpase="\n"))

	y <- as.numeric(dados[, vary])
	x <- as.matrix(as.data.frame(dados[, varx]))
	formulacompleta <- as.formula(paste(c(paste(vary, " ~ 1"), varx), collapse=" + "))

	Lm <- lm(formulacompleta, data = dados)
	Lm <- try(model.select(Lm, sig = sig, verbose = verbose))

	newvarx <- NA
	if (class(Lm)!="try-error") {
		newvarx <- as.character(names(Lm$model))[-1]
		if (length(newvarx) > 0) {
			newvarx <- as.character(newvarx)
			cat("--------------------------------------------------------\n")
			cat("Operacao finalizada com sucesso: Variaveis fortes escolhidas:\n", paste(newvarx, collpase="\n"))
			cat("--------------------------------------------------------\n")
			return(c(as.character(vary), as.character(newvarx)))
		} else {
			cat("Log: Nenhuma preditora forte foi encontrada!\n")
			return(c(as.character(vary, NA)))
		}
	} else {
		cat("Log: Nenhuma preditora forte foi encontrada!\n")
		return(c(as.character(vary, NA)))
	}
}


wald_reduction <- function (dados, vary, varx, sig = 0.10, ...){

	require("rms")
	varx <- as.character(varx)
	vary <- as.character(vary)

	cat("Log: Reducao de variaveis por Wald!\n")
	cat("--------------------------------------------------------\n")
	cat("Variaveis de entrada:\n", paste(varx, collpase="\n"))

	y <- as.numeric(dados[, vary])
	x <- as.matrix(as.data.frame(dados[, varx]))
	formulacompleta <- as.formula(paste(c(paste(vary, " ~ 1"), varx), collapse=" + "))

	Ols <- try(ols(formulacompleta, data=dados))

	newvarx <- NA
	if (class(Ols)[1] != "try-error") {
		z <- fastbw(Ols, sls = sig, eps=1e-30)
		newvarx <- unlist(dimnames(dados[,z$names.kept]))

		if (length(newvarx) > 0) {
		newvarx <- as.character(newvarx)
		cat("--------------------------------------------------------\n")
		cat("Operacao finalizada com sucesso: Variaveis fortes escolhidas:\n", paste(newvarx, collpase="\n"))
		cat("--------------------------------------------------------\n")
		return(c(as.character(vary), as.character(newvarx)))
		} else {
			cat("Log: Nenhuma preditora forte foi encontrada!\n")
			return(c(as.character(vary, NA)))
		}
	} else {
		cat("Log: Nenhuma preditora forte foi encontrada!\n")
		return(c(as.character(vary, NA)))
	}
}


elasticnet_reduction <- function (dados, vary, varx, nfolds = 5, ...){
	require("glmnet")
	varx <- as.character(varx)
	vary <- as.character(vary)

	cat("Log: Reducao de variaveis por Elastic Net!\n")
	cat("--------------------------------------------------------\n")
	cat("Variaveis de entrada:\n", paste(varx, collpase="\n"))

	y <- as.numeric(dados[, vary])
	x <- as.matrix(as.data.frame(dados[, varx]))

	#formulacompleta <- as.formula(paste(c(paste(vary, " ~ 1"), varx), collapse=" + "))

	## Fixar semente aleat?ria para possibilidade de replica de codes
	set.seed(nfolds)
	nfolds <- nfolds

	cvfit <- try(cv.glmnet(x, y, alpha=1, type.measure = "deviance", nfolds=nfolds))
	newvarx <- NA
	if (class(cvfit) !="try-error"){
		co <- as.matrix(coef(cvfit, s = "lambda.1se"))
		co <- data.frame(coefnames = rownames(co), coef = as.numeric(co))
		co <- co[-1, ] # remove a primeira linha (do intercepto)
		rownames(co) <- NULL
		newvarx <- as.character(co[co$coef > 0, ]$coefnames)

		if (length(newvarx) > 0) {
			newvarx <- as.character(newvarx)
			cat("--------------------------------------------------------\n")
			cat("Operacao finalizada com sucesso: Variaveis fortes escolhidas:\n", paste(newvarx, collpase="\n"))
			cat("--------------------------------------------------------\n")
			return(c(as.character(vary), as.character(newvarx)))
		} else {
			cat("Log: Nenhuma preditora forte foi encontrada!\n")
			return(c(as.character(vary, NA)))
		}
	} else {
		cat("Log: Nenhuma preditora forte foi encontrada!\n")
		return(c(as.character(vary, NA)))
	}
}

LjungBtest_Acuracia <- function(model,...) {
	bt <- Box.test(residuals(model), lag=10, type="Ljung", fitdf=length(coef(model)))

	# Teste de autocorrela??o
	LJungBox <- data.frame(LJB_X_quad = round(bt[[1]], 5), LJB_gl = bt[[2]], LJB_p.valor = round(bt[[3]], 5))
	rownames(LJungBox) <- c()
	acc <- data.frame(round(accuracy(model), 5))

	## Correla??o Serial por Durbin-Watson
	D.Watson <- round(dwtest(model$x~residuals(model))$statistic, 5)

	res <- cbind(acc, LJungBox, D.Watson)
	res
}

plot.cmisdyn <- function(obj, ...) {
	if (class(obj) != "cmisdyn") stop(cat("O objeto informado deve ser da classe 'cmis'\n"))

	forecasts <- obj[[2]]
	forey  <- forecasts$forey
	foreyx <- forecasts$foreyx
	acc <- c(11, 10, 6, 5)
	fun <- function(X){if(X %in% c(-Inf, Inf)) NA else X}
	linhas <- function(fit) {
		lines(fitted(fit), col='green')
		legend("topleft", c("observado","predito"), pch = '_', col = c("black","green"))
	}
	
	if(!is.logical(foreyx) & !is.logical(forey)) {

		## MAPE = Mean absolute percentage error
		## MASE = Mean absolute scaled error

		acy  <- as.matrix(obj[[3]]$acuraciay[acc])
		ac_y <- sapply(acy, fun)
		names(ac_y) <- colnames(acy)

		acyx <- as.matrix(obj[[3]]$acuraciayx[acc])
		ac_yx <- sapply(acyx, fun)
		names(ac_yx) <- colnames(acyx)

		layout(matrix(c(1,1,2,3,3,4), 2, 3, byrow = TRUE))

		plot(foreyx, ylab="ARIMA com regressoras")
		linhas(foreyx)
		
		yy <- barplot(round(ac_yx, 4),  col=3, horiz=TRUE, main="Acuracia", las=1, border = NA)
		text(x = round(ac_yx, 2), y = yy, labels=as.character(round(ac_yx, 2)), xpd=TRUE)

		plot(forey,  ylab="ARIMA")
		linhas(forey)
		
		yy <- barplot(round(ac_y, 4),  col=5, horiz=TRUE, main="Acuracia", las=1, border = NA)
		text(x = round(ac_y, 2), y = yy, labels=as.character(round(ac_y, 2)), xpd=TRUE)
		par(mfrow=c(1,1))
	} else if(!is.logical(forey)) {
		layout(matrix(c(1,1,2,2,3,3), 1, 3, byrow = TRUE))
		acy  <- as.matrix(obj[[3]]$acuraciay[acc])
		ac_y <- sapply(acy, fun)
		names(ac_y) <- colnames(acy)

		plot(forey,  ylab="ARIMA")
		linhas(forey)
		yy <- barplot(round(ac_y, 4),  col=5, horiz=TRUE, main="Acuracia", las=1, border = NA)
		text(x = round(ac_y, 2), y = yy, labels=as.character(round(ac_y, 2)), xpd=TRUE)
		par(mfrow=c(1,1))
	} else {
		cat("Log: Nao foram encontrados forecasts em ", obj)
	}
}

summary.cmisdyn <- function(obj, plot=TRUE, digits = 10){
if (class(obj) != "cmisdyn") stop(cat("O objeto informado deve ser da classe 'cmis'\n"))
	options(digits = digits)
	forey <- obj[[2]]$forey
	foreyx<- obj[[2]]$foreyx
	ac_y  <- obj[[3]]$acuraciay
	ac_yx <- obj[[3]]$acuraciayx

	acuracia <- rbind(ac_yx, ac_y)
	rownames(acuracia) <- c("Y~X", "Y")

	if(plot) plot.cmisdyn(obj)
	print(acuracia)
	summary(foreyx)
	summary(forey)
}

elastic_net_reduction <-
function (dados, vresposta, explicativas, nfolds = 10, ...){
	y <- as.numeric(dados[, vresposta])
	x <- as.matrix(as.data.frame(dados[, explicativas]))
	
	## Fixar semente aleatoria para possibilidade de replica de codes
	set.seed(nfolds)
	nfolds <- nfolds
	cvfit <- try(cv.glmnet(x, y, alpha=1, type.measure = "deviance", nfolds=nfolds))
	
	if (class(cvfit) !="try-error"){	
		co <- as.matrix(coef(cvfit, s = "lambda.1se"))
		co <- data.frame(coefnames = rownames(co), coef = as.numeric(co))
		co <- co[-1, ] # remove a primeira linha (do intercepto)
		rownames(co) <- NULL
		dd <- co[co$coef > 0, ]			
		
		# Caso haja alguma variavel explicativa com coeficiente significativo (ou seja, que foi estimado via CrossValidation)
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
