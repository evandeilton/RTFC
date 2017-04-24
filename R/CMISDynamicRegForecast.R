cmis_dynamic_reg_forecast <- function(dados, vresposta, explicativas = NULL,  n_passos_frente, nfolds = 5, sig = 0.05, nivel=95, p_aumento=NULL, trace = FALSE, automatic = TRUE, freq = 7, nvif = Inf, ...) {
# Chamada pacotes
#Rpacks()

	varx <- as.character(explicativas)
	vary <- as.character(vresposta)

	#formulacompleta <- as.formula(paste(c(paste(vary, " ~ 1"), varx), collapse=" + "))

	if (automatic) {# Escolha automatica
	
		r1 <- try(elasticnet_reduction(dados, vary, varx, nfolds = nfolds))

		if (class(r1) == "try-error" | length(r1) == 1) {
			cat("Log: Reducao por Elastic nao foi possivel. Tentando Wald!\n")
			r1 <- try(wald_reduction(dados, vary, varx, sig = sig))
			if (class(r1) == "try-error" | length(r1) == 1) {
				cat("Log: Reducao por Wald nao foi possivel. Tentando STEPWISE!\n")
				r1 <- try(step_reduction(dados, vary, varx, sig = sig, verbose = trace))
				if (class(r1) == "try-error" | length(r1) == 1) {
					cat("Log: Reducao por STEPWISE nao foi possivel. Sem preditoras fortes a retornar!\n")
				}
			}
			varxy <- r1
		} else {
			varxy <- r1 ## Apenas a variavel resposta
		}
	} else {
		cat("Log: Decisao automatica desabilitada!\n")
		cat("Log: Verificando multicolinearidade...\n")
				
		Vif <- my_glm_vif(dados,  vary,  varx, nvif = nvif)
		
		for(i in 1:length(Vif$vif)){
			 cat(paste("Log:", names(Vif$vif[i]), "tem VIF:", round(Vif$vif[i],2),collpase="\n"))
		}	
		
		if (any(round(Vif$vif[i],2) > 20)) {
			cat("\nLog: ATENCAO! Valores acima de 20 podem indicar problemas de multicolinearidade.\nIsto pode gerar modelos fracos ou mal estimados!\n")
		}
		
		varxy <- as.character(c(Vif$vresposta, Vif$vexplicativas))
	}
	
	if(length(varxy) < 1) {
		cat("Log: Forecast sem covariavel. Vegetativo apenas!!\n")
		# Dados para forecast da s?rie resposta
		Dy <- dados[, varxy]

		# ARIMA mais Forecast para a s?rie explicativa sem s?ries preditoras
		arimay  <- auto.arima(Dy, trace = trace)
		forey   <- forecast(arimay , h = n_passos_frente, level=nivel)
		acuraciay = LjungBtest_Acuracia(arimay)

		# Objetos de sa?da da fun??o.
		out <- list(models, forecasts, acuracia)
		class(out) <- "cmisdyn"
		return(out)
	} else {

		Dy <- ts(dados[, varxy[1]], frequency = freq) # Dados para forecast da s?rie resposta
		Dx <- ts(dados[,varxy[-1]], frequency = freq) # Dados para forecast das s?ries preditoras

		## Forecast para x colocando o resultado em um objeto lista e aplicando aumento % antes do forecast
		forex <- lapply(Dx, function(X, ...) {
			fit <- try(auto.arima(X))
			if (class(fit)[1]!="try-error") fit  <- fit
			fit <- try(forecast(fit, h=n_passos_frente, level=nivel)$mean)
			if (class(fit)[1]!="try-error") fit else NA
			}
		)

		# Remove NA's da lista
		forex <- forex[!is.na(forex)]

		# Agrupar as s?ries em um objeto mts para forecasts
		forex <- do.call("cbind", forex)

		# ARIMA mais Forecast para a s?rie explicativa sem s?ries preditoras
		arimay  <- auto.arima(Dy, trace = trace)
		forey   <- forecast(arimay , h = n_passos_frente, level=nivel)
		acuraciay = LjungBtest_Acuracia(arimay)
		# ARIMA mais Forecast para a s?rie explicativa com as s?ries preditoras
		arimayx <- try(auto.arima(Dy, xreg = Dx, trace = trace, ic="aicc"))

		if (class(arimayx)[1] != "try-error") {
			acuraciayx = LjungBtest_Acuracia(arimayx)
			## Aumento percentual em X antes do ser projetado.
			if (!is.null(p_aumento)) {
				forex <- forex*(1+p_aumento)
			}
			## Forecast Y em fun?ao de X1, X2, X3, ..., Xn
			foreyx  <- forecast(arimayx, xreg = forex, h=n_passos_frente, level=nivel)

			# Objetos de sa?da da fun??o.
			models    <- list(arimay = arimay, arimayx = arimayx)
			forecasts <- list(forey = forey, foreyx = foreyx)
			acuracia  <- list(acuraciay = acuraciay, acuraciayx=acuraciayx)
			out <- list(models, forecasts, acuracia)
			class(out) <- "cmisdyn"
			return(out)

		} else {
			cat("Log: Impossivel ARIMA(p,d,q) com preditoras. Retornando forecast simples!\n")
					# Dados para forecast da s?rie resposta
			Dy <- dados[, vresposta]

			# ARIMA mais Forecast para a s?rie explicativa sem s?ries preditoras
			arimay  <- auto.arima(Dy, trace = trace)
			forey   <- forecast(arimay , h = n_passos_frente, level=nivel)
			acuraciay = LjungBtest_Acuracia(arimay)

			# Objetos de sa?da da fun??o.
			models    <- list(arimay = arimay, arimayx = NA)
			forecasts <- list(forey = forey, foreyx = NA)
			acuracia  <- list(acuraciay = acuraciay, acuraciayx=NA)
			out <- list(models, forecasts, acuracia)
			class(out) <- "cmisdyn"
			return(out)
		}
	}
}
