cmis_forecast <- 
function(dados, n_passos_frente, freq, foreplot = FALSE, id_metrica_frequencia=0, primeiro_modelo = 0, parametro_1 = NA, parametro_2 = NA, parametro_3 = NA, parametro_4 = NA, parametro_5 = NA, parametro_6 = NA, parametro_7 = NA, modelo = 0, saida_geral = FALSE, ore_save = FALSE, ...) {

	## Confere tipo dos dados
	#if(any(class(dados) %in% c("data.frame","xts","ts"))) {
	#	stop("Exigido dados no formato (data='dd/mm/aaaa HH:MM:SS', valor=numerico ou \ndata-valor-nome_variavel) em forma de 'data.frame'!\n")
	#}

	# Require packages
	#Rpacks()
	nr <- nrow(as.data.frame(dados))
	if (nr < 2*freq) {
		cat("======================================================================\n| A metrica possui (", nr, ") elementos < que 2*(", freq, ")!\n| Mais dados = melhores projecoes!\n======================================================================\n")
	}

	options(digits = 5)
	################################################################################################################
	# CÁLCULOS DAS PROJEÇÕES E ETAPAS DE SERVIDOR
	################################################################################################################

	out <- c()

	# Trata as datas para gerar séries temporais
	etapa1 <- ts.dados(dados, n_passos_frente, freq)

	# Trata outliers para as métricas
	etapa2 <- trata.outliers.ts(etapa1)

	# Ajusta modelos iniciais com os dados tratados e seleciona o melhor modelo dentre os cinco tipos (ARIMA, ETS, e Holt três tipos) com base em 4 tipos de estatísticas de bondade do ajuste.
	etapa3 <- ajusta_modelos_iniciais(etapa2)

	# Com base no modelo escolhido, salva os resultados para o CMIS
	etapa4 <- salva.resultados_v2(etapa3, id_metrica_frequencia, primeiro_modelo, parametro_1, parametro_2, parametro_3, parametro_4, parametro_5, parametro_6, parametro_7, modelo, foreplot)

	# Objeto final
	out <- etapa4

	# Elementos de interesse
	# Outliers:  tratar data em caracter
	outliers   <- out$outliers
	outliers   <- transform(outliers, data = as.character(data))

	# Projeção: tratar data e valor
	projecao   <- out$projecao
	projecao <- transform(projecao, data = as.character(data), valor = as.numeric(valor))

	# Parametros do modelo escolhido
	par_modelo <- out$par_modelo
	pars <- unlist(c(par_modelo))

	# Juntar dados da projecao com dados de outliers
	temp <- data.frame(matrix(NA, nrow = nrow(outliers), ncol = ncol(projecao)))
	names(temp) <- names(projecao)
	temp[, "data"]  <- outliers$data
	temp[, "valor"] <- outliers$realizado

	projoutliers <- rbind(projecao, temp)	
	projoutliers$parametros <- paste(pars, collapse=":")

	## Salva outliers e parâmetros do modelo no servidor

	if (ore_save == TRUE) {
		ore.save(outliers,   name = paste("ds_metr_hist_",   id_metrica_frequencia, sep=''), overwrite = TRUE)
		ore.save(par_modelo, name = paste("ds_metrica_mod_", id_metrica_frequencia, sep=''), overwrite = TRUE)
	}

	if (saida_geral == TRUE) {
		return(out)
	} else {
		return(projoutliers)
	}
}
