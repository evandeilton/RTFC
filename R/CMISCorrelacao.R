cmis_correlacao <- 
function(dados, n_passos_frente = 1, freq=7, nvif = 100, nfolds = 5, nmod = 10, r2 = 0.50, residteste=0.8, sig=0.15, prefixo.expl = "neg", prefixo.resp = "hdw", verbose = FALSE, saida_geral = FALSE, regsimples = TRUE, ...) {

# Require packages
#  Rpacks()

## N?mero de digitos para o ambiente da fun??o.
options(digits=10)

inf_micro <- function() {
	cat(
	"\n======================================================================= \n",
	"Informacoes da maquina \n\n",
	"Sistema Operacional     :", paste(Sys.info()[[1]], Sys.info()[[2]], Sys.info()[[3]], sep=""),"\n",
	"Nome de rede da maquina :", Sys.info()[[4]], "\n",
	"Usuario executando o R  :", Sys.info()[[6]], "\n",
	"Versao do R instalado   :", R.Version()$version.string,
	"\n======================================================================= \n")
}
inf_micro()

## data.frame vazio (contorno de problema no sql
colClasses = c("factor", "factor", rep("numeric", 12), "integer")
col.names  = c("v.resposta","explicativa", "coeficiente", "ic_2.5", "ic_97.5", "p.valor", "aic","r2","parcial.r2","rmse","resid.test","f.test", "press", "press.r2","id")
df_out <- read.table(text = "",colClasses = colClasses, col.names = col.names)
df_out

################################################################################
## Etapa 1: entrada e tratamento dos dados
################################################################################
inicio <- as.character(Sys.time())

etapa1 <- my_glm_data(dados, n_passos_frente, freq, prefixo.expl, prefixo.resp)

explicativas <- as.character(etapa1[[1]])
vresposta <- as.character(etapa1[[2]])
dados_mod <- etapa1[[3]]
if (length(explicativas) < 1) {stop("Nao ha metricas explicativas neste conjunto de dados!")}
final <- as.character(Sys.time())

cat("\n=======================================================================",
"\nEtapa 1: Tratamento dos dados",
"\n=======================================================================\n",
"Processo iniciado por          :", Sys.getenv("USERNAME"),"\n",
"Inicio                         :", inicio, "\n",
"Final                          :", final, "\n",
"Total de metricas explicativas :", length(explicativas), "\n",
"Total de metricas explicadas   :", length(vresposta),
#"Total de modelos possiveis     :", sum(nmods), "\n",
"\n=======================================================================\n")

################################################################################

etapa2 <- modelos <- obj1 <- modelos <- c()
inicio <- as.character(Sys.time())

	cat("\n=======================================================================",
		"\nEtapa 2: Analise exaustiva dos modelos candidatos!\n"
	)

	if (regsimples) {
	cat(
	"\nApenas Regressoes simples y = b0 + b1*X1 serao avaliadas e retornadas!\n",
	"\n======================================================================= \n")
	} else {
	cat(
	"\nRegressoes multiplas y = b0 + b1*X1 + b2*X2 + bn*Xn serao avaliadas.\nModelos simples e tambem multiplos podem retornar!\n",
	"\n======================================================================= \n")
	}
	if (regsimples==TRUE){
		nvif <- Inf
	} else {
		nvif <- nvif
	}
	for (i in 1:length(vresposta)){
		vres <- vresposta[i]
		lambdacv <- NA
		# Avalia variaveis com VIF < vif e remove variaveis com VIF < nvif
		gtemp <- try(my_glm_vif(dados_mod, vres, explicativas, nvif=nvif))
		if(class(gtemp)!="error" & nrow(gtemp$dados) > 0) {
			dtemp <- gtemp$dados
			etemp <- gtemp$vexplicativas
			if (length(etemp) == 1) {
				cat("Log: Analise da Metrica (", vres, ") concluida sem Cross!\n")
				temp <- trycatch_w_e(exaustive_lm(vres, etemp, dtemp, r2, nmod, residteste, sig, verbose, regsimples))$value
			} else if (length(etemp) > 1) {
			# Faz an?lise de CrossValidation via LASSO
			# http://statweb.stanford.edu/~tibs/lasso/lasso.pdf
			# The LASSO (Least Absolute Shrinkage and Selection Operator) is a regression method that involves
      # penalizing the absolute size of the regression coefficients.
			#cat("Log: Analisando regressao multiplas, com CrossValidation!\n")

			obj0 <- try(my_glm_cv(dtemp, vres, etemp, nfolds))
				if(class(obj0)!="try-error") {
					dado <- obj0$dados
					expl <- obj0$vexplicativas
					temp <- trycatch_w_e(exaustive_lm(vres, expl, dado, r2, nmod, residteste, sig, verbose, regsimples))$value
					lambdacv <- trycatch_w_e(obj0$lambdacv)$value
					if(!is.null(lambdacv)){
					cat("Log: Analise da Metrica (", vres, ") concluida!", "Lambda CV:", round(lambdacv,5), "\n")
					}
				}
			} else {
				cat("\n")
			}
			etapa2[[i]]   <- temp$modelos
			modelos[[i]]  <- temp$objects
		}
	}

	################################################################################
	## Etapa 3: Conslida??o dos modelos e estat?sticas de sa?da!
	################################################################################

	# Remove poss?veis valores missings
	etapa2  <- Filter(Negate(function(x) is.null(unlist(x))), etapa2)
	modelos <- Filter(Negate(function(x) is.null(unlist(x))), modelos)

	#ff <- function(x) {if class}
	etapa3 <- trycatch_w_e(do.call("rbind", etapa2))$value

  if (class(etapa3)[1]=="data.frame") {
  	if (ncol(etapa3) > 1) {

  		## Reordena modelos pelo ID para evitar duplicatas na m?trica resposta autal
  		etapa3$id <-  cumsum(!duplicated(etapa3[, c("v.resposta","r2")]))
  		etapa3full <- list(etapa3, modelos)
  		nmelhores <- length(unique(etapa3$aic))
  		final <- as.character(Sys.time())

  		cat("\n=======================================================================\n",
  		"Etapa 2 concluida!\n\n",
  		"Inicio do algoritmo            :", inicio, "\n",
  		"Fim algoritmo                  :", final, "\n",
  		#"Total modelos factiveis  (app) :", sum(nmods), "\n",
  		"Total modelos escolhidos (app) :", max(etapa3$id), "com ", "R2 >", r2, "e residuos <", residteste, "\n",
		"Total reamostragens (CrossVal) :", nfolds,
  		"\n======================================================================= \n\n")

  		# Se desejar a sa?da como caractere pra todas m?tricas.
  		#etapa4 <- data.frame(lapply(etapa3, as.character), stringsAsFactors=FALSE)

  		if (saida_geral == TRUE) {
  			return(etapa3full)
  		} else {
  			return(etapa3)
  		}
  	} else {
  		cat("Log: Nao foram encontrados modelos significativos para este aplicativo!\n")
		return(df_out)
  	}
  } else {
    cat("Log: Nao foram encontrados modelos significativos para este aplicativo!\n")
	return(df_out)
  }
}
