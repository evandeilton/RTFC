\name{exaustive_lm}
\alias{exaustive_lm}
\title{
Ajuste e escolha de modelos lineares simples e multiplos automaticamente
}
\description{
Utiliza diversos critérios estatísticos para ajustar e escolher os melhore modelos lm de através de um conjunto de séries temporais de forma automatizada
}
\usage{
exaustive_lm(vresposta, explicativas, dados, r2, nmod, residteste, sig, verbose,
			regsimples = FALSE, ...)
}
\arguments{
    \item{vresposta}{
  Variável reposta
  }
    \item{explicativas}{
  Vetor de variáveis explicativas
  }
    \item{dados}{
  Data.frame no formato data.frame(data, var1, var2, var3, ..., varn)
  }
    \item{r2}{
  Threshold do coeficiente de determinação ou r-quadrado
  }
    \item{nmod}{
  Número de melhores modelos a retornar para a variável explicada em função das explicativas
  }
    \item{residteste}{
  Threshold do indicador de qualidade do resíduo gerado por \code{analise_residuos_glm}
  }
    \item{sig}{
  Threshold de decisão sobre o coeficiente, se ele entra ou não no modelo. O Padrão é 0.05 mas atualmente usa-se 0.10.
  }
    \item{verbose}{
  Se TRUE, exibe log de todas as iterações
  }
    \item{regsimples}{
  Se TRUE, ajusta apenas modelos de regressão linear simples, caso contrário, simples e múltiplos
  }
    \item{\dots}{
  Passagem de parâmetros internos
  }
}
\details{
Esta função utiliza a aliminação stepwise por p-valor para chegar aos melhores modelos para cada variável explicativa e retorna alguns deles segundo existência ou definição do threshold nmod (número de modelos). Esta função integra várias outras do processo de ajuste, análise de resíduos e decisão final. A saída contém os melhores modelos em formato de objeto lm, assim como data.frame com as principais estatísticas no formato para o banco de dados. O número de modelos cresce na escala de 2^n-1 para cada variavel explicativa adicional. Devido a isto, este procedimento se torna muito lento e consome muita memória para conjuntos de variáveis muito grande.
}

\value{
Lista contendo os seguintes elementos
  \item{objects}{Lista dos objetos lm resultados da análise, apenas os melhores escolhidos}
  \item{formulas}{As formulas matemáticas e cada modelo}
  \item{vresposta}{Variável resposta}
  \item{modelos}{data.frame com os modelos em formado de tabela para o banco de dados}
}

\references{
Alan Miller "Subset Selection in Regression" Chapman \& Hall

Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, http://www.stanford.edu/~hastie/Papers/glmnet.pdf

Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
http://www.jstatsoft.org/v33/i01/

Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13
http://www.jstatsoft.org/v39/i05/

Buckland (1997) Model Selection: an Integral Part of Inference. Biometrics 10:41

Burnham & Anderson (2002) Model Selection and Multimodel Inference: an Information Theoretic Approach Calcagno \& de Mazancourt 2010 J. Stat. Soft. v34 i12. See http://www.jstatsoft.org/v34/i12
}
\author{
LOPES, J. L.
}

\seealso{
Esta função depende de \code{\link{model.select}}, \code{\link{analise_residuos_glm}} e  \code{\link{glm_extract_list_ci}}
}
\examples{
# Objeto com dados para modelagem (demorado!)

#data(mensal)
#glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
#damod <- glmdata$ds_modelagem
#vresposta <- glmdata$respostas[1]
#explicativas <- glmdata$explicativas

## Infraestrutura em função de negócio
#modelos <- exaustive_lm(vresposta, explicativas, damod, r2=0.50, nmod=10,
#			residteste = 0.3, sig = 0.05, verbose=FALSE)
}