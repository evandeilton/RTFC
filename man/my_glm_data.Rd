\name{my_glm_data}
\alias{my_glm_data}
\title{
Trata dados para regressoes
}
\description{
Trata dados para entrada no algoritimo das regressões
}
\usage{
my_glm_data(dados, n_passos_frente, freq, prefixo.expl, prefixo.resp, ...)
}

\arguments{
  \item{dados}{
data.frame com variáveis seguindo o modelo (data.frame (DATA, var1, var2, var3, ..., varn)
}
  \item{n_passos_frente}{
Número de passos utilizados na projeção, necessário para a funçao ts.dados trabalhar com datas
}
  \item{freq}{
Frequência dos dados. Mensal (12), diário(7), horário(24)
}
  \item{prefixo.expl}{
String tamanho 3 com o prefixo das variáveis explicativas na forma (aaa). Ex: negocio = (neg), aplicação = (apl). Estes nomes precisam ser definidos no data.frame antes da entrada na função
}
  \item{prefixo.resp}{
String tamanho 4 com o prefixo das variáveis respostas na forma (bbb). Ex: infraestrutura = (inf), hardware = (hdw). Estes nomes precisam ser definidos no data.frame antes da entrada na função.
}
  \item{\dots}{
Passagem de parâmetros
}
}
\details{
Esta função faz uso da função \code{dataframe.outlier} para tratar os dados e retornar um data.frame limpo pronto para modelagem.
}
\value{
Lista com variaveis resposta, explicativas e data.frame com os dados limpos para modelagem.
\item{explicativas}{vetor de variaveis explicativas}
\item{respostas}{vetor de variáveis resposta}
\item{ds_modelagem}{data.frame limpo e pronto para modelagem}
}

\author{
LOPES, J. L.
}

\examples{
data(mensal)
glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
names(glmdata)
}
