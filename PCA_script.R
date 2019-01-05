####################################################################
### O Objetivo do script é gerar gráficos para PCA
###
### Arthur V. Rodrigues - criado em dez/2018 

################################
################################
## Carregue pacote e funções ##
################################
################################
library(vegan)

arw.config <- function(cut.arw, arw.size, show.axis){
  arw.name <- sqrt(pca$rotation[,show.axis[1]]^2) >= cut.arw | 
    sqrt(pca$rotation[,show.axis[2]]^2) >= cut.arw
  bip <- scores(pca, choices = show.axis, display = "species")[arw.name,]
  mul <- ordiArrowMul(bip , fill = arw.size)
  bip.scl <- bip*mul
  labs <- row.names(bip)
  bip.lab <- ordiArrowTextXY(bip.scl, choices = show.axis, 
                           rescale = FALSE, labels = labs)
  
  list( bip.scl = bip.scl, labs = labs, bip.lab = bip.lab)
  
}

def.lim <- function(show.axis){
  r1 <- range(pca$x[,show.axis[1]])
  r1 <- max(abs(r1))*1.5
  r2 <- range(pca$x[,show.axis[2]])
  r2 <- max(abs(r2))*1.5
  
  list(x = c(-r1,r1), y = c(-r2,r2))
}

#######################
## Carregue os dados ##
pca.data <- read.csv("species-level_trait.csv",
                     header=T, # tem cabeçalho
                     sep=",", # qual separador de colunas
                     dec=".", # qual separador de decimal
                     row.names = 1) # a primeira coluna contem nome das espécies

##################################
## Encurtar nomes das variaveis ##
s.split <- strsplit(names(pca.data), "_")
names(pca.data) <- sapply(s.split, function(x) x[1])

####################################################################
## Para nome de espécie use a função makecepnames do pacote vegan ##
## names(pca.data) <- make.cepnames(names(pca.data))

##################################
## Remova variaveis indesejadas ##
pca.data <- pca.data[,-15] #Retirar variavel binária "DS"

############################
## Padronize as Variáveis ##
pca.data.std <- decostand(na.omit(pca.data), #na.omit() remove linhas que contenham algum NA
                          "stand") 

###################
## Execute a PCA ##
###################
pca <- prcomp(pca.data.std)


#####################################
## Há agrupamentos para os pontos? ##

## Se sim:
# insira os grupos como factor em 'gr'
# defina as cores em 'color'

#gr <- as.factor(rep(c("g1", "g2", "g3"), each = 26)) # insira os grupos'
#color <- c("green", "black", "yellow") # escolha as cores

## Se não
gr <- as.factor(rep("g1", nrow(pca.data.std)))
color <- "green" ## escolha as cores

#############################################
## Quais eixos da PCA irão para o gráfico? ##
show.axis <- c(x=1, y=5)
lims <- def.lim(show.axis) # definindo limites do gráfico

######################################
## Plote os pontos da PCA
plot(pca$x[,show.axis], 
     pch = 16, ## Altere o simbolo
     xlim = lims$x, 
     ylim = lims$y, 
     col = color[gr], ## altere a cor
     xlab = paste("PCA", show.axis[1]),
     ylab = paste("PCA", show.axis[2]), 
     cex = 0.8) ## Altere o tamanho dos pontos

###################
## Gerando setas ##

# Selecione as variaveis a serem mostradas no 
# gráfico a partir de um limite de correlação
# da variável com os eixos mostrados
arw.text <- arw.config(cut.arw = 0.30, ## limite minimo de correlação com os eixos PCA
                       arw.size = 0.5, ## proporção das setas
                       show.axis)

############################################
## Desenhe as setas e os nomes no gráfico ##
arrows(0, 0, arw.text$bip.scl[,1], arw.text$bip.scl[,2], 
       length = 0.1, # Altere o tamanho da ponta da seta
       col = 2, # Altere a cor da seta
       lwd = 1) # espessura da seta
text(arw.text$bip.lab, labels = arw.text$labs, 
     col = 1, # Altere a cor do texto
     font = 2, # 1- normal, 2-Negrito, 3-Itálico
     cex = 0.7) # Altere tamanho do texto

legend("bottomright", # posição da legenda
       legend = levels(gr), 
       pch = 16, ## Altere o simbolo
       col = color,
       cex = 1, ## Altere o tamanho da legenda
       bty = "n") ## com borda bty = "o" / sem borda bty = "n"
