####################################################################
### O Objetivo do script é gerar gráficos para RDA
###
### Arthur V. Rodrigues - criado em dez/2018 

################################
################################
## Carregue pacote e funções ##
################################
################################

library(vegan)

arw.config.rda <- function(arw.size, show.axis){
  
  bip <- arws.var[,show.axis]
  mul <- ordiArrowMul(bip , fill = arw.size)
  bip.scl <- bip*mul
  labs <- row.names(bip)
  bip.lab <- ordiArrowTextXY(bip.scl, choices = show.axis, 
                             rescale = FALSE, labels = labs)
  
  list( bip.scl = bip.scl, labs = labs, bip.lab = bip.lab)
  
}

def.lim.rda <- function(show.axis){
  r1 <- range(pt.sites[,show.axis[1]])
  r1 <- max(abs(r1))*1.1
  r2 <- range(pt.sites[,show.axis[2]])
  r2 <- max(abs(r2))*1.1
  
  list(x = c(-r1,r1), y = c(-r2,r2))
}


#######################
## Carregue os dados ##
spp <- read.csv("teste_spp.csv", row.names = 1)
env <- read.csv("teste_env.csv", row.names = 1)

###################################################
## Transforma dados de composição para hellinger ##
spp.hel<-decostand(spp,"hell") 
dim(spp.hel) 

##########################################
## Diminui tamanho do nome das espécies ##
names(spp.hel) <- make.cepnames(names(spp.hel))

#######################################
## Padronize as variaveis ambientais ##
env.std <-decostand(env,"stand") 

###################
## Execute a PCA ##
###################
rda.res <- rda(spp.hel, env.std)
rda.summ <- summary(rda.res) ## resultados resumidos
pt.sites <- rda.summ$sites ## coordenadas gráficas para pontos dos sitios
arws.var <- rda.summ$biplot ## coordenadas gráficas para variaveis ambientais

#############################################
## Quais eixos da RDA irão para o gráfico? ##
show.axis <- c(x=1, y=3)
lims <- def.lim.rda(show.axis) # definindo limites do gráfico

###################################################################################
## Define a prioridade para mostrar o nome da espécie de acordo com a abundancia ##
stems <- colSums(spp)

#####################################
## Há agrupamentos para os sitios? ##

## Se sim:
# insira os grupos como factor em 'gr'
# defina as cores em 'color'

#gr <- as.factor(rep(c("g1", "g2", "g3"), each = 26)) # insira os grupos'
#color <- c("green", "black", "yellow") # escolha as cores

## Se não
gr <- as.factor(rep("g1", nrow(spp)))
color <- "lightgreen" ## escolha as cores
##############################
## Plote as espécies da RDA ##
plot(rda.res, dis = "sp", 
     type = "n", 
     xlim = lims$x, 
     ylim = lims$y,
     xlab = paste("RDA", show.axis[1]),
     ylab = paste("RDA", show.axis[2]))
orditorp(rda.res, dis = "sp", lab = names(spp.hel), # Nomes a serem mostrados
         priority = stems, # define a prioridade 
         pcol = "gray", # defini cor do simbolo para espécies
         pch = "+") # define simbolo para espécies

################################
# Plote pontos para os sítios ##
points(pt.sites[,1:2], 
       pch = 16,
       col = color[gr],
       cex = 0.7)

###################
## Gerando setas ##

# Selecione as variaveis a serem mostradas no 
# gráfico a partir de um limite de correlação
# da variável com os eixos mostrados
arw.text <- arw.config.rda(arw.size = 0.7, ## proporção das setas
                       show.axis)

arrows(0, 0, arw.text$bip.scl[,1], arw.text$bip.scl[,2], 
       length = 0.1, # Altere o tamanho da ponta da seta
       col = "blue", # Altere a cor da seta
       lwd = 1) # espessura da seta
text(arw.text$bip.lab, 
     labels = arw.text$labs, # Altere o nome das variáveis
     col = "black", # Altere a cor do texto
     font = 2, # 1- normal, 2-Negrito, 3-Itálico
     cex = 0.7) # Altere tamanho do texto
legend("bottomright", # posição da legenda
       legend = levels(gr), 
       pch = 16, ## Altere o simbolo
       col = color,
       cex = 0.7, ## Altere o tamanho da legenda
       bty = "n") ## com borda bty = "o" / sem borda bty = "n"

