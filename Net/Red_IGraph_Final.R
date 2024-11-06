############################################################
## Tema: Análisis de redes de coexpresión                 ##
## Autor: Olga Andrea Hernandez Miranda, Miranda H        ##
## Fecha: 31/10/2024                                      ##
## Nota: Redes de co-expresión en Igraph y modulos WGCNA  ##      ##
############################################################

# Cargar los paquetes necesarios
library(tidyverse)
library(igraph)
library(plotly)
library(BoolNet)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# Establecer el directorio de trabajo
directorio <- "C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial5/RedVainilla"
setwd(directorio)

# Cargar los datos de expresión diferencial
diff.expr <- fread("Redh_Tot2.csv", header = TRUE)

# Verificar las dimensiones del data frame
cat("Dimensiones de diff.expr:", dim(diff.expr), "\n")

# Convertir a matriz de expresión
allexpr <- as.matrix(diff.expr[, -1])  # Excluir la primera columna (nombres de filas)
row.names(allexpr) <- diff.expr$V1      # Establecer nombres de las filas

# Filtrar por etapas (opcional, si necesitas subconjuntos específicos)
Prepol <- allexpr[grep("Pre-pol", row.names(allexpr)), ]
Pol <- allexpr[grep("Pol", row.names(allexpr)), ]
Postpol <- allexpr[grep("Post-pol", row.names(allexpr)), ]
Fer <- allexpr[grep("Fer", row.names(allexpr)), ]

# Seleccionar la matriz de expresión completa
expr <- allexpr

# Calcular la matriz de correlación, ignorando NA
gene.correlation <- cor(expr, use = "pairwise.complete.obs")

# Eliminar filas y columnas con NA
gene.correlation <- gene.correlation[complete.cases(gene.correlation), complete.cases(t(gene.correlation))]

thresholds<-seq(0.5,0.99,0.01)
mean.connectivities<-NULL
scale.free.R2<-NULL
## Recorrer todos los valores de umbral
for(i in 1:length(thresholds))
{
  ## Se construye una red que corresponde al umbral de correlaci?n especifico que se evalua en este paso
  
  ## Matriz de adyacencia
  current.adjacency <- gene.correlation > thresholds[i]
  ## Red
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## C?lculo de grados de los nodos
  node.degrees <- degree(threshold.network)
  
  ## Se guarda la conectividad promedio de esta red
  mean.connectivities[i] <- mean(node.degrees)
  
  ## Evaluaci?n de la propiedad de escalamiento libre
  h <- hist(node.degrees,main=paste("Histograma de grado de los nodos para el umbral de correlacion",thresholds[i]))
  ## C?lculo de la frecuencia de grados
  degree.frequencies <- table(node.degrees)
  ## Determinaci?n por regresi?n lineal para la transformaci?n logar?tmica de las frecuencias de grado
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  ## Obtenci?n de R cuadrados como una estimaci?n del ajuste a la propiedad de escalamiento libre
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

#El siguiente paso consiste en evaluar el mejor punto de corte para el grado de correlaci?n que mejor ajusta a 
#la red generada a los modelos de redes observados comunmente en la naturaleza. Es importante tener una idea de 
#las propiedades deseadas en la red, por ejemplo, si se sabe cual es la conectividad promedio observada en redes de 
#regulaci?n. As? mismo el criterio del ajuste al modelo de red de libre escala, comunmente observado en la naturaleza,
#puede ser un par?metro ideal para determinar la mejor red.

plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Umbral de correlacion",ylab="Conectividad promedio")
plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.50,0.99),xlab="Umbral de correlacion",ylab=expression('Ajuste al modelo de escala libre (R  '^2*')'))

#---Construcci?n de la red ?ptima---
#Se selecciona un valor para trabajar con esa red tmax va a ser la variable que contiene nuestro umbral de 
#correlaci?n. A continuaci?n se obtiene la red tomando en consideraci?n una conexi?n por cada pareja de genes
#con una correlaci?n mayor al umbral que hemos determinado. Nuestra red quedar? guardada en nuestra variable gene.
#coexpression.network.

tmx<-which(max(scale.free.R2)==scale.free.R2)
tmax<-thresholds[tmx[length(tmx)]]

print(tmax)

adjacency.tmax <- gene.correlation > tmax

for(i in 1:ncol(adjacency.tmax))
{
  #  print(i)
  adjacency.tmax[i,i] <- FALSE
}

gene.coexpression.network <- graph.adjacency(adjacency.tmax, mode="undirected")

#Como una alternativa, se puede exportar una red en formato GML, el cual es importable en diferentes tipos
#de software de visualizaci?n.
gene.coexpression.network
#write.graph(gene.coexpression.network,file="gene_coexpression_network2_TFF_ID.gml",format="gml")

#--Visualizaci?n--
#El siguiente paso es darle un vistazo r?pido a la red generada

plot(gene.coexpression.network,vertex.size=node.degrees, vertex.label.color="black",
     vertex.label.dist=1, vertex.label.cex=0.25, layout=layout_nicely,
     vertex.color="blue", edge.color="grey")

#Y revisamos nuevamente la estructura de la red. Para ello se revisa la distribuci?n del grado en los nodos
#y revisaremos los ajustes lineales que se obtienen.


network.degrees <- degree(gene.coexpression.network)
degree.histogram <- hist(network.degrees,col="royalblue",
                         xlab="Grado de nodo (K)",
                         ylab="Numero de nodos con K enlaces [P(K)]",
                         main = NA)

# Calcular el grado de los nodos
network.degrees <- degree(gene.coexpression.network)

# Convertir los grados a un data frame
degree.data <- data.frame(Degree = network.degrees)

# Crear el histograma usando ggplot2
ggplot(degree.data, aes(x = Degree)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "#440154", alpha = 0.8) + # Ajuste de transparencia aquí
  labs(x = "Grado de nodo (K)",
       y = "Número de nodos con K enlaces [P(K)]",
       title = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 11),  # Ajustar el tamaño de la fuente de los ejes
        axis.title = element_text(size = 11), # Ajustar el tamaño de la fuente de los títulos de los ejes
        plot.title = element_text(size = 11),  # Ajustar el tamaño de la fuente del título del gráfico
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20)) # Agregar márgenes (top, right, bottom, left)

#write.table(data.frame(Name = names(network.degrees), Degree = network.degrees), 
 #           file = "grados_de_nodo_TFF_ID.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

degrees.frequencies <- degree.histogram[["counts"]]
node.degrees <- degree.histogram[["mids"]]
log10.degrees.frequencies <- log10(degrees.frequencies[-3])
log10.node.degrees <- log10(node.degrees[-3])
lm.r <- lm(log10.degrees.frequencies[1:6] ~ log10.node.degrees[1:6])
summary(lm.r)

plot(log10.node.degrees,log10.degrees.frequencies, col = "black", pch=16,
     xlab="log10 grado de nodo (K)",
     ylab="log10 numero de nodos con K enlaces [P(K)]",
     abline(lm.r, col = "blue"))

# Data frame de los datos
data <- data.frame(
  Log10Degree = log10.node.degrees,
  Log10Frequency = log10.degrees.frequencies
)

# Predicciones del modelo para la línea de regresión
predictions <- data.frame(Log10Degree = log10.node.degrees[1:6])
predictions$Log10PredictedFrequency <- predict(lm.r, newdata = predictions)

# Crear el gráfico
ggplot(data, aes(x = Log10Degree, y = Log10Frequency)) +
  geom_point(colour = "black", shape = 16) +  # Dibuja los puntos
  geom_smooth(method = "lm", colour = "#440154", fill = "#440154", alpha = 0.2, se = TRUE) +  # Añade la línea de regresión con sombra morada y transparencia
  xlab("log10 grado de nodo (K)") +
  ylab("log10 número de nodos con K enlaces [P(K)]") +
  theme_classic() +
  theme(axis.text = element_text(size = 11),  # Ajustar el tamaño de la fuente de los ejes
        axis.title = element_text(size = 11), # Ajustar el tamaño de la fuente de los títulos de los ejes
        plot.title = element_text(size = 11), # Ajustar el tamaño de la fuente del título del gráfico
plot.margin = margin(t = 20, r = 20, b = 20, l = 20))
#Prueba de bondad de ajuste para una exponenecial negativa

network.degree.distribution <- degree.distribution(gene.coexpression.network)
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]

#---Parametros de la red---
#A continuaci?n revisamos algunos atributos de la red como el score de hub

## La mayoría de los nodos en la redes libre de escala presentan un número pequeño de vecinos. 
## Sin embargo existen unos pocos nodos destacados que tiene un alto número de vecinos. Este 
## último tipo de nodos se denominan hubs. La función hub.score de igraph que recibe como 
## entrada una red calcula y almacena en el valor vector una puntuación entre uno y cero para 
## cada nodo de la red. Cuánto mayor sea esta puntuación de un nodo más se ajustan sus 
## características a las de un hub de la red.

network.hub.scores <- hub.score(gene.coexpression.network)
hub.score.attributes <-network.hub.scores[["vector"]]

#write.table(data.frame(Name = names(hub.score.attributes), `Hub score` = hub.score.attributes), 
 #          file = "hub_score_TFF_ID.txt", 
  #         sep = "\t", 
   #        row.names = FALSE, 
    #       col.names = TRUE, 
     #      quote = FALSE)

plot_ly(y = network.hub.scores$vector[which(network.hub.scores$vector>0.1)], text=names(network.hub.scores$vector[which(network.hub.scores$vector>0.1)]),type="bar")
barplot(network.hub.scores$vector[which(network.hub.scores$vector>0.1)],las=2, col = "light green")

#La transitividad y el camino promedio
#Nodos y aristas
gene.coexpression.network
#619 5998

transitivity(gene.coexpression.network)
average.path.length(gene.coexpression.network)

## Para comprobar si el coeficiente de agrupamiento en lo suficientemente alto y 
## la longitud media del camino mínimo entre nodos es lo suficientemente 
## pequeña como para considerarla de mundo pequeño es común generar redes libres de escala del 
## mismo orden y tamaño de la estudiada para estimar la probabilidad de que por pura 
## aleatoriedad se obtenga una red similar a la estudiada pero con una longitud media del 
## camino mínimo entre nodos inferior. La función barabasi.game permite generar redes libres 
## de escala con el número de nodos proporcionado en el argumento n.

number.of.added.edges<-10
clustering.coefficients<-NULL
for(i in 1:10000)
{
  if(i%%100==0){
    print(i)
  }
  random.scale.free.graph <- barabasi.game(n=dim(expr)[2],directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

sum(clustering.coefficients > transitivity(gene.coexpression.network,type="global")) / 10000

#Es una red de mundo peque?o

#---Buesqueda de patrones---

## Para la identificación de clusteres o grupos de genes co-expresados en la red necesitamos 
## instalar los siguientes paquetes.

#install.packages("WGCNA")
library("WGCNA")
library("cluster")
allowWGCNAThreads()

## La identificación de clústeres o grupos de elementos se basa en una medida de similitud.
## La medida de similitud seleccionada en este estudio es 1 - correlacion
similarity.matrix <- 1 - gene.correlation

## La función hclust usando la matriz de similitudes como distancias y el método promedio
## para recalcular distancias calcula el clustering jerárquico.
hierarchical.clustering <- hclust(as.dist(similarity.matrix),method="average")

## La función cutree permite cortar el árbol generado en el clustering jerárquico a distintas
## alturas para producir distintos números de clústeres.
hclust.2 <- cutree(hierarchical.clustering,k=2)
hclust.3 <- cutree(hierarchical.clustering,k=3)
hclust.4 <- cutree(hierarchical.clustering,k=4)
hclust.5 <- cutree(hierarchical.clustering,k=5)
hclust.6 <- cutree(hierarchical.clustering,k=6)
hclust.7 <- cutree(hierarchical.clustering,k=7)
hclust.8 <- cutree(hierarchical.clustering,k=8)
hclust.9 <- cutree(hierarchical.clustering,k=9)
hclust.10 <- cutree(hierarchical.clustering,k=10)

## La función pam usa la matriz de similitudes como distancias para determinar clústeres
## según el método de partición entorno a medoides.
pam.2 <- pam(as.dist(similarity.matrix),k=2,diss=TRUE)
pam.3 <- pam(as.dist(similarity.matrix),k=3,diss=TRUE)
pam.4 <- pam(as.dist(similarity.matrix),k=4,diss=TRUE)
pam.5 <- pam(as.dist(similarity.matrix),k=5,diss=TRUE)
pam.6 <- pam(as.dist(similarity.matrix),k=6,diss=TRUE)
pam.7 <- pam(as.dist(similarity.matrix),k=7,diss=TRUE)
pam.8 <- pam(as.dist(similarity.matrix),k=8,diss=TRUE)
pam.9 <- pam(as.dist(similarity.matrix),k=9,diss=TRUE)
pam.10 <- pam(as.dist(similarity.matrix),k=10,diss=TRUE)

## La función silhouette nos permite calcular la silueta de un clustering que sirve de medida
## para la bondad de dicho clustering.
sil2 <- silhouette(hclust.2,dist=similarity.matrix)
sil3 <- silhouette(hclust.3,dist=similarity.matrix)
sil4 <- silhouette(hclust.4,dist=similarity.matrix)
sil5 <- silhouette(hclust.5,dist=similarity.matrix)
sil6 <- silhouette(hclust.6,dist=similarity.matrix)
sil7 <- silhouette(hclust.7,dist=similarity.matrix)
sil8 <- silhouette(hclust.8,dist=similarity.matrix)
sil9 <- silhouette(hclust.9,dist=similarity.matrix)
sil10 <- silhouette(hclust.10,dist=similarity.matrix)

plot(sil3, col = c("#F67E14","#DB5C68","#4CC26C"))

# Graficar el perfil de la silueta con letra más grande
plot(sil3, col = c("#F67E14", "#DB5C68", "#4CC26C"), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)


hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

sil2 <- silhouette(pam.2)
sil3 <- silhouette(pam.3)
sil4 <- silhouette(pam.4)
sil5 <- silhouette(pam.5)
sil6 <- silhouette(pam.6)
sil7 <- silhouette(pam.7)
sil8 <- silhouette(pam.8)
sil9 <- silhouette(pam.9)
sil10 <- silhouette(pam.10)

pam.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]])

## Representamos para los dos métodos de clustering, jerárquico y pam, y para diferentes números
## de clústeres la silueta correspondiente para elegir la mejor combinación de método de 
## clustering y número de clústeres.
plot(2:10,pam.sil.values,type="o",col="#680280",pch=0,ylim=c(0.3,0.8),xlab="Numero de clusters",ylab="Ancho promedio de la silueta",lwd=3)
lines(2:10,hclust.sil.values,type="o",col="royalblue",pch=1,xlab="",ylab="",lwd=3)
legend("topright",legend=c("PAM","HCLUST"),col=c("#680280","royalblue"),pch=c(0,1),lwd=3)

# Assuming you've already loaded ggplot2 and viridis libraries

# Prepare the data
clusters <- 2:10
data <- data.frame(
  Cluster = rep(clusters, 2),
  Silhouette = c(pam.sil.values, hclust.sil.values),
  Method = factor(rep(c("PAM", "HCLUST"), each = length(clusters)))
)

# Ajustamos manualmente los colores
colors <- c("PAM" = "#009C8C", "HCLUST" = "#440154")

ggplot(data, aes(x = Cluster, y = Silhouette, group = Method, color = Method)) +
  geom_line(size = 1.5, aes(color = Method)) +  # Asegura que color = Method para mapear los colores
  geom_point(size = 3, aes(color = Method)) +
  scale_color_manual(values = colors) +  # Usa scale_color_manual para asignar los colores
  labs(x = "Numero de clusters", y = "Ancho promedio de la silueta") +
  theme_classic(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 11),    
    axis.text = element_text(size = 11),     
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20))

## Visualización de clústeres
pam.3[["clustering"]]

clustering.pam.3 <- pam.3[["clustering"]]

#write.table(data.frame(Name = names(clustering.pam.3), `Module` = clustering.pam.3), 
 #           file = "Clusters_3.txt", 
  #          sep = "\t", 
   #         row.names = FALSE, 
    #        col.names = TRUE, 
     #       quote = FALSE)

