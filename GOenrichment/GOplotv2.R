######################################################
## Tema: Gráficas de burbuja para enriquecimiento   ##
## Autor: Olga Andrea Hernandez Miranda, Miranda H  ##
## Fecha: 10/11/2023                                ##
## Nota: Enriquesimiento funcional en con GO        ##
######################################################

library(ggplot2)
library(forcats)
library(tidyverse) 
library(dplyr)
library(viridis)
library(ggthemes)

setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutoral3/Enriquesimiento/EnriquesimientoHM")


mydat <- read.table("FecBubbleActualizada.csv", sep = ",", header = T)

head(mydat)


mydat <- mydat[order(mydat$Value, decreasing = TRUE), ]

primeras_20_filas <- head(mydat, 10)

head(primeras_20_filas)
# Suponiendo que tu tabla se llama "mi_tabla"
#mi_tabla_filtrada <- mydat[mydat$Value > 50000, ]


#Con Log10PValue
 
ggplot(primeras_20_filas, aes(y = GO_term, x = LogSize, size = LogSize, color = Value)) +
  geom_point()

ggplot(primeras_20_filas, aes(x= LogSize,y= GO_term)) +
    geom_point(aes(color= Value,size= LogSize)) +
    scale_color_gradientn(colours = rainbow(5)) +
                            ggtitle("Biologicarl process")


# Reorder the y-axis (GO_term) based on the negative of Log10Pvalue
primeras_20_filas$GO_term <- reorder(primeras_20_filas$GO_term, -primeras_20_filas$Value)

# Reverse the order of levels to order from largest to shortest Log10Pvalue
primeras_20_filas$GO_term <- factor(primeras_20_filas$GO_term, levels = rev(levels(primeras_20_filas$GO_term)))

#Create the bubble plot with the reordered y-axis
ggplot(primeras_20_filas, aes(x = LogSize, y = GO_term)) +
  geom_point(aes(color = Value, size = LogSize)) +
  scale_color_gradientn(colours = rainbow(5)) +
  ggtitle("Biological process")

# Create the bubble plot with the viridis color palette
ggplot(primeras_20_filas, aes(x = LogSize, y = GO_term)) +
  geom_point(aes(color = Value, size = LogSize)) +
  scale_color_viridis(option = "viridis") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14))  # Aumentar el tamaño del texto en el eje Y




