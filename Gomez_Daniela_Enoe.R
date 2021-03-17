################ Primer parcial de Genomica
##### P2 Red de IPP
library(igraphdata)
library(igraph)
data("yeast")
class(yeast)


# 10 proteinas mas conectadas 
proteinas_mas_conectadas<-sort(degree(yeast),decreasing = TRUE)[1:10]
names(proteinas_mas_conectadas)
# Grafica de conectividades
hist(degree(yeast))
# Diametro de la red
diameter(yeast)
# Coeficiente de clusterizacion de las 20 proteinas mas conectadas 
p20_mas_conectadas<-sort(degree(yeast),decreasing = TRUE)[1:20]
names(p20_mas_conectadas)

transitivity(yeast, vids = c("YPR110C", "YPL131W", "YNL178W", "YIL021W","YOL127W","YJL063C","YBR283C","YLR378C","YMR260C","YBR251W","YGL123W",
                             "YNL284C","YGL103W","YOL040C","YGR034W","YDL191W","YDL136W","YGR220C","YNR037C","YOR063W") ,type="local")

# Encontrar proteinas de coeficiente de clusterizacion de 1
proteinas_ci_1<-which(transitivity(yeast, type= "local")== 1)
proteinas_ci_1 # Te dice las posiciones de donde se encuentran las proteinas del objeto yeast

# Interpretación : Si tenemos un coeficiente de clusterizacion de 1 en un nodo o vertice, es que esta completamente conectada con sus vecinos,
# donde la escala de valores de 0(no esta conectado con sus vecinos) a 1, donde muestra el porcentaje de que tan conectados esta el nodo con sus vecinos.

###El porcentaje de conexiones respecto al total
# Primero hay que calcular :El degree cuando una red esta completamente conectada
n_nodos<-length(V(yeast))
n_nodos
links_red_comp_conectada<-2617*(2617-1)/2
links_red_comp_conectada
n_links<-length(E(yeast))
n_links
porcentaje_conexiones<-(n_links/links_red_comp_conectada)*100
porcentaje_conexiones

# El promedio de conectividades
mean(degree(yeast))

# las 10 proteínas más importantes con  3 metodos de centralidad 
## Excentricidad
eccentricity<-sort(eccentricity(yeast),decreasing = TRUE)[1:10] 
names(eccentricity)
## Closeness     
Closeness<-sort(closeness(yeast),decreasing = TRUE)[1:10] 
names(eccentricity)
## Betweennes
Betweennes<-sort(betweenness(yeast),decreasing = TRUE)[1:10] 
names(Betweennes)

## Encuentre los caminos entre las proteínas más alejadas 
diameter(yeast)

## Seleccione 1 proteína al azar y la la elimine de la red, 
## que  calcule el promedio de las distancias después de quitar la proteína. Hacer esto 100 veces más y compare los promedios
for(i in 0:99){
  yeast_proteinas_nueva<-delete.vertices(yeast, sample(1:(100-i),1))
  print(mean_distance(yeast_proteinas_nueva))
}

#### Discusion: Una red biologica esta construida asi porque cada vez que quita un nodo al azar, afecta o no la arquitectura de la
# red si el nodo que se removio de la red poseia muchas conexiones o era un hub, entonces si quitas a los hubs de la red, 
# en una red biologica es distribucion free scale, entonces es vulnerable si atacas los hubs. 



##### P1 Insulinas fasta
# Librerias
library(Biostrings)
library(msa)
library(seqinr)
#  Leer el archivo
sec_proteinas<-readAAStringSet("C:/Users/Daniela/Downloads/Insulinas.fasta")
sec_proteinas
#2 alineamientos múltiples por dos métodos distintos
ClustalW<-msa(sec_proteinas)
ClustalW
Muscle<-msa(sec_proteinas,method = "Muscle")
Muscle

# Construir una matriz de distancia por cada método de alineamiento
secuencias_Alineadas<-msa(sec_proteinas)
secuencias_Alineadas_2<-msaConvert(secuencias_Alineadas, type = "seqinr::alignment")# Convertir este alineamiento a un objeto para otros analisis
secuencias_Alineadas_2

secuencias_Al<-msa(sec_proteinas, method= "Muscle")
secuencias_Al_2<-msaConvert(secuencias_Al, type = "seqinr::alignment")# Convertir este alineamiento a un objeto para otros analisis
secuencias_Al_2
# Previamente instalar y cargar la libreria seqinr
d1<-dist.alignment(secuencias_Alineadas_2,"identity")
d2<-dist.alignment(secuencias_Alineadas_2,"identity")
# Objeto para calcular la distancia del alineamiento, objeto anterior, identity= matriz de distancia por coincidencia, esta sea por identidad
d1<-as.matrix(d1) # Lea como una matriz
d1

d2<-as.matrix(d2)
d2
#Construya dos redes pesadas a partir de cada matriz de distancia
red_pesasa1<-graph_from_adjacency_matrix(d1, mode="undirected", weight=T)
red_pesada1

red_pesasa2<-graph_from_adjacency_matrix(d2, mode="undirected", weight=T)
red_pesada2

# Grafique las redes tomando en cuenta su peso
plot(d1,edge.width-E(d1)$weight*2.5, edge.color="pink")
plot(d2,edge.width-E(d2)$weight*2.5, edge.color="orange")

# Clusterizar la red con al menos tres métodos de agrupamiento
propagacion_1<-label.propagation.community(d1)
propagacion_2<-label.propagation.community(d2)

edge_betweenness_1<-cluster_edge_betweenness(d1)
edge_betweenness_2<-cluster_edge_betweenness(d2)

cluster_spinglass_1<-cluster_spinglass(d1)
cluster_spinglass_2<-cluster_spinglass(d2)



