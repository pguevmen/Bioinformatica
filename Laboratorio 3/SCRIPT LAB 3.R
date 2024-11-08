# Instalación de librerías si aún no las tienes
install.packages("ape")
install.packages("phangorn")
install.packages("phytools")

# Carga de librerías
library(ape)
library(phangorn)
library(phytools)
update.packages(ask = FALSE, checkBuilt = TRUE)

#para ver en qué directorio estoy
getwd

#ir a una carpeta
setwd("C:/Usuarios/TuUsuario/MiCarpeta") 

#creamos un data frame con las secuencias alineadas y utilizamos la función read.phyDat para convertirlo en un formato llamado phyDat, que es adecuado para análisis filogenéticos.
#el resto es para indicar que el archivo está en formato FASTA y que las secuencias son de aminoácidos
fraxatin <- read.phyDat(file = "clustalo-I20241024-152127-0297-79465327-p1m.aln-fasta", 
                        format = "FASTA", type = "AA")
fraxatin

#Hay 11 secuencias con 310 caracteres y 180 columnas en las que hay diferentes posiciones. 
#Los estados son los aminoácidos A (alanina), R (arginina), etc. 


#Se crea un data frame nuevo en el que cambia el formato del archivo fraxatin a AAbin para poder hacer la matriz de distancia
matrizdist <- as.AAbin(fraxatin)

#Una matriz de distancias es una tabla que muestra entre pares de objetos o elemntos en un conjunto de datos.En biologia, una matriz se utiliza para comparar secuencias de ADN, ARN o proteínas donde la distancia refleja cuán diferentes son dos secuencias entre sí.


#Utiliza dist.aa para hacer la matriz de distancia. Los números corresponden a la distancia evolutiva, es decir, a la diferencia entre las secuencias de la proteína fraxatina de las dos especies. 
#Cuanto menor sea el numero mas parecidas son las secuencias
matrizdist <- dist.aa(matrizdist)
matrizdist

#Ahora creamos un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) usando la matriz de distancia que acabamos de calcular.
#En este tipo de arbol las ramas no indican la distancia evolutiva entre especies
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)

#Ahora hacemos un árbol con el método de unión de vecinos (NJ) usando la misma matriz de distancias:
#En este arbol las ramas sí que indican la distancia evolutiva entre las especies
arbolNJ <- nj(matrizdist)
plot(arbolNJ)

#Para personalizar los árboles podemos utilizar:
#type se utiliza para el tipo de árbol, siendo, p para filograma, c para cladograma, f para un arbol en círculo, u para arbol sin raíces, r para arbol radial
#cex se utiliza para el tamaño de la letra
#edge.width se utiliza para el grosor de las ramas
#edge.color se utiliza para el color de las ramas
#font se utiliza para el estilo de la letra, siendo 1: Fuente regular, 2: Fuente en negrita, 3: Fuente en cursiva, 4: Fuente en negrita y cursiva.
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)

#Además de plot podemos graficar árboles con el método plotTree del paquete phytools, el cual es compatible con ape y con phangorn.
#plotTre permite personalizaciones mucho más detalladas y avanzadas del árbol.
#ftype sirve para el estilo de la letra, siendo, "i": Fuente en cursiva, "b": Fuente en negrita, "r": Fuente regular, "bi": Fuente en negrita y cursiva.
#fsize sirve para el tamaño de la letra
#offset sirve para el espacio entre la rama y el texto
#color sirve para el color de las ramas
#lwd sirve para el grosor de las ramas
#plottree a difernecia de plot del paquete ape proporciona mayor flexibilidad y más detallada para visualizar el árbol
plotTree(arbolNJ)
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)

#En los árboles, sin cambiar la topología, se puede cambiar el orden en que los grupos son visualizados.
#La función ladderize() reorganiza el árbol en "escalera", es decir, coloca los nodos de forma ordenada, ubicando las ramas más largas a la derecho.
plotTree(ladderize(arbolNJ))

#Para guardar un árbol podemos usar el comando write.tree(tree, file = "file_name.nex"). 
#Este archivo puede ser leído usando read.tree(file = file_name.nex).
write.tree(arbolUPGMA, file = "arbol_1")
read.tree(file = "arbol_1")

#Enraizar el árbol es establecer un punto de referencia evolutivo en el árbol colocando un grupo externo (outgroup) como base del árbol. 
#r=TRUE, reorganiza las ramas del árbol para que cuando el grupo externo aparezca en la raíz, se siga manteniendo un orden para mejorar la visualización
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)

#Además podemos visualizar los dos árboles a la vez con los siguientes comandos:
#layout sirve para visualizar dos árboles a la vez
#matrix(c(1,2)) para establecer que los árboles se visualicen en una columna y dos filas (encima uno de otro)
#height=c(10,10) para establecer que la altura de los dos árboles sea la misma
layout(matrix(c(1,2)), height=c(10,10))

#ajusta el espacio de los 4 márgenes
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)

#La función parsimony calcula el grado de parsimonia del árbol. 
#Cuanto mas bajo sea el valor obtenido, menos cambios evolutivos tiene y más parsimonioso es el árbol.
parsimony(arbolUPGMAraiz, fraxatin)
parsimony(arbolUPGMA, fraxatin)

#utilizamos la funcion optim.parsimony para evaluar entre los posibles árboles y obtener el que tenga más parsimonia, es decir, menos cambios evolutivos tenga 
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)

#Otra estrategia para hacer el proceso de búsqueda de árbol con mayor parsimonia es con el algoritmo de búsqueda pratchet. 
#nos da el mejor valor de parsimonia y el número de árboles posibles con ese valor de parsimonia
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
fraxatin_parsimonia

#con lo anterior obtenemos que el mejor valor de parsimonia es 307 y hay 4 árboles posibles con ese valor
#para poder compararlos hay que enraizarlos
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)

#con la función consensus obtenemos un árbol consenso a partir de los 4 árboles anteriores
#p = 1 significa que se está usando un consenso estricto (del 100%), es decir, solo incluyen los grupos monofiléticos o clados que estén presentes en todos los árboles.
#Un clado es un conjunto de organismos que incluye a un ancestro común y a todos sus descendientes
#cex = .6 reduce el tamaño del texto al 60% del valor por defecto.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)

#utilizamos p=0,3 para generar un árbol de consenso del 30% en lugar de un árbol de consenso estricto.
#estamos estableciendo que solo se incluyan en el árbol consenso los clados comunes en al menos el 30% de los árboles
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)

#Genera múltiples árboles a partir de variaciones aleatorias en los datos y evalúa cuán confiables son ciertos clados (grupos de especies relacionadas) en esos árboles.
#FUN indica la función que se va a usar para construir los árboles. 
#pratchet realiza una búsqueda para encontrar el árbol filogenético más probable utilizando el algoritmo de "neighbor-joining"
#se van a generar 10 árboles mediante NJ a partir de 10 réplicas del conjunto de datos.
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
plot(arbolesbootstrap, cex = .6)

#se hace un árbol consenso del 60%, donde solo se muestran los clados que están presentes en al menos el 60% de los árboles de bootstrap.
# consensus crea un árbol consenso a partir de un conjunto de árboles
#arbolesbootstrap es el conjunto de árboles generados por bootstrap 
#P=0,6 indica el umbral de consenso que es del 60%, es decir, sólo los clados que estén presentes en al menos el 60% de los árboles de bootstrap se incluirán en el árbol
# plot se usa para graficar el arbol consenso 
#cex= .6 ajusta el tamaó del texto 
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

#objeto se llama arbolazar
# rtree genera un árbol filogenético aleatorio con un número específico de hojas (especies)
#n=11 indica que el árbol tendrá 11 hojas
#tip label=names permite asignar etiqueas a las hojas del árbol
#plot dibuja el árbol generado
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)


#root enraiza el árbol especificando un grupo externo o outgroup
#phy: indica el arbol que se quiere enraizar (arbolazar)
# outgroup: establece que la raíz del árbol estará en la rama que conecta el resto del árbol con el grupo externo "Ornitorrinco". Esto es útil cuando quieres orientar el árbol y basarlo en un grupo conocido que esté más alejado evolutivamente del resto de las especies.
# ladderize reorganiza el árbol de forma que las ramas más largas queden hacia un lado 
#add.scale.bar() agrega una barra de escala al gráfico, lo cual es útil para representar distancias filogenéticas en unidades proporcionales
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()


#pml() es una función que ajusta un modelo de máxima verosimilitud a un árbol filogenético y datos de seuencia 
#arbolazarR es el árbol que se ajustará 
#fraxatin representa los datos de secuencia o de caracteres asociados a las hojas del árbol. 
ajustado <- pml(arbolazarR, fraxatin)
ajustado

#optimizar un modelo de máxima verosimilitud para un árbol filogenético, aplicando el modelo de sustitución de Dayhoff 
#ptim.pml() optimiza los parámetros de un modelo de máxima verosimilitud previamente definido
#model = "Dayhoff" indica que quieres aplicar el modelo de sustitución de Dayhoff, un modelo comúnmente usado para datos de proteínas, el cual define probabilidades de sustitución específicas para aminoácidos basadas en datos empíricos
#rearrangement = "ratchet" define el método de reordenamiento del árbol durante la optimización. La opción "ratchet" es una estrategia de búsqueda intensiva que aplica perturbaciones aleatorias para salir de óptimos locales y buscar una mejor topología, mejorando la probabilidad de encontrar el mejor árbol posible.
ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

#$tree accede al elemento específico llamado tree dentro del objeto ajustadoconDay, que representa el árbol filogenético optimizado.
ajustadoconDay$tree

#este codigo enraiza el árbol usando el grupo externo y reorganiza el árbol con ladderize de manera que las ramas más largas quedan hacia un lado lo cual facilita la visualización
#root() enraiza el árbol filogenético. Al aplicar esta función al árbol optimizado (ajustadoconDay$tree), le estás añadiendo una raíz especificando un grupo externo o outgroup.
#outgroup = "Ornitorrinco" indica que el grupo externo para enraizar el árbol será "Ornitorrinco". Esto significa que el árbol se orientará de tal forma que el "Ornitorrinco" se ubique en una de las ramas más distales, ayudando a dar contexto evolutivo a las relaciones de los otros taxones.
#ladderize(ajustadoconDayraíz) reorganiza el árbol enraizado. ladderize() es una función que ajusta la disposición de las ramas para que las más largas estén hacia un lado,
#plot(...) grafica el árbol reorganizado, y cex = .5 reduce el tamaño de las etiquetas al 50% para una presentación más compacta.
#add.scale.bar() agrega una barra de escala al gráfico, la cual muestra las distancias evolutivas o el número de sustituciones por unidad de longitud de rama.
ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()


#usa la función optim.pml otra vez pero con un modelo de sustitución diferente y el mismo método de ordenamiento
#optim.pml() ajusta un modelo de máxima verosimilitud a un árbol filogenético dado, tal como ya has hecho antes.
#object = ajustado: Esto especifica que el objeto ajustado, que fue previamente generado con pml(), es el árbol que deseas optimizar.
#model = "Blosum62": Aquí estás cambiando el modelo de sustitución. El modelo Blosum62 es una matriz de sustitución popular para secuencias de proteínas, basada en las frecuencias observadas de sustitución entre aminoácidos. Este modelo es adecuado para análisis filogenéticos con datos de proteínas.
#rearrangement = "ratchet": Al igual que en el ajuste previo con el modelo Dayhoff, estás especificando que se use el método de reordenamiento "ratchet". Este es un método de búsqueda intensiva que intenta mejorar la topología del árbol al realizar pequeñas perturbaciones aleatorias, ayudando a escapar de óptimos locales y buscando una mejor solución.
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")


#realiza la optimización del árbol usando el modelo de sustitución de JTT y el método de ordenamiento de ratchet
#optim.pml(): Esta función optimiza un árbol filogenético bajo un modelo de máxima verosimilitud.
#object = ajustado: Aquí, ajustado es el objeto que contiene el árbol filogenético inicial que deseas optimizar. Este árbol probablemente se haya generado previamente con la función pml.
#model = "JTT": Especificas que deseas usar el modelo JTT (Jones-Taylor-Thornton). Este es un modelo de sustitución de aminoácidos muy utilizado en bioinformática, basado en las frecuencias observadas de sustituciones entre 20 aminoácidos. El modelo JTT es adecuado para datos de proteínas y es uno de los más comunes en análisis filogenéticos.
#rearrangement = "ratchet": El método ratchet es una técnica de reordenamiento que realiza perturbaciones aleatorias para escapar de óptimos locales, lo que permite encontrar una mejor topología del árbol filogenético.
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")


#AIC se usa para calcular el criteio de información de AKaike para modelos estadísticos 
#En este caso, estás usando AIC() para comparar tres árboles optimizados (ajustadoconDay, ajustadoconBlo, y ajustadoconJTT), cada uno ajustado con diferentes modelos de sustitución (Dayhoff, Blosum62, y JTT, respectivamente).
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
#df hace referencia a los grados de libertad del modelo, que generalmente está relacionado con el número de parámetros ajustados en el modelo



#está optimizando el árbol ajustadoconDay utilizando el modelo JTT y el método de reordenamiento ratchet para mejorar la topología
#optim.pml(): Es una función que realiza la optimización de un árbol filogenético utilizando un modelo de sustitución específico y un método de reordenamiento para mejorar la topología y las longitudes de las ramas del árbol.
#object = ajustadoconDay: Aquí, el objeto ajustadoconDay es el árbol que fue optimizado previamente con el modelo Dayhoff. Este es el árbol de entrada para esta nueva optimización.
#model = "JTT": El modelo "JTT" (Jones-Taylor-Thornton) se utiliza para los análisis de proteínas, basado en las frecuencias observadas de sustitución entre los 20 aminoácidos. En este caso, estás cambiando el modelo de sustitución de Dayhoff a JTT.
#rearrangement = "ratchet": El parámetro "ratchet" se refiere a un método de reordenamiento que realiza perturbaciones aleatorias en la topología del árbol para mejorar su ajuste. Este método busca escapar de óptimos locales y encontrar una mejor estructura para el árbol.

mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")

mejorarbol

#está enraizando el árbol optimizando mejorarbol utilizando el outgroup "Ornitorrinco" y luego lo visualizas con una disposición ordenada y una barra de escala.
#root(): Esta función toma un árbol filogenético y lo enraíza utilizando un grupo externo (outgroup). En este caso, el outgroup es "Ornitorrinco", que sirve como referencia para ubicar el punto de enraizamiento del árbol.
#mejorarbol$tree: Este es el árbol optimizado, previamente generado con optim.pml y utilizando el modelo JTT.
#outgroup = "Ornitorrinco": Especifica el taxón "Ornitorrinco" como el grupo externo, que se usará para enraizar el árbol. Este outgroup generalmente se elige porque se considera más lejano evolutivamente a los demás taxones en el árbol.
#ladderize(): Ordena el árbol en una disposición "de escalera" que mejora la legibilidad del árbol al agrupar las ramas largas y cortas en cada nivel.
#mejorarbolRlot(): Muestra el árbol en la ventana gráfica de R.
#cex = 0.5: Ajusta el tamaño del texto en el gráfico, haciéndolo más pequeño. <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()

