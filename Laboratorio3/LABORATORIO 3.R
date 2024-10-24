library(ape)
library(phangorn)
library(phytools)
#preparo para poder meter los algotitmos 
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                        format = "FASTA", type = "AA")
fraxatin
#creamos la matriz
matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist
#creamos un tipo de árbol de distancia tipo upgma
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)
#Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.
#Ahora hacemos un árbol con el método de unión de vecinos (NJ) usando la misma matriz de distancias
arbolNJ <- nj(matrizdist)
plot(arbolNJ)
#Para personalizar los árboles podemos agregar argumentos a parámetros como cex, para el tamaño de la letra, edge.color, para el grosos de las ramas, etc. También se puede escoger entre diferentes visualizaciones de árbol como filograma, cladograma, radial y demás.
# font es negrita(2), cursiva (3), normal (1), cursiva y negrita (4), otro idioma (5)-> letra tipo
#type es el tipo de esquema
#width anchura rama 
#cex tamaño letra
# node.pos mueve la linea
#edge.lty tipo linea continua , dicontinua
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="lightblue", font=2)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=3, node.pos=1, cex=0.8, edge.width=2, edge.color="lightblue", font=3)
#Además de plot podemos graficar árboles con el método plotTree del paquete phytools, el cual es compatible con ape y con phangorn.
plotTree(arbolNJ)
#offset espacios de la letra a las ramas
#lwd anchura rama
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="darkred", lwd=2)
#En los árboles, sin cambiar la topología, se puede cambiar el orden en que los grupos son visualizados. Por ejemplo, se pueden ordenar las puntas de manera alfabética (en la medida de lo posible), o con los grupos más derivados hacia uno de los lados del árbol. Para escalerizar hacia la derecha
plotTree(ladderize(arbolNJ))
#para guardarlo
write.tree(arbolUPGMA, file = "arbolUPGMA.nex")
read.tree(file = "arbolUPGMA.nex")#sale informacion del arbol
#Hasta ahora los árboles que hemos construido no están enraízados. Para enraizarlos podemos usar la función root del paquete ape.
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)#enraizado es diagrama que muestra las relaciones evolutivas entre especies. Todo parte de lo mismo porque es la misma linea pero se va difernciando segun sean más diferentes
#hacemos lo mismo pero con el arbol UPMA
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)
#Además podemos visualizar los dos árboles a la vez con los siguientes comandos
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)
#en parsimonia se evalúan múltiples árboles.La parsimonia no usa todos los caracteres (aminoácidos de la secuecia)
#si comparamos dos árboles nos va a coger el de menos cambios evolutivos y nos los va a contar
parsimony(arbolUPGMAraiz, fraxatin)
#nos da 313 que esto significan los pasos que cambios evolutivos.
parsimony(arbolUPGMA, fraxatin)#aunque esté con raíz o no el número de pasos debe ser el mismo, efectivamente nos da igual
#encontramos arbol con mejor parsimonia
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
#Final p-score 307 after  2 nni operations. Esto quiere decir 
#Intercambio de vecino más cercano (Nearest Neighbor Interchange, NNI): Consiste en intercambir ramas adyacentes hacia una rama interna.
#Podado de subárbol y reinjerto (Subtree pruning and regrafting, SPR): Se corta un pedazo de rama y se inserta en otra rama de manera azarosa.
#Bisección de árbol y reconexión (Tree bisection and reconnection, TBR): Se corta el árbol en la mitad y una porción del árbol que queda de coloca en alguna parte que queda de las ramas del otro árbol.
#queremos lo que tenga menos operaciones en este caso seria el mejorNJ porque nos ha dado lo mismo 307 y encima solo 1 nni operations
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)#te enseña lo que hace para generar la respuesta de 307
fraxatin_parsimonia#nos da el número de árboles que tienen la mejor parsimonia, la de 307
#para compararlos usamos, asi vemos todos 
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)
#Para escoger sólo un árbol con igual parsimonia hay que conseguir un árbol de consenso, qeu tenga todos los grupos que estan presentes en todos.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)
#Para un árbol menos estricto podemos cambiar el valor del parámetro p, quiere decir qeu nos ataja nos saca todo de uno, es menos especifico , menos ramas
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)
#generamos varias pseudoreplicas, esto es una nueva version de la matriz, donde se pueden obtener matrices con caracteres repetidos o alguno que no se use.
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)#Nos da valores en los que nos puntua los arboles con menos cambios como los mejores 
plot(arbolesbootstrap, cex = .6)#La rutina anterior genera 10 árboles pseudoréplicas.
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)#como hemos hecho antes, uno donde haya muchos en común
#Se calcula la verosimilitud (probabilidad de obtener un dato según un modelo) de un árbol de acuerdo a un alineamiento de secuencias usando un modelo de sustitución de aminoácidos. También se tienen en cuenta la frecuencia de de ocurrencia de los aminoácidos.
#Se calcula la verosimilitud (probabilidad de obtener un dato según un modelo) de un árbol de acuerdo a un alineamiento de secuencias usando un modelo de sustitución de aminoácidos. También se tienen en cuenta la frecuencia de de ocurrencia de los aminoácidos.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)# rtree nos genera un arbol aleatorio, en este caso me crea debido a n=11 uin arbol aleatorio con 11 especies.Con tiplabel metemos las especies qwe van a aparecer.
#lo enraizamos por las secuencias de Ornitorinco.Además los «escalerizamos» hacia la derecha y le agregamos escala; aquí la longitud de la rama sí es significativa, indica cantidad de cambio en cuanto a sustituciones de aminoácidos.
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()
#iniciamos la búsqueda del mejor árbol por máxima verosimilitud. Lo primero que se hace es calcular la verosimilitud del árbol dadas las secuencias.
ajustado <- pml(arbolazarR, fraxatin)
ajustado
# nos da:verosimilitud del árbol al azar que habíamos creado, que es -4348.064.-> loglikelihood: -4450.32 
#Lo que hay que hacer es encontrar un árbol que optimice la verosimilitud usando un modelo de sustitución; para esto vamos a usar el método optim.pml del paquete phangorn, el cual computa la verosimilitud de un árbol filogenético dado un alineamiento múltiple de secuencias y un modelo de evolución de AA. Toma como argumentos un objeto de clase pml, el tipo de modelo que se quiere usar así como el tiempo de rearreglo para los árboles.
ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")
#Para ver el árbol oculto usamos $tree. También lo enraizamos.
ajustadoconDay$tree
ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")#  hemos usado Dayhoff
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()
#otro método
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
#otro
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")
#Podemos comparar los modelos calculando el Criterio de información de Akaike AIC
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")
mejorarbol
mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()