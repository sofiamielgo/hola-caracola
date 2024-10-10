#¿Cuántos registros cumplen las condiciones finales?
  #En las condiciones de sólo humanos: 35
  #En las condiciones de sólo humanos del planeta Tattoine: 87
	#En las condiciones de todas las especies menos los Droides: 77
#¿Cómo calcularías la desviación estándar (sd) de esos parámetros?
  #Primero ponemos nombre a nuestros parametros:
  coco=starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))
  #Depués eliminamos la columna de nombres:
    fruta=coco%>%select(contains("_"))
  #Luego creamos un vector de la línea de height:
    mandarina=fruta%>%pull(mean_height)
  #Finalmente conseguimos la desviación para height: que nos da: [1] 42.24324
  sd(mandarina) 
  #Hacemos lo mismo para mass, pero utilizamos cuando queramso saber la desviación: [1] 229.3549
  sd(pera, na.rm=TRUE) #Porque hay un valor que pone NA y no lo queremos. 
  #Al crear los gráficos puedes observar que hay un punto que corresponde a un personaje con una masa muy grande. Inspecciona el datset, filtra usando las funciones de tidyverse, crea un nuevo dataframe sin ese personaje y crea de nuevo el gráfico final. (Exporta el gráfico con la opción exportar en el panel derecho y adjúntalo en el pdf)
  #Para no aparezca he eliminado la fila utilizando el comando:
    starwars_filtered <- starwars[-16, ]
  #Inspecciona el dataset, haz un resumen de la media (mean) de las variables (Peso, Altura,IMC, IAS, CCintura). Agrupando por sexo.
    Mujeres<-toy%>%filter(Sex=="Women")
    Hombres<-toy%>%filter(Sex=="Men")toy
    toy %>% group_by(Sex) %>% summarise(mean_peso = mean(Weight_Kg, na.rm = TRUE),mean_altura = mean(Height_cm, na.rm = TRUE),mean_imc = mean(IMC, na.rm = TRUE),mean_ias = mean(IAS, na.rm = TRUE),mean_ccintura = mean(Ccintura, na.rm = TRUE))
  #Haz una tabla sólo con los pacientes femeninos ¿Cuántos registros cumplen las condiciones? ¿De estos cuantos tienen Sobrepeso (Overweight)?
  femenino <- toy %>% filter(Sex == "Women")
  toy %>% filter(Sex == "Women") %>% tally()
  #58
  femenino %>% filter(IMC_clas == "Overweight") %>% tally()
  #9
  
  #Haz un gráfico usando ggplot relacionando el IMC (Indice de masa corporal) con el peso (Weight_Kg) de todos los pacientes.
  ggplot(toy, aes(IMC, Weight_Kg)) + geom_point(colour = "red") + theme_light()
  	#Repítelo filtrando sólo los pacientes categorizados como "Overweight" y "Obesity".
obesidad <- toy %>% filter(IMC_clas!="Normal")
ggplot(obesidad, aes(IMC, Weight_Kg)) + geom_point(colour = "red") + theme_light()
