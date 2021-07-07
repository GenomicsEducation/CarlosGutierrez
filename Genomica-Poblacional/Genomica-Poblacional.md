# Introducción a la genómica poblacional y ancestría

## **Autor**
### Carlos Gutierrez Ferreira  
- Chileno
- Magíster en Biotecnología

## Conexión al servidor Pomeo  

![PUTTY3](https://user-images.githubusercontent.com/80927233/119919416-67792b00-bf38-11eb-8e85-ffe2a8c69777.jpg)

## Configuración de bioconda e instalación de software  

```
# Se configuró el canal de Bioconda para su uso.
conda config --add channels bioconda
# Se instaló el software plink.
conda install -c bioconda plink
# Se instaló el software admixture.
conda install -c bioconda admixture
```

## Creación del directorio de trabajo “population” y preparación de los archivos para el análisis poblacional.  

[EU_OC_US.vcf](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/EU_OC_US.vcf)  
[Admixture_plot.R](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Admixture_plot.R)  

```
# Se creó la carpeta "population" para agregar los archivos "EU_OC_US.vcf" y "Admixture_plot.R".
mkdir population
cd population
# Se copiaron los archivos "EU_OC_US.vcf" y "Admixture_plot.R" mediante WinSCP y se revisó el contenido de la carpeta.
ls -l -h

```

- ### El archivo [EU_OC_US.vcf](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/EU_OC_US.vcf) contiene las muestras provenientes de tres poblaciones de Salmon del Atlántico (Salmo salar).  

- Europa: 2_WG0341511-DNA_A02_5408, 3_WG0341511-DNA_A03_5416, 5_WG0341511-DNA_A05_5450.  

- Oceanía: FR07958834, FR07958842, FR07958850.  

- Norteamérica: GNB12-1, GNB12-10, GNB12-11.  

## Análisis de diversidad  

### Se estimó la cantidad de sitios heterocigotos para cada individuo y la heterocigosidad observada y esperada para cada marcador

```
vcftools --vcf EU_OC_US.vcf --het --out EU_OC_US
vcftools --vcf EU_OC_US.vcf --hardy --out EU_OC_US

```  

- Analizando los archivos de salida con "less" se pudo observar que estos contienen el cromosoma y la posición del gen, la cantidad de individuos homocigotos, heterocigotos y la proporción de esta distribución.  

### Se calculó la diversidad en una ventana no superpuesta de 200 kb para cada individuo de las tres poblaciones

```
# Para la poblacion de Europa.
vcftools --vcf EU_OC_US.vcf --window-pi 200000 --indv 2_WG0341511-DNA_A02_5408 --indv 3_WG0341511-DNA_A03_5416 --indv 5_WG0341511-DNA_A05_5450 --out EU
# Para la poblacion de Oceania.
vcftools --vcf EU_OC_US.vcf --window-pi 200000 --indv FR07958834 --indv FR07958842 --indv FR07958850 --out OC
# Para la poblacion de Norteamerica.
vcftools --vcf EU_OC_US.vcf --window-pi 200000 --indv GNB12-1 --indv GNB12-10 --indv GNB12-11 --out US

```

### Se Calculó el desequilibrio de ligamiento (LD) para las tres poblaciones

```
# Para la poblacion de Europa.
vcftools --vcf EU_OC_US.vcf --geno-r2 --chr 1 --ld-window-bp 100000 --min-r2 0.001 --indv 2_WG0341511-DNA_A02_5408 --indv 3_WG0341511-DNA_A03_5416 --indv 5_WG0341511-DNA_A05_5450 --out EU
# Para la poblacion de Oceania.
vcftools --vcf EU_OC_US.vcf --geno-r2 --chr 1 --ld-window-bp 100000 --min-r2 0.001 --indv FR07958834 --indv FR07958842 --indv FR07958850 --out OC
# Para la poblacion de Norteamerica.
vcftools --vcf EU_OC_US.vcf --geno-r2 --chr 1 --ld-window-bp 100000 --min-r2 0.001 --indv GNB12-1 --indv GNB12-10 --indv GNB12-11 --out US

```
- Mientras más cercano a 1 es el valor R^2 del desequilibrio de ligamiento, menores son las probabilidades de que ambos genes, evaluados mediante LD, segreguen de forma independiente.  
- Cada uno de estos comandos filtró parte de los datos según las ventanas indicadas, con el objetivo de estandarizar las muestras poblacionales, obteniendo a su vez, un [texto de salida](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Output.txt) después de cada comando.  

## Análisis de los datos en RStudio.

- Se accedió a la página de [RStudio.Cloud](https://rstudio.cloud/projects) para realizar los análisis gráficos de los datos previamente obtenidos.  

### Carga de los archivos al proyecto en RStudio.Cloud e instalación de las librerías.  

- Se ingresaron todos los archivos de la carpeta "[population.zip](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/population.zip)" en formato zip, al nuevo proyecto RStudio.  
- Se cargaron las siguientes librerías utilizando el [Sript](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Genomica.Rmd) en formato R-Markdown.

```{r setup}

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(cowplot)
```
- La librería "tidyverse" demora en su instalación debido a la cantidad de comandos que presenta en su repertorio.

### Promedio de los datos y corte de la muestra según las ventanas indicadas

```{r}

EU <- read_delim("EU.geno.ld", delim = "\t")
OC <- read_delim("OC.geno.ld", delim = "\t")
US <- read_delim("US.geno.ld", delim = "\t")

EU$dist <- ceiling((EU$POS2 - EU$POS1)/1000)*1000
OC$dist <- ceiling((OC$POS2 - OC$POS1)/1000)*1000
US$dist <- ceiling((US$POS2 - US$POS1)/1000)*1000

EU2 <- group_by(EU,dist) %>% summarise(meanR2 = mean(`R^2`))
OC2 <- group_by(OC,dist) %>% summarise(meanR2 = mean(`R^2`))
US2 <- group_by(US,dist) %>% summarise(meanR2 = mean(`R^2`))  

dd <- bind_rows(EU2,OC2,US2)
dd$pop <- rep(c("EU","OC","US"),c(nrow(EU2),nrow(OC2),nrow(US2))) 
write_csv(dd,"EU_OC_US.windowed.ld.csv")

```

### Heterocigosidad individual.  

```{r}

het <- read_delim("EU_OC_US.het",delim = "\t")
het

het$Heterozygosity <- 1-(het$`O(HOM)`/het$N_SITES) 
het$Population <- c(rep("EU",3),rep("OC",3),rep("US",3))
Aplot <- ggplot(het,aes(x = Population, y = Heterozygosity, col = Population)) +
  geom_point()+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")
Aplot

```

### Diversidad de nucleótidos.  

```{r}

pi_EU <- read_delim("EU.windowed.pi",delim = "\t")
pi_EU
pi_OC <- read_delim("OC.windowed.pi",delim = "\t")
pi_OC
pi_US <- read_delim("US.windowed.pi",delim = "\t")
pi_US

pi_all <- bind_rows(pi_EU,pi_OC,pi_US)
pi_all$Population<-c(rep("EU",nrow(pi_EU)),rep("OC",nrow(pi_OC)),rep("US",nrow(pi_US)))

Bplot <- ggplot(pi_all,aes(x = Population, y = PI, col = Population))+
      geom_jitter(col = "grey",width = 0.1)+ 
      geom_boxplot(notch = T, alpha = 0,outlier.shape = NA)+ 
      theme_bw()+
      theme(legend.position = "none")+
      xlab("")+
      ylab(expression(pi))
Bplot

```

### Desequilibrio de ligamiento.  

```{r}

ld <- read_csv("EU_OC_US.windowed.ld.csv")
ld

Cplot <- ggplot(ld,aes(x = dist/1000, y = meanR2, col = pop)) +
      geom_point()+
      geom_line()+
      theme_bw()+
      xlab("Distance (kb)")+
      ylab(expression(R^2))+
      scale_colour_discrete(name = "Population")
Cplot

```

### Gráfico de paneles múltiples.  

```{r}

top_row <- plot_grid(Aplot,Bplot,labels = "AUTO")
plot_grid(top_row,Cplot,nrow = 2,labels = c("","C"))

```

![Genomica-Poblacional](https://user-images.githubusercontent.com/80927233/124712104-e153eb00-decc-11eb-89b2-68af763ce0db.jpg)

- Grafico 1: Análisis de los datos genómicos provenientes de las muestras poblacionales de Europa (EU), Oceanía (OC) y EE.UU. (US). A) Heterocigocidad esperada y promedio, B) Diversidad de nucleótidos y C) Desequilibrio de ligamiento (LD) en las tres poblaciones.

## Análisis de estructura poblacional.  

- Volviendo al servidor Pomeo en PuTTY (Windows) se ejecutaron los siguientes comandos para extraer los datos filtrados y podados, obteniendo un [texto de salida](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Output2.txt) después de cada comando.

```
# Se generó el archivo de entrada en formato plink.
plink --vcf EU_OC_US.vcf --recode --out EU_OC_US --double-id --allow-extra-chr --chr-set 29
Se generó el archivo de entrada en formato plink binario.
plink --file EU_OC_US --make-bed --out EU_OC_US --allow-extra-chr --chr-set 29
# Filtrado basado en equilibrio de Hardy-Weinberg y frecuencia del alelo menor.
plink --bfile EU_OC_US --hwe 0.01 --maf 0.05 --make-bed --out EU_OC_US.Filtered --allow-extra-chr --chr-set 29
# Filtrado y exclusión de marcadores por desequilibrio de ligamiento.
plink --bfile EU_OC_US.Filtered --indep-pairwise 50 10 0.05 --make-bed --out EU_OC_US.Filtered --allow-extra-chr --chr-set 29
plink --bfile EU_OC_US.Filtered --extract EU_OC_US.Filtered.prune.in --make-bed --out EU_OC_US.FilteredPruned --allow-extra-chr --chr-set 29
# Filtrado para remover individuos relacionados.
plink --bfile EU_OC_US.FilteredPruned --rel-cutoff 0.4 --out EU_OC_US.FilteredPruned --allow-extra-chr --chr-set 29
plink --bfile EU_OC_US.FilteredPruned --keep EU_OC_US.FilteredPruned.rel.id --make-bed --out EU_OC_US.FilteredPrunedUnrel --allow-extra-chr --chr-set 29
# Preparación de los archivos para PCA (Principal Component Analysis). 
plink --bfile EU_OC_US.FilteredPrunedUnrel --pca 4 --out EU_OC_US.FilteredPrunedUnrel --allow-extra-chr --chr-set 29

```

### Graficos de PCA con R

- Volviendo a [RStudio.Cloud](https://rstudio.cloud/projects) se ejecutaron los siguientes comandos presentes en [Sript](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Genomica.Rmd) indicado.

```{r}

pca1 <- read_delim("EU_OC_US.FilteredPrunedUnrel.eigenvec", delim = " ",col_names = F)

colnames(pca1) <- c("Population","Individual",paste0("PC",c(1:4)))

mycols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")

Dplot <- ggplot(pca1,aes(x = PC1,y = PC2,col = Population))+
      geom_point()+
      theme_bw()+
      scale_colour_manual(values = mycols)
Dplot

```

![PCA](https://user-images.githubusercontent.com/80927233/124719087-c5ecde00-ded4-11eb-9731-c76a185110ac.jpg)

- Gráfico 2: Análisis de los componentes principales (PCA) de los datos filtrados y podados desde la población.  

## Análisis mediante admixture (PuTTY).

```
# Se seleccionaron al azar el 1% de los marcadores presentes en la muestra filtrada y podada de los individuos no relacionados.
plink --bfile EU_OC_US.FilteredPrunedUnrel --thin 0.01 --make-bed --out EU_OC_US.Thinned --allow-extra-chr --chr-set 29
# Análisis de ancestría de 2 a 6 poblaciones.

for K in `seq 2 6`;
do
admixture EU_OC_US.Thinned.bed $K;
done

```

- Admixture generó 2 archivos: Los archivos .Q que contienen asignaciones de grupos para cada individuo y los archivos .P que contienen las frecuencias alélicas de la población para cada SNP.
- [Texto de salida](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Archivos/Output3.txt) obtenido posterior a los comandos.

### Graficos generados por Admixture para 2, 4 y 6 poblaciones (RStudio.Cloud).

```{r}

source("Admixture_plot.R")

pops <- read_delim("EU_OC_US.Thinned.fam", delim = " ",col_names = F)

K2 <- read_delim("EU_OC_US.Thinned.2.Q", delim = " ",col_names = F)
Eplot <- admixtureplot(str_out = K2,k = 2, pops = pops, xaxis = F)
Eplot

K4 <- read_delim("EU_OC_US.Thinned.4.Q", delim = " ", col_names = F)
Gplot <- admixtureplot(str_out = K4,k = 4, pops = pops, xaxis = F)
Gplot

K6 <- read_delim("EU_OC_US.Thinned.6.Q", delim = " ", col_names = F)
Hplot <- admixtureplot(str_out = K6,k = 6, pops = pops, xaxis = T)
Hplot

top_row <- plot_grid(Eplot,Gplot,labels = "AUTO")
plot_grid(top_row,Hplot,nrow = 2,labels = c("","C"))

```

![Admixture](https://user-images.githubusercontent.com/80927233/124721754-62b07b00-ded7-11eb-81e3-277dfeebd404.jpg)

- Gráfico 3: Simulación de cruza genética entre 2, 4 y 6 poblaciones (A, B y C respectivamente) observando la heredabilidad de los marcadores correspondientes a cada una.

