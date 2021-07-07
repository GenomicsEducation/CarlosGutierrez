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
# Se instaló el software admixture
conda install -c bioconda admixture
```

## Creación del directorio de trabajo “population” y preparación de los archivos para el análisis poblacional con plink y admixture  

[EU_OC_US.vcf](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/EU_OC_US.vcf)  
[Admixture_plot.R](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/Admixture_plot.R)  

```
# Se creó la carpeta "population" para agregar los archivos "EU_OC_US.vcf" y "Admixture_plot.R".
mkdir population
cd population
# Se copiaron los archivos "EU_OC_US.vcf" y "Admixture_plot.R" mediante WinSCP y se reviso el contenido de la carpeta.
ls -l -h

```

- ### El archivo [EU_OC_US.vcf](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Genomica-Poblacional/EU_OC_US.vcf) contiene las muestras provenientes de tres poblaciones de salmon del Atlantico (Salmo salar).

- Europa: 2_WG0341511-DNA_A02_5408, 3_WG0341511-DNA_A03_5416, 5_WG0341511-DNA_A05_5450.

- Oceania: FR07958834, FR07958842, FR07958850.

- Norteamerica: GNB12-1, GNB12-10, GNB12-11.

```{r setup}

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(cowplot)
```

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

### Heterocigosidad individual

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

```{r}

pi_EU <- read_delim("EU.windowed.pi",delim = "\t")
pi_EU

pi_OC <- read_delim("OC.windowed.pi",delim = "\t")
pi_OC

pi_US <- read_delim("US.windowed.pi",delim = "\t")
pi_US

```

```{r}

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

```{r}

top_row <- plot_grid(Aplot,Bplot,labels = "AUTO")
plot_grid(top_row,Cplot,nrow = 2,labels = c("","C"))

```



```{r}

pca1 <- read_delim("EU_OC_US.FilteredPrunedUnrel.eigenvec", delim = " ",col_names = F) %>% head(pca1)

colnames(pca1) <- c("Population","Individual",paste0("PC",c(1:4))) %>% head(pca1)

mycols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")

Dplot <- ggplot(pca1,aes(x = PC1,y = PC2,col = Population))+
      geom_point()+
      theme_bw()+
      scale_colour_manual(values = mycols)
Dplot

```

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

```

