# GWAS y Selección Genómica  

## **Autor**  
### Carlos Gutierrez Ferreira  
- Chileno  
- Magíster en Biotecnología

## Conexión y acceso a RStudio Cloud

- Se accedió a la página [RStudio.Cloud](https://rstudio.cloud/projects) para generar el [Proyecto](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/GWAS_Seleccion_Genomica/Archivos/Practica_GWAS_asociacion_genomica.Rmd) en formato RMarkdown.

## instalación de las librerías y configuración de los trozos de código.

- Se colocaron las siguientes líneas de código para configurar la sesión.

```{r, echo=FALSE, message = FALSE, warning = FALSE}
knitr::opts_chunk$set(message = FALSE)
install.packages("rrBLUP")
install.packages("ggplot2")
```
```{r, echo=FALSE}
library(utils)
library(rrBLUP)
library(ggplot2)
```

## Importado y exportado de los archivos correspondientes a genotipos y fenotipos.

### a) Se importaron los archivos [*geno*](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/GWAS_Seleccion_Genomica/Archivos/geno.txt) y [*pheno*](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/GWAS_Seleccion_Genomica/Archivos/pheno.txt) al entorno de RStudio.

```{r}
geno<-read.delim("geno.txt",sep="\t",dec=",",header=T)
pheno<-read.delim("pheno.txt",sep="\t",dec=",",header=T)
```

### b) Se exportaron y analizaron los datos.

```{r}
dim(geno)
dim(pheno)
head(geno[1:9,1:9])
head(pheno)
```

![1](https://user-images.githubusercontent.com/80927233/125408714-2d99a200-e389-11eb-8d0b-ea9699da6b2f.jpg)

```{r}
hist(pheno$y,main="Histograma del fenotipo",xlab="Fenotipo observado")
```

![2](https://user-images.githubusercontent.com/80927233/125408737-35f1dd00-e389-11eb-98fc-58ce24c16fa6.jpg)

- Cada SNP se encuentra codificado como presente (1) o ausente (0) en el genotipo de cada animal.  
- Los heterocigotos son todos aquellos individuos que se encuentran al rededor del valor 0 en el Histograma del fenotipo, dado que no presentan rasgos homocigotos de ambos alelos.

## Matriz de relación aditiva (A.mat {rrBLUP}) y análisis GWAS.

### a) Se investigó la función *A.mat* y se calculó la matriz de parentesco genómico para el set de datos *geno*.

```{r}
help(A.mat)
A<-A.mat(geno[4:203]) 
dim(A)
head(A[1:6,1:6])
```
```{r}
hist(A,main="Histograma de la matriz de relación aditiva")
```
![3](https://user-images.githubusercontent.com/80927233/125408926-676aa880-e389-11eb-8325-edf79b348a2d.jpg)

### b) Grafica de la diagonal de la matriz.

```{r}
endogamia<-diag(A)
mean(endogamia)
hist(endogamia,main="Histograma de endogamia")
```

![4](https://user-images.githubusercontent.com/80927233/125409046-836e4a00-e389-11eb-9a8d-2104a64057d3.jpg)

- El nivel de endogamia promedio de esta población es de 1.00026.  
- Un valor de endogamia sobre 1 (1.1) representa a un coeficiente de consanguinidad positivo en los individuos analizados, con una correlación gamética intraindividual positiva en locus único, utilizado cuando la densidad de marcadores es baja.  
- Asi mismo, un valor de endogamia bajo 1 (0.9) representa a un coeficiente de consanguinidad negativo.

### c) Se investigó la función *GWAS* y se realizó un análisis de asociación genómica GWAS.

```{r}
help(GWAS)
score<-GWAS(pheno,geno,plot=TRUE)
class(score)
```

![5](https://user-images.githubusercontent.com/80927233/125409106-9254fc80-e389-11eb-95f2-7f5eb398b2da.jpg)

- Al observar la grafica generada por la función *GWAS* se detectaron dos QTLs significativos en el análisis, uno en el cromosoma 3 y otro en el cromosoma 10.

## Efecto de los QTLs detectados en el análisis GWAS.

### a) Se exploró el objeto *score* con los comandos *head* y *View*

```{r}
head(score)
View(score)
dplyr::filter(score,y>5)
exp(-8.508100)
exp(-7.5047236)
```

![6](https://user-images.githubusercontent.com/80927233/125409166-a0a31880-e389-11eb-983a-cefadddc0bd1.jpg)

- Tal como se observo en el análisis GWAS, solo se encontraron dos SNPs significativos, correspondientes al snp300 y snp1000 presentes en los cromosomas 3 y 10 respectivamente.  
- Entre ambos, el snp300 presenta el valor de p más alto, con un valor de 5.5x10^-4  
- Cabe destacar que el analisis GWAS considera por defecto la significancia de los marcadores con un FDR de 0.05 calculado mediante el q-value.

### b) Regresión lineal del fenotipo en función del genotipo para cada SNP significativo.

```{r}
qtl300<-t(geno[300,4:203])+1
qtl1000<-t(geno[1000,4:203])+1
qtls<-data.frame(qtl300,qtl1000,pheno$y)
head(qtls)
```

![7](https://user-images.githubusercontent.com/80927233/125409217-b1ec2500-e389-11eb-8951-f6ddda2f637f.jpg)

```{r}
qtl1<-ggplot(qtls,aes(x=X300,y=pheno.y))
qtl1+geom_point()+xlab("snp 300")+ylab("Pheno")+geom_smooth(method=lm)
```

![8](https://user-images.githubusercontent.com/80927233/125409234-b57fac00-e389-11eb-9ef0-d6336f9b5184.jpg)

```{r}
qtl2<-ggplot(qtls,aes(x=X1000,y=pheno.y))
qtl2+geom_point()+xlab("snp 1000")+ylab("Pheno")+geom_smooth(method=lm)
```

![9](https://user-images.githubusercontent.com/80927233/125409249-b9133300-e389-11eb-9e6c-95dfeb4f7b62.jpg)

### c) Estimacion del efecto (beta o pendiente) de los QTLs con mayor score.

```{r}
lm.qtl.300<-lm(pheno.y~X300,data=qtls)
summary(lm.qtl.300)
```

- El efecto del snp300 sobre el rasgo cuantitativo es de 1.8121 respecto al eje x.

```{r}
lm.qtl.1000<-lm(pheno.y~X1000,data=qtls)
summary(lm.qtl.1000)
```

- El efecto del snp1000 sobre el rasgo cuantitativo es de 1.8549 respecto al eje x.

![10](https://user-images.githubusercontent.com/80927233/125409286-c03a4100-e389-11eb-8d95-f2187c9917ec.jpg)

