# Introducción al análisis de secuencias NGS - Análisis de control de calidad, filtrado y poda  

## **Autor**
### Carlos Gutierrez Ferreira  
- Chileno
- Magíster en Biotecnología

## Conexión al servidor Pomeo

![PUTTY3](https://user-images.githubusercontent.com/80927233/119919416-67792b00-bf38-11eb-8e85-ffe2a8c69777.jpg)

## Configuración de bioconda e instalación de software

```
# Se configuró el canal de Bioconda para su uso.
- conda config --add channels bioconda 
- conda search -c bioconda fast-qc
- conda search -c bioconda fastqc
# **Se buscaron los softwares de fastqc y trimmomatic en Bioconda, el primer comando no funciona debido a que no existe.**
- conda search -c bioconda trimmomatic 
- conda install -c bioconda fastqc
# **Se instalaron los softwares previamente buscados.**
- conda install -c bioconda trimmomatic 
# **Se creó y accedió a la carpeta SRA_samples para trabajar en ella.**
- mkdir SRA_samples
- cd SRA_samples 
```

## Descarga de la biomuestra desde SRA

```
# Se cargó Nano con el script donwload, ingresando los comandos indicados, guardar y cerrar.
- nano [download.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/SCRIPT/download.sh) 
# Se ejecutó el script download.sh que descargara y validara las secuencias de la biomuestra SRR2006763
- bash download.sh 
# Revisa los contenidos del directorio SRA_samples.
- ls -l -h 
# Se cargó a Nano con el script fdump, ingresando los comandos indicados, guardar y cerrar.
- nano [fdump.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/SCRIPT/fdump.sh)
# El script fdump Permite obtener los archivos fastq de la biomuestra SRR2006763
- bash fdump.sh 
```

### Al finalizar resultara lo siguiente:  

spots read : 2,856,007  
reads read : 5,712,014  
reads written : 5,712,014  

## Comprobación de la integridad de los archivos

```
# md5sum verifica los archivos y redirecciona el resultado. 
- md5sum SRR2006763_1.fastq SRR2006763_2.fastq > md5_samples
```

- Entregando los valores:  
dd0bdf8c722226ea34611941f2391774  SRR2006763_1.fastq  
1c63ca4a6e14de4f93f7621e3e990ec9  SRR2006763_2.fastq  

```
# Se comprobó la integridad de ambas biomuestras, indicando que "La suma coincide" en caso verdadero.
- md5sum -c md5_samples 
```

## Análisis del control de calidad

```
# Se cargó a Nano con el script fastqc, ingresando los comandos indicados, guardar, cerrar y ejecutar. 
- nano [fastqc.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/SCRIPT/fastqc.sh) 
- bash fastqc.sh
```

[Carpeta](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/)  
[SRR2006763_1_fastqc](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR2006763_1_fastqc.html)  
[SRR2006763_2_fastqc](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR2006763_2_fastqc.html)  

Los comandos dentro del Script fastqc.sh procesan los archivos fastq y entrega un informe en formato HTML y un archivo .zip con el reporte y los datos.  
Se puede acceder al directorio home2 desde el puerto [8787](http://200.54.220.141:8787/) donde se ingresa con el mismo usuario y contraseña.  
Desde ahí se pueden ver todos los archivos utilizados y los reportes obtenidos con una interfaz de RStudio como se observa en la imagen.

![fastqc](https://user-images.githubusercontent.com/80927233/121597512-cf3a7600-ca0e-11eb-8c0f-803c4dab20d1.png)

## Filtrado y poda

```
# Se cargó a Nano con el script trimm, ingresando los comandos indicados, guardar y cerrar.
- nano [trimm.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/SCRIPT/trimm.sh) 
# Este script poda a las secuencias con un tamaño menor a 60bp y entrega los archivos filtrados.
- bash trimm.sh 
# Con este comando se pueden descomprimir los archivos obtenidos, no obstante, fastqc puede trabajar sobre archivos comprimidos.
- gunzip SRR20067634_filtered_1P.fastq.gz 
# Se realizó un análisis de calidad de las muestras.
- fastqc  *.fastq.gz 
```

[SRR20067634_filtered_1P](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR20067634_filtered_1P_fastqc.html)  
[SRR20067634_filtered_2P](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR20067634_filtered_2P_fastqc.html)  
[SRR20067634_filtered_1U](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR20067634_filtered_1U_fastqc.html)  
[SRR20067634_filtered_2U](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/FastQC/SRR20067634_filtered_2U_fastqc.html)  

Al comparar los distintos filtrados con los datos originales se puede observar que han sido podados exitosamente.
El comando eliminó todos los fragmentos inferiores a 60bp, no obstante, siguen existiendo muchos fragmentos inferiores a 89bp en el set de datos y por tanto el programa entrega un signo de exclamación en ese resultado. 

