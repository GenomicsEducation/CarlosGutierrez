# Introducción al análisis de secuencias NGS

## **Autor**
### Carlos Gutierrez Ferreira  
- Chileno
- Magíster en Biotecnología

## Conexión al servidor Pomeo

![PUTTY3](https://user-images.githubusercontent.com/80927233/119919416-67792b00-bf38-11eb-8e85-ffe2a8c69777.jpg)

## Configuración de bioconda e instalación de software

- conda config --add channels bioconda # Configura el canal de Bioconda para su uso.
- conda search -c bioconda fast-qc
- conda search -c bioconda fastqc
- conda search -c bioconda trimmomatic
- # Busca los softwares de fastqc y trimmomatic en Bioconda, el primer comando no funcionara debido a que no existe.
- conda install -c bioconda fastqc
- conda install -c bioconda trimmomatic
- # Instala los softwares previamente buscados.
- mkdir SRA_samples
- cd SRA_samples
- # Crea y accede a la carpeta SRA_samples para trabajar en ella.

## Descarga de la biomuestra desde SRA

- nano [download.sh]() # Carga Nano con el script donwload, ingresar los comandos indicados, guardar y cerrar.
- bash download.sh # Ejecuta el script download.sh que descargara y validara las secuencias de la biomuestra SRR2006763
- ls -l -h # Revisa los contenidos del directorio SRA_samples.
- nano [fdump.sh]() # Carga a Nano con el script fdump, ingresar los comandos indicados, guardar y cerrar.
- bash fdump.sh # El script fdump Permite obtener los archivos fastq de la biomuestra SRR2006763

### Al finalizar resultara lo siguiente:  

spots read : 2,856,007  
reads read : 5,712,014  
reads written : 5,712,014  

## Comprobación de la integridad de los archivos

- md5sum SRR2006763_1.fastq SRR2006763_2.fastq > md5_samples # md5sum verifica los archivos y redirecciona el resultado entregando los valores:  
dd0bdf8c722226ea34611941f2391774  SRR2006763_1.fastq  
1c63ca4a6e14de4f93f7621e3e990ec9  SRR2006763_2.fastq  

- md5sum -c md5_samples # Comprueba la integridad de ambas biomuestras indicando que "La suma coincide" en caso verdadero.

## Análisis del control de calidad

- nano [fastqc.sh]() # Carga a Nano con el script fastqc, ingresar los comandos indicados, guardar y cerrar.  
- bash fastqc.sh
  
Este comando procesa los archivos fastq y entrega un informe en formato HTML y un archivo .zip con el reporte y los datos.  
Se puede acceder al directorio home2 desde el puerto [8787](http://200.54.220.141:8787/) donde se ingresa con el mismo usuario y contraseña.  
Desde ahí se pueden ver todos los archivos utilizados y los reportes obtenidos con una interfaz de RStudio como se observa en la imagen.



## Filtrado y poda

- nano [trimm.sh]() # Carga a Nano con el script trimm, ingresar los comandos indicados, guardar y cerrar.
- bash trimm.sh # Este script poda a las secuencias con un tamaño menor a 60bp y entrega los archivos filtrados.
- gunzip SRR20067634_filtered_1P.fastq.gz # Con este comando se pueden descomprimir los archivos obtenidos, no obstante, fastqc puede trabajar sobre archivos comprimidos.
- fastqc  *.fastq.gz # Realiza un análisis de calidad de las muestras.

[SRR20067634_filtered_1P]()
[SRR20067634_filtered_2P]()
[SRR20067634_filtered_1U]()
[SRR20067634_filtered_2U]()

Al comparar los distintos filtrados con los originales se puede observar que han sido podados exitosamente.
El comando elimino todos los fragmentos inferiores a 60bp, no obstante, siguen existiendo muchos fragmentos inferiores a 89bp en el set de datos y por tanto el programa entrega un signo de exclamación en ese resultado. 

