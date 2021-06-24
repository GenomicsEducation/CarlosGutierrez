# Introducción al análisis de secuencias NGS - Alineamiento

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
# Se instaló el software BWA.
conda install -c bioconda bwa

# Se instaló SamTools, proceso que demora entre 5 a 10 minutos.
conda install -c bioconda samtools
conda config --add channels bioconda 
conda config --add channels conda-forge 
conda install samtools==1.11

# Para verificar a los directorios que contienen a cada software se utilizó el comando "whereis"
whereis sratoolkit 
whereis samtools
whereis bwa 
# El cual entrega la ruta completa de instalación.
```

## Creación de directorio de trabajo y descarga de datos para alineamiento

```
# Se creó una carpeta denominada "alineamiento" en el directorio home2 del usuario
cd
mkdir alineamiento
# Se ingresó a la carpeta “alineamiento” y se transfirieron los archivos fastq obtenidos en la práctica anterior
cd alineamiento
mv /home2/"**USUARIO**"/SRA_samples/SRR2006763/SRR2006763_1.fastq /home2/"**USUARIO**"/alineamiento/
mv /home2/"**USUARIO**"/SRA_samples/SRR2006763/SRR2006763_2.fastq /home2/"**USUARIO**"/alineamiento/
# Finalmente se utilizó el comando "ls" para listar los contenidos del directorio
ls
```

Se descargo el genoma mitocondrial de *Salmo salar* desde [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=salmo+salar)
como indica la siguiente imagen
