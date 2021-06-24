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
# Se creó una carpeta denominada "alineamiento" en el directorio home2 del usuario.
cd
mkdir alineamiento
# Se ingresó a la carpeta “alineamiento” y se transfirieron los archivos fastq obtenidos en la práctica anterior.
cd alineamiento
mv /home2/"**USUARIO**"/SRA_samples/SRR2006763/SRR2006763_1.fastq /home2/"**USUARIO**"/alineamiento/
mv /home2/"**USUARIO**"/SRA_samples/SRR2006763/SRR2006763_2.fastq /home2/"**USUARIO**"/alineamiento/
# Finalmente se utilizó el comando "ls" para listar los contenidos del directorio.
ls
```

- Se descargo el genoma mitocondrial de referencia proveniente de *Salmo salar* desde [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=salmo+salar) en formato [Fasta](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Fasta/mt.fasta)
- El archivo de texto obtenido se cargo al directorio "alineamiento" mediante la aplicación [WinSCP](https://winscp.net/eng/download.php)
- Ingresando al servidor con los datos indicados en la imagen:

![Winscp](https://user-images.githubusercontent.com/80927233/123209976-5ff16700-d48f-11eb-9183-9e165dc07b4f.png)

## Indexación del genoma de referencia y alineamiento  
[aln_mt.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/SCRIPT/aln_mt.sh)
[muestra_stat.txt](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Fasta/muestra_stat.txt)
[Controles "less"](https://www.thegeekstuff.com/2010/02/unix-less-command-10-tips-for-effective-navigation)

```
# Se realizó la indexación del genoma mitocondrial de referencia con el software BWA.
bwa index mt.fasta
# Se utilizo el script "aln_mt.sh" para alinear las secuencias fastq previamente analizadas, con el genoma mitocondrial.
nano aln_mt.sh
bash aln_mt.sh

# Se observó el resultado con el comando "less" recordando los controles para su visualización.
less SRR2006763.sam
# Finalmente se realizó un análisis estadístico con salida en un archivo de texto.
samtools flagstat SRR2006763.bam > muestra_stat.txt
```

## Visualización del alineamiento con IGV

Se utilizó el software IGV para visualizar el alineamiento obtenido frente al genoma mitocondrial de referencia, para esto se descargaron los archivos [SRR2006763.sort.bam](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Fasta/SRR2006763.sort.bam) y [SRR2006763.sort.bam.bai]( https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Fasta/SRR2006763.sort.bam.bai) obtenidos desde WinSCP y se cargaron al programa IGV junto con el genoma de referencia para *Salmo salar* y su sector mitocondrial como se observa en la siguiente imagen:

![IGVSSA](https://user-images.githubusercontent.com/80927233/123213039-86b19c80-d493-11eb-9535-3a11104b2daa.png)
