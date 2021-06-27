# Introducción al análisis de secuencias NGS - Llamado de variantes

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
# Se instaló el software gatk4.
conda install -c bioconda gatk4
# gatk4 tiene un tamaño de 277 MB aproximadamente, esto demoró su instalación en aproximadamente 3 minutos.
# Se utilizó el software VCFtools, previamente instalado y utilizado en prácticas anteriores.
conda install -c bioconda vcftools
```

## Creación del directorio de trabajo “variant_call” y preparación de los archivos para el llamado de variantes.

```
# Se creó un nuevo directorio para analizar las variantes con el comando “mkdir” ingresando al mismo mediante “cd”
mkdir variant_call
cd variant_call
# Para obtener el programa Picard se utilizó el comando “wget”
wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
# Este programa corresponde a un grupo de comandos en Java para manipular grandes sets de datos.
# Se copiaron los archivos de alineamiento .bam y .bam.bai analizados previamente, en el nuevo directorio.
cp /home2/home2/"**USUARIO**"/alineamiento/SRR2006763.sort.bam /home2/"**USUARIO**"/variant_call
cp /home2/home2/"**USUARIO**"/alineamiento/SRR2006763.sort.bam.bai /home2/"**USUARIO**"/variant_call
# Se obtuvo el genoma de referencia ref_genome.fna
Wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/8030/100/GCF_000233375.1_ICSASG_v2/
# Indexando el genoma de referencia con SamTools para obtener el archivo ref_genome.fna.fai
samtools faidx ref_genome.fna
```

## Exploración del genoma de referencia
[Controles "less"](https://www.thegeekstuff.com/2010/02/unix-less-command-10-tips-for-effective-navigation)  

```
# Con el comando "less" se puede explorar el genoma de referencia y su indice generado con SamTools.
less ref_genome.fna
less ref_genome.fna.fai
# Así como "head" y "tail" entregan los primeras y últimas líneas de los archivos analizados.
head -n 20 ref_genome.fna
head -n 30 ref_genome.fna.fai 
tail -n 20 ref_genome.fna
tail -n 20 ref_genome.fna.fai  
  
# Se investigaron los cromosomas del salmón con el comando "grep"
grep 'NC_' ref_genome.fna
grep -c 'NC_' ref_genome.fna
grep 'NC_' ref_genome.fna.fai
grep -c 'NC_' ref_genome.fna.fai
# Notar que al utilizar el índice del genoma de referencia el proceso es más rápido.

# Se Investigaron los contigs no mapeados del genoma de referencia de salmón, obteniendo 232125 contigs.
grep -c 'NW_' ref_genome.fna
grep -c 'NW_' ref_genome.fna.fai
```

## Llamado de variantes

```
# Posteriormente se ejecutó el programa Picard para crear el diccionario de secuencias, Obteniendo el archivo ref_genome.dict
java -jar picard.jar CreateSequenceDictionary R=ref_genome.fna O=ref_genome.dict
# Se añadieron los grupos de reads al archivo sort.bam obteniendo el archivo SSRR2006763_sorted_RG.bam
java -jar picard.jar AddOrReplaceReadGroups I=SRR2006763.sort.bam O=SSRR2006763_sorted_RG.bam ID=sample LB=Paired-end PL=Illumina PU=Unknown SM=sample
# Indexando el archivo anterior con SamTools.
samtools index SSRR2006763_sorted_RG.bam
# Finalmente se utilizó el comando HaplotypeCaller del software GATK para obtener el archivo raw_variants.vcf
gatk HaplotypeCaller -R ref_genome.fna -I SSRR2006763_sorted_RG.bam -O raw_variants.vcf
```
- Notar que este proceso demora aproximadamente 54 minutos, debido a que se están comparando las variantes ingresadas con el genoma de referencia para llamar a los haplotipos del genoma en estudio, el numero indicado entre el nombre del contig y la localización de la variante en bp, corresponde al tiempo transcurrido en minutos, ej: 53.8

[raw_variants.vcf](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Analisis-variantes/raw_variants.vcf)
[raw_variants.vcf.idx](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Analisis-variantes/raw_variants.vcf.idx)

```
# Al terminar se obtuvieron los archivos raw_variants.vcf y raw_variants.vcf.idx
# Se puede explorar el archivo con los comandos “less, head y tail”
less raw_variants.vcf
head -n 30 raw_variants.vcf
tail -n 30 raw_variants.vcf
# Se utilizó el comando "grep" para contar el número de líneas en el encabezado y la cantidad de variantes detectadas en el archivo.
grep "^#" -c  raw_variants.vcf
grep "^#" -c -v raw_variants.vcf
# Obteniendo 232180 líneas y 57820 variantes respectivamente.
# Se utilizó el sigueinte comando para obtener el nombre de la muestra utilizada en el llamado de las variantes.
grep "^#CHROM" raw_variants.vcf | cut -f 10-
# Para comprender los datos obtenidos se analizaron las columnas.
grep "^#CHROM" raw_variants.vcf
# Luego las primeras 10 variantes.
grep "^#" -v raw_variants.vcf | head
# Finalmente se analizo la codificacion de las columnas INFO y FORMAT que contine la información del genotipo de la muestra.
grep "##INFO" raw_variants.vcf
grep "##FORMAT" raw_variants.vcf
```

- Para finalizar se realizó una extracción de las variantes con alta calidad, la cual está determinada por el número de reads que presenta cada una, para esto se utilizó un comando en formato pipeline para extraer las variantes que presentaran calidad sobre 100 (en la columna 6), imprimiendo el resultado en el archivo hq_variant.txt
[hq_variant.txt](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Analisis-variantes/hq_variant.txt)

```
grep -v "#" raw_variants.vcf | awk '{if ($6 > 100 ) print }' > hq_variant.txt
# Se analizo el archivo con las variantes de alta calidad con el comando "grep"
grep "NC_" -c -v hq_variant.txt
grep "NW_" -c -v hq_variant.txt
# Obteniendo 1736 contigs mapeados y 10669 sin mapear.
```

## Análisis de variantes con vcftools
[hq.freqs.txt](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Analisis-variantes/hq.freqs.txt)
[Resumen VCFtools](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Analisis-secuencias-NGS/Analisis-variantes/out.log)

```
# Se contaron los individuos y variantes en el archivo raw_variants.vcf obteniendo 1 individuo y 57820 variantes.
vcftools --vcf raw_variants.vcf
# Se determinaron las frecuencias de todos los alelos, obteniendo el archivo hq.freqs.txt
vcftools --vcf raw_variants.vcf --freq -c > hq.freqs.txt
# Es posible filtrar los datos entregados en vcftools por algún cromosoma particular incluyendo el argumento –chr
# o podemos excluir el genoma mitocondrial con –not-chr
vcftools --vcf raw_variants.vcf --chr NC_027300.1
vcftools --vcf raw_variants.vcf --not-chr NC_001960.1
# el argumento --freq entrega la frecuencia de las variantes en el cromosoma analizado.
vcftools --vcf raw_variants.vcf --freq -c --chr NC_027300.1
vcftools --vcf raw_variants.vcf --freq -c --not-chr NC_001960.1
# Finalmente se extrajeron los INDELS con el argumento –keep-only-indel y los SNP con –remove-indels
vcftools --vcf raw_variants.vcf --freq -c --chr NC_027300.1 --keep-only-indels
vcftools --vcf raw_variants.vcf --freq -c --chr NC_027300.1 --remove-indels
```

## Visualización de variantes con IGV

![IGV variante](https://user-images.githubusercontent.com/80927233/123532587-822bf480-d6dc-11eb-8794-eca7e2933c32.jpg)
