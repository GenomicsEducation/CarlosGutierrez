# Práctica elaboración de un proyecto de genómica aplicada.

## **Autor**
- Carlos Gutierrez Ferreira
- Chileno
- Magíster en Biotecnología

## **Descripción:**  

### **GRCz11** Genome Reference Consortium Zebrafish Build 11

https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6

- Organismo: Danio rerio (zebrafish)
- Nombre infraespecífico: Raza: Tuebingen
- Muestra biológica: SAMN06930106
- Projecto biológico: PRJNA11776
- Remitente: Genome Reference Consortium
- Fecha: 2017/05/09
- Tipo de ensamblaje: Haploide con loci alternativo
- Nivel de ensamblaje: Cromosoma
- Representación del Genoma: Completa
- Categoría RefSeq: Genoma de referencia
- Sinónimos: danRer11
- Acceso de ensamblaje GenBank: GCA_000002035.4 (último)
- Acceso de ensamblaje RefSeq: GCF_000002035.6 (último)
- IDs: 1104621 [UID] 4482478 [GenBank] 4514068 [RefSeq]

La secuencia GRCz11 está compuesta por una secuencia genómica, principalmente clones que fueron secuenciados como parte del Proyecto del Genoma del Pez Cebra en el Instituto Wellcome Trust Sanger. Los productos de PCR y la secuencia shotgun WGS, provenientes principalmente del conjunto WGS31 (CABZ00000000.1), pero también WGS29 (CAAK00000000.1) y WGS32 (CZQB00000000.1) se han agregado cuando es necesario para llenar los vacíos en el genoma. Además, se agregaron secuencias no colocadas desde WGS31 si presentaban al menos 5 kb de secuencia no repetitiva, o si se podían detectar alineaciones con el cDNA que no estaban ya cubiertas por las secuencias colocadas.

  |Característica|Tamaño|
  |----------------|--------------------|
  |Regiones con loci alternativos o parches|607|
  |Largo total de la secuencia|1,373,454,788|
  |Largo total sin vacíos|1,368,765,506|
  |Espacio entre andamios|925|
  |Número de andamios|1,917|
  |Andamio N50|7,379,053|
  |Andamio L50|44|
  |Número de contigs|19,725|
  |Contig N50|1,422,317|
  |Contig L50|219|
  |Número total de cromosomas y plásmidos|25|
  |Número de secuencias componentes (WGS o clone)|31,634|

### SRX6747970: GSM4041277: MPs200vm.3; Danio rerio; OTROS
- 1 ION_TORRENT (Ion S5 XL) Ejecución: 52,087 puntos, 21,5M bases, 16.5Mb

https://www.ncbi.nlm.nih.gov/sra/SRX6747970[accn]

**Remitido por:** NCBI (GEO)  
**Estudio:** rDNA 16S Análisis de secuenciación por amplicones del intestino de pez cebra expuesto a microplásticos de poliestireno  
**Muestra:** MPs200vm.3  
- Organismo: Danio rerio
  
**Librería:**  
- Instrumento: Ion S5 XL
- Estrategia: OTRO
- Fuente: GENOMICA
- Selección: OTRO
- Diseño: ÚNICO

Protocolo de construcción: Se extruyó el contenido del intestino extraído del sujeto. El ADN total en el contenido intestinal se aisló con DNeasy Blood & Tissue Kit (QIAGEN, Alemania). La concentración de ADN se midió usando un espectrofotómetro NanoDrop (NanoDrop One, Thermo, EE. UU.). Después del análisis de la calidad del ADN por electroforesis, el DNA extraído se amplificó utilizando cebadores (cebador sentido: 5'-CCTAYGGGRBGCASCAG-3 '; cebador antisentido: 5'-GGACTACNNGGGTATCTAAT-3') para apuntar a las regiones V3-V4 del gen bacteriano del ARNr 16S . Los productos de PCR obtenidos se purificaron mediante el kit de extracción de gel GeneJET (Thermofisher, EE. UU.) Y se aplicaron para construir la biblioteca de secuenciación. La biblioteca se secuenció en una plataforma Ion S5TM XL.
Atributos del experimento:
Adhesión GEO: GSM4041277

- Ejecuciones: 1 Ejecución, 52.087 puntos, 21,5M bases, 16,5 MB

  |Ejecución |Puntos |Bases| Tamaño| Fecha|
  |--------|--------|-------|--------|-----|
  |SRR10009745 |52,087 |21.5M |16.5Mb |2019-10-01|



## Instalación y configuración de software para acceso remoto y transferencia de archivos, programas WinSCP y PuTTY

### **PuTTY**
![PUTTY](https://user-images.githubusercontent.com/80927233/119920352-34d03200-bf3a-11eb-815e-ce236832d618.jpg)

### **WinSCP**
![WinSCPpage](https://user-images.githubusercontent.com/80927233/119920551-84aef900-bf3a-11eb-8c0f-fb8a2d486099.jpg)

### **Nano**
![NANO](https://user-images.githubusercontent.com/80927233/119920375-3dc10380-bf3a-11eb-885f-92805dd9d2b1.jpg)

## Acceso remoto a servidor Pomeo

![PUTTY3](https://user-images.githubusercontent.com/80927233/119919416-67792b00-bf38-11eb-8e85-ffe2a8c69777.jpg)

## Instalación y configuración Conda, Nano y SRA Toolkit

### **Anaconda en Pomeo (PuTTY) y Nano :: Anaconda.org**
![Conda](https://user-images.githubusercontent.com/80927233/119927124-d493bd00-bf46-11eb-9bf2-0dfac07f129b.jpg)

### **SRA Toolkit**
![SRA](https://user-images.githubusercontent.com/80927233/119927129-d6f61700-bf46-11eb-9d69-f38b276c9a26.jpg)

- https://www.nano-editor.org/dist/
- https://anaconda.org/conda-forge/nano
- https://hpc.nih.gov/apps/sratoolkit.html

## Práctica de Shell y Linux

- bash --version # Versión bash.
- pwd # Nombre del directorio.
- df -hP # Espacio total en el sistema.
- mkdir # Crea un directorio.
- cd # Cambia al directorio.
- cat; less; wc # Lee un archivo e imprime su contenido; recorre el archivo; cuenta líneas, palabras y caracteres.
- ls -l -h # Listado de objetos en un directorio, en detalle y legible.
- rm -r # Remover un fichero o directorio forzando la acción.
- exit # Cierra PuTTY de forma correcta.
- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # Descargar el repositorio de anaconda.
- bash Miniconda3-latest-Linux-x86_64.sh # Ejecutar anaconda.
- source ~/.bashrc # Activacion de miniconda.
- conda list # Contenido de conda.
- conda --version # Version de conda.
- conda install -c conda-forge nano # Instalación de nano en conda.
- nano script1.sh | # !/bin/bash | # Mi primer script
- echo Curso de Genómica # Prueba de script en Nano, guardar con Ctrl+O y cerrrar con Ctrl+X
- bash script1.sh # Ejecutar el script en la terminal.

## Práctica descarga de secuencias NGS con SRA Toolkit

- nano script2.sh | # !/bin/bash | # Descarga y descomprime SRA Toolkit
- http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
- tar -xzf sratoolkit.current-centos_linux64.tar.gz # Script en Nano para descargar e instalar SRA Toolkit, ejecutar con Bash.
- bin/vdb-config --interactive # Ejecutar en el directorio: /sratoolkit.2.10.5-centos_linux64
- fastq-dump --stdout SRR390728 | head -n 8 # Probar que SRAToolkit está trabajando correctamente.
- fastq-dump -X 5 -Z SRR6019464 # Descarga y muestra el contenido de las 5 primeras secuencias del archivo SRR6019464.
- fastq-dump -X 5 SRR6019464 # Descarga el contenido de las 5 primeras secuencias y las almacena en un archivo con formato fastq.
- fastq-dump --gzip --split-3 SRR6019464 # Descarga la biomuestra completa, detener la ejecución luego de unos momentos debido a que son demaciados datos.
- zcat SRR6019464.fastq.gz | echo $((`wc -l`/4)) # Explora y entrega el número de reads descargados.
