# Práctica Introducción a Linux para genómica

## **Autor**
### Carlos Gutierrez Ferreira  
- Chileno
- Magíster en Biotecnología

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

### Práctica creación de script usando el editor de texto NANO
- nano [script1.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Linux-Genomica/SCRIPT/script1.sh) # Prueba de script en Nano, guardar con Ctrl+O y cerrrar con Ctrl+X
- bash script1.sh # Ejecutar el script1 en la terminal.

## Práctica descarga de secuencias NGS con SRA Toolkit

- nano [script2.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Linux-Genomica/SCRIPT/script2.sh) # Carga Nano con el script2, guardar y cerrar
- bash script2.sh # Ejecuta el script2 en Nano para descargar e instalar SRA Toolkit.
- cd sratoolkit.2.11.0-centos_linux64 # Cambiar el directorio para ejecutar SRA Toolkit.
- bin/vdb-config --interactive # Ejecuta vdb-config en el directorio bin de forma interactiva, Tab y Enter para salir
- bin/fastq-dump --stdout SRR390728 | head -n 8 # Probar que SRAToolkit está trabajando correctamente.
- bin/fastq-dump -X 5 -Z SRR6019464 # Descarga y muestra el contenido de las 5 primeras secuencias del archivo SRR6019464.
- bin/fastq-dump -X 5 SRR6019464 # Descarga el contenido de las 5 primeras secuencias y las almacena en un archivo con formato fastq.
- bin/fastq-dump --gzip --split-3 SRR6019464 # Descarga la biomuestra completa, detener la ejecución luego de unos momentos debido a que son demasiados datos.
- zcat SRR6019464_1.fastq.gz | echo $((`wc -l`/4)) # Explora y entrega el número de reads descargados en el 1er archivo.
- # Ej: 23772

