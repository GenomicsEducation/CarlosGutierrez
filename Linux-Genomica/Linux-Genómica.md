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

```
# Versión bash disponible.
- bash --version 
# Para visualizar el nombre del directorio.
- pwd 
# Espacio total en el sistema.
- df -hP 
# Crea un directorio.
- mkdir 
# Cambia al directorio indicado, si se encuentra en blanco procede al principal, en este caso home2.
- cd 
# Lee un archivo e imprime su contenido; recorre el archivo; cuenta líneas, palabras y caracteres.
- cat; less; wc 
# Listado de objetos en un directorio, en detalle y legible.
- ls -l -h 
# Remover un fichero o directorio forzando la acción.
- rm -r 
# Cierra PuTTY de forma correcta.
- exit 
# Se descargó el repositorio de anaconda.
- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
# Se ejecutó anaconda.
- bash Miniconda3-latest-Linux-x86_64.sh 
# Se activó miniconda.
- source ~/.bashrc 
# Contenido de conda.
- conda list 
# Versión de conda.
- conda --version 
# Se instaló nano en conda.
- conda install -c conda-forge nano 
```

### Práctica creación de script usando el editor de texto NANO

```
# Prueba de script en Nano, Se ingresaron los comandos indicados en el archivo, guardar con Ctrl+O y cerrrar con Ctrl+X
- nano [script1.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Linux-Genomica/SCRIPT/script1.sh) 
# Se ejecutó el script1 en la terminal.
- bash script1.sh 
```

## Práctica descarga de secuencias NGS con SRA Toolkit

```
# Se ejecutó Nano con el script2, ingresando los comandos indicados, guardar y cerrar.
- nano [script2.sh](https://github.com/GenomicsEducation/CarlosGutierrez/blob/main/Linux-Genomica/SCRIPT/script2.sh) 
# Se ejecutó el script2 en Nano para descargar e instalar SRA Toolkit.
- bash script2.sh 
# Se cambió el directorio para ejecutar SRA Toolkit.
- cd sratoolkit.2.11.0-centos_linux64 
# Ejecuta vdb-config en el directorio bin de forma interactiva, Tab y Enter para salir.
- bin/vdb-config --interactive 
# Probar que SRAToolkit está trabajando correctamente, lo cual entrega las primeras 8 lineas del bioproyecto.
- bin/fastq-dump --stdout SRR390728 | head -n 8 
# Se descargaron las primeras 5 secuencias del archivo SRR6019464, mostrandolas en la interfaz.
- bin/fastq-dump -X 5 -Z SRR6019464 
# Descarga el contenido de las 5 primeras secuencias y las almacena en un archivo con formato fastq.
- bin/fastq-dump -X 5 SRR6019464 
# Descarga la biomuestra completa, detener la ejecución luego de unos momentos debido a que son demasiados datos.
- bin/fastq-dump --gzip --split-3 SRR6019464 
# Explora y entrega el número de reads descargados en el 1er archivo.
- zcat SRR6019464_1.fastq.gz | echo $((`wc -l`/4)) 
```

- ### Ejemplo: 23772
