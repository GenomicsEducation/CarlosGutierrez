#!/bin/bash
 #SBATCH - J fdump_"**USUARIO**"
 /home2/"**USUARIO**"/sratoolkit.2.11.0-centos_linux64/bin/fasterq-dump /home2/"**USUARIO**"/SRA_samples/SRR2006763/*.sra -O /home2/"**USUARIO**"/SRA_samples/SRR2006763/