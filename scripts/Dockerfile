FROM continuumio/miniconda3
WORKDIR /usr/src/app
RUN conda install -c bioconda -c conda-forge samtools bedtools numpy "pysam>=0.16" pysamstats
COPY mapq0_cnv.py .
ENTRYPOINT ["python", "./mapq0_cnv.py"]