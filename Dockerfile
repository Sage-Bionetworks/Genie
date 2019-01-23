FROM ubuntu:16.04

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && apt-get install -y \
	bedtools \ 
	dos2unix \
	wget \ 
	python3 \
	python3-pip \
	#python-pandas \
	git \
	r-base \
	r-base-dev \
	dpkg-dev \
	zlib1g-dev \
	libssl-dev \
	curl \
	libcurl3 \
	libcurl3-dev \ 
	libffi-dev \
	libmariadb-client-lgpl-dev \
	libxml2-dev

RUN pip3 install --upgrade pip
RUN pip install synapseclient httplib2 pycrypto
RUN pip install pandas numexpr --upgrade

RUN rm /usr/bin/python 
RUN ln -s /usr/bin/python3 /usr/bin/python 

COPY docker/installPackages.R /installPackages.R
RUN Rscript /installPackages.R

#install pandoc 1.19.2.1 (dashboard use)
#RUN wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
#RUN dpkg -i pandoc-1.19.2.1-1-amd64.deb	

WORKDIR /root/
#These are only necessary for now
RUN git clone -b develop https://github.com/Sage-Bionetworks/synapser.git
RUN R CMD build synapser/ --no-build-vignettes
RUN R CMD INSTALL synapser_0.0.0.tar.gz

RUN git clone https://github.com/cBioPortal/cbioportal.git

#Only copy necessary files, so no data is copied ever and most recent changes 
#in code are always installed
WORKDIR /root/Genie
COPY ./genie/*.py ./genie/
COPY ./genie/*.sh ./genie/
COPY setup.py ./
COPY MANIFEST.in ./
#No need for oncotree_code_converter script
#COPY ./analyses/clinicalData/oncotree_code_converter.py ./analyses/clinicalData/
COPY ./analyses/genomicData/MAFinBED.R ./analyses/genomicData/
COPY ./analyses/mergeFlag/mergeCheck.R ./analyses/mergeFlag/
RUN python3 setup.py sdist
RUN python3 setup.py develop

WORKDIR /root/Genie/genie
#RUN wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
#RUN gunzip Homo_sapiens.GRCh37.75.gtf.gz
#RUN awk '$3 == "exon" {print}' Homo_sapiens.GRCh37.75.gtf > exon.gtf
#RUN awk '$3 == "gene" {print}' Homo_sapiens.GRCh37.75.gtf > gene.gtf
