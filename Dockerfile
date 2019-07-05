FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive 

# Must install this because gpg not installed
RUN apt-get update && apt-get install -y gnupg2

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && apt-get install -y \
	bedtools \ 
	dos2unix \
	wget \ 
	python3 \
	python3-pip \
	git \
	r-base-core \
	r-base-dev \
	curl 
#synapser client dependencies
RUN apt-get install -y \
    dpkg-dev \
	zlib1g-dev \
	libssl-dev \
	libffi-dev \
	libcurl4-openssl-dev \
# VariantAnnotation dependency
	libxml2-dev
# 	libcurl3 \
# 	libcurl3-dev \ 
# 	libmariadb-client-lgpl-dev \

RUN pip3 install --upgrade pip
RUN pip install synapseclient httplib2 pycrypto PyYAML
RUN pip install pandas numexpr --upgrade

# RUN rm /usr/bin/python 
# RUN ln -s /usr/bin/python3 /usr/bin/python 

#install pandoc 1.19.2.1 (dashboard use)
RUN wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
RUN dpkg -i pandoc-1.19.2.1-1-amd64.deb	

COPY docker/installPackages.R /installPackages.R
RUN Rscript /installPackages.R

# Only copy most recent changes in code are always installed
# Do not build from local computer
WORKDIR /root/Genie
COPY ./ ./
RUN python3 setup.py sdist
RUN python3 setup.py develop

WORKDIR /root/
# Must move this git clone to after the install of Genie,
# because must update cbioportal
RUN git clone https://github.com/cBioPortal/cbioportal.git

RUN pip install synapseclient --upgrade
WORKDIR /root/Genie/genie
