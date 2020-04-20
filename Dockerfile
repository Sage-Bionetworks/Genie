FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive 

# Must install this because gpg not installed
RUN apt-get update && apt-get install -y gnupg2 software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get update && apt-get install -y --allow-unauthenticated \
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
RUN apt-get install -y --allow-unauthenticated \
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

# Supports data guide creation
RUN apt-get install -y texlive \
	texinfo \
	texlive-generic-recommended \
	texlive-latex-extra

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

WORKDIR /root/Genie/genie
