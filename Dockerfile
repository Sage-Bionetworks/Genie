FROM ubuntu:hirsute-20211107
ENV DEBIAN_FRONTEND=noninteractive 

# Must install this because gpg not installed
RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		software-properties-common \
		dirmngr \
		wget && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# General sys dependencies
RUN apt-get update && apt-get install -y --allow-unauthenticated --no-install-recommends \
		bedtools \
		dos2unix \
		python3 \
		python3-pip \
		python3-dev \
		git \
		r-base \
		r-base-dev \
		curl \
		# synapser client dependencies
		dpkg-dev \
		zlib1g-dev \
		libssl-dev \
		libffi-dev \
		libcurl4-openssl-dev \
		# VariantAnnotation dependency
		libxml2-dev \
		# Supports data guide creation
		texlive \
		texinfo \
		# texlive-generic-recommended \
		texlive-latex-extra \
		# genome nexus
		openjdk-8-jre && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

#install pandoc 1.19.2.1 (dashboard use)
RUN wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
RUN dpkg -i pandoc-1.19.2.1-1-amd64.deb	

# Only copy most recent changes in code are always installed
# Do not build from local computer
WORKDIR /root/Genie
COPY . .

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=true
RUN Rscript R/install_packages.R

RUN python3 -m pip install --no-cache-dir cython
RUN python3 -m pip install --no-cache-dir -r requirements.txt
RUN python3 -m pip install -e .
# RUN python3 setup.py sdist
# RUN python3 setup.py develop

WORKDIR /root/
# Must move this git clone to after the install of Genie,
# because must update cbioportal
RUN git clone https://github.com/cBioPortal/cbioportal.git
RUN git clone https://github.com/Sage-Bionetworks/annotation-tools.git

WORKDIR /root/Genie
