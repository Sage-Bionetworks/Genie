# TODO: update this
FROM ubuntu:jammy-20250714
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
		r-base-core=4.3.3-1.2204.0 \
		r-base-dev=4.3.3-1.2204.0 \
		cmake \
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
		openjdk-11-jre \
		# This is for reticulate
		# TODO: update this
		python3.10-venv && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

#install pandoc 1.19.2.1 (dashboard use)
# RUN wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
# RUN dpkg -i pandoc-1.19.2.1-1-amd64.deb
RUN wget https://github.com/jgm/pandoc/releases/download/3.0.1/pandoc-3.0.1-1-amd64.deb
RUN dpkg -i pandoc-3.0.1-1-amd64.deb



# Only copy most recent changes in code are always installed
# Do not build from local computer
WORKDIR /root/Genie
# Copy install packages first, because R installation takes
# a long time and unless there are changes to the actual
# R packages used.  So the only files copied over are
# renv/ renv.lock and the installation R script
COPY R/install_packages.R R/install_packages.R
COPY renv/ renv/
COPY renv.lock renv.lock
RUN Rscript R/install_packages.R

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=true
COPY . .
RUN echo "source('renv/activate.R')" >> .Rprofile

# TODO Must include R/ and templates/ within the
# genie/ directory to use MANIFEST.in
# For now, install using develop parameter so that
# the package is called from the directory
RUN pip3 install --no-cache-dir -r requirements.txt
RUN python3 setup.py develop
# RUN pip3 install --no-cache-dir .

WORKDIR /root/
# Must move this git clone to after the install of Genie,
# because must update cbioportal
RUN git clone https://github.com/cBioPortal/cbioportal.git -b v5.3.19
RUN git clone https://github.com/Sage-Bionetworks/annotation-tools.git -b 0.0.6


WORKDIR /root/Genie
