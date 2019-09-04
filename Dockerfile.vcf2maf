FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y gnupg2 software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'


# GENIE dependencies
RUN apt-get update && apt-get install -y --allow-unauthenticated \
	git \
	r-base-core \
	r-base-dev \
	python3 \
	python3-pip \
	bedtools \ 
	dos2unix \
	g++ \
	gcc 

# VCF2MAF dependencies
RUN apt-get install -y --allow-unauthenticated \
	build-essential \
	cpanminus \
	curl \
	wget \
	unzip \
	automake \
	samtools \
	tabix 

RUN apt install -y --allow-unauthenticated \
	libmysqlclient-dev \
	libncurses5-dev \
	zlib1g-dev \
	libgsl0-dev \
	libexpat1-dev \
	libgd-dev

RUN pip3 install --upgrade pip
RUN pip install synapseclient httplib2 pycrypto PyYAML
RUN pip install pandas numexpr --upgrade

RUN rm /usr/bin/python 
RUN ln -s /usr/bin/python3 /usr/bin/python 

RUN cpanm --notest LWP::Simple DBI DBD::mysql Archive::Zip Archive::Extract HTTP::Tiny Test::Simple
RUN cpanm --notest File::Copy::Recursive Perl::OSType Module::Metadata version TAP::Harness CGI Encode
RUN cpanm --notest CPAN::Meta JSON DBD::SQLite Set::IntervalTree Archive::Tar Time::HiRes Module::Build
RUN cpanm --notest Bio::Root::Version

WORKDIR /root
RUN wget https://github.com/mskcc/vcf2maf/archive/v1.6.17.zip
RUN unzip v1.6.17.zip

WORKDIR /root/Genie
COPY ./ ./
RUN python3 setup.py sdist
RUN python3 setup.py develop

WORKDIR /root
# Must move this git clone to after the install of Genie,
# because must update cbioportal
RUN git clone https://github.com/cBioPortal/cbioportal.git

RUN pip install synapseclient --upgrade
WORKDIR /root/Genie/genie
