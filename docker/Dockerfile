FROM ubuntu:16.04

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && apt-get install -y \
	autoconf \
	automake \
	make \
	g++ \
	gcc \
	build-essential \ 
	zlib1g-dev \
	libgsl0-dev \
	perl \
	curl \
	git \
	wget \
	unzip \
	tabix \
	libncurses5-dev \
	bedtools \ 
	dos2unix \
	python3 \
	python3-pip \
	r-base \
	r-base-dev \
	dpkg-dev \
	libcurl3 \
	libcurl3-dev \ 
	libffi-dev 

RUN pip3 install --upgrade pip
RUN pip install synapseclient httplib2 pycrypto
RUN pip install pandas numexpr --upgrade

RUN rm /usr/bin/python 
RUN ln -s /usr/bin/python3 /usr/bin/python 

RUN apt-get install -y cpanminus
RUN apt-get install -y libmysqlclient-dev
RUN cpanm CPAN::Meta \
	Archive::Zip \
	DBI \
	Devel::CheckLib \
	DBD::mysql \ 
	JSON \
	DBD::SQLite \
	Set::IntervalTree \
	LWP \
	LWP::Simple \
	Archive::Extract \
	Archive::Tar \
	Archive::Zip \
	CGI \
	Time::HiRes \
	Encode \
	File::Copy::Recursive \
	Perl::OSType \
	Module::Metadata version \
	Bio::Root::Version \
	TAP::Harness \
	Module::Build --force

WORKDIR /root
RUN wget https://github.com/mskcc/vcf2maf/archive/v1.6.14.zip
RUN unzip v1.6.14.zip
RUN git clone https://github.com/cBioPortal/cbioportal.git
RUN git clone https://github.com/Sage-Bionetworks/Genie.git

RUN wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
RUN tar jxf samtools-1.3.tar.bz2
WORKDIR /root/samtools-1.3
RUN make
RUN make install

WORKDIR /root/Genie
RUN python3 setup.py sdist
RUN python3 setup.py develop

WORKDIR /root/Genie/genie
