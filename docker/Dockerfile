FROM ubuntu:16.04

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && apt-get install -y \
	aptitude \
	autoconf \
	automake \
	bedtools \ 
	build-essential \
	curl \
	dos2unix \
	dpkg-dev \
	g++ \
	gcc \
	git \
	gzip \
	libcairo2-dev \ 
	libcurl4-openssl-dev \
	libexpat1-dev \
	libgsl0-dev \
	libmysqlclient-dev \ 
	libcurl3 \
	libcurl3-dev \ 
	libffi-dev \
	libncurses-dev \
	libncurses5-dev \
	libssl-dev \
	libxt-dev \
	make \
	python3 \
	python3-pip \
	tabix \
	unzip \
	r-base \
	r-base-dev \
	tabix \
	wget \ 
	zlib1g-dev \
   	libxml2-dev \
   	r-cran-rjava 

RUN pip3 install --upgrade pip
RUN pip install synapseclient httplib2 pycrypto aacrgenie
RUN pip install pandas numexpr --upgrade

RUN rm /usr/bin/python 
RUN ln -s /usr/bin/python3 /usr/bin/python 

COPY ./installPackages.R /installPackages.R
RUN Rscript /installPackages.R

ENV HOME=/root
ENV PERL_BASE=$HOME/perl
RUN mkdir -p $PERL_BASE

WORKDIR $PERL_BASE

# Download source tarball into a subfolder named src, and untar:
RUN curl --create-dirs -L -o src/perl-5.22.2.tar.gz http://www.cpan.org/src/5.0/perl-5.22.2.tar.gz
WORKDIR $PERL_BASE/src
RUN tar -zxf perl-5.22.2.tar.gz

# Configure and build, making sure that this non-standard location is baked into the perl binary:

WORKDIR $PERL_BASE/src/perl-5.22.2
RUN ./Configure -des -Dprefix=$PERL_BASE/perl-5.22.2 -Dotherlibdirs=$PERL_BASE/perl-5.22.2/lib/perl5
RUN make
RUN make install

# Install dependencies for Ensembl API, VEP, OncoKB, and other tools we use here at MSKCC:
RUN curl -sL https://cpanmin.us | $PERL_BASE/perl-5.22.2/bin/perl - --notest -l $PERL_BASE/perl-5.22.2 App::cpanminus JSON::Parse XML::Parser XML::Simple LWP LWP::Simple LWP::Protocol::https Archive::Extract Archive::Tar Archive::Zip CGI DBI Time::HiRes DBD::mysql Encode File::Copy::Recursive Perl::OSType Module::Metadata Statistics::Lite Tie::Autotie Tie::IxHash Log::Log4perl FindBin::Real Getopt::Long Catalyst::Runtime Catalyst::Devel List::Util Test::XML::Simple Test::XPath IO::String Bio::Perl version

# Set $PERL5LIB to find these libraries, and set $PATH to use this Perl instead of the system Perl:
ENV PERL5LIB=$PERL_BASE/perl-5.22.2/lib/perl5:$PERL_BASE/perl-5.22.2/lib/perl5/x86_64-linux:$PERL5LIB
ENV PATH=$PERL_BASE/perl-5.22.2/bin:$PATH

ENV VEP_PATH=$HOME/vep
ENV VEP_DATA=$HOME/.vep
RUN mkdir $VEP_PATH $VEP_DATA $VEP_PATH/htslib $VEP_PATH/samtools $VEP_PATH/samtools/bin

ENV PERL5LIB=$VEP_PATH:$PERL5LIB
ENV PATH=$VEP_PATH/htslib:$PATH
ENV PATH=$HOME/vep/samtools/bin:$VEP_PATH/htslib:$PATH

WORKDIR /root
RUN wget https://github.com/mskcc/vcf2maf/archive/v1.6.14.zip
RUN unzip v1.6.14.zip

RUN git clone https://github.com/Sage-Bionetworks/Genie.git
RUN git clone https://github.com/cBioPortal/cbioportal.git

WORKDIR /root/Genie/genie
