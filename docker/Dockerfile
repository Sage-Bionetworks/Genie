FROM ubuntu:14.04

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
	libncurses5-dev

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

#ENV VEP_PATH=$HOME/vep
#ENV VEP_DATA=$HOME/.vep
#RUN mkdir $VEP_PATH $VEP_DATA $VEP_PATH/htslib $VEP_PATH/samtools $VEP_PATH/samtools/bin

#ENV PERL5LIB=$VEP_PATH:$PERL5LIB
#ENV PATH=$VEP_PATH/htslib:$PATH
#ENV PATH=$HOME/vep/samtools/bin:$VEP_PATH/htslib:$PATH

WORKDIR /root
RUN wget https://github.com/mskcc/vcf2maf/archive/v1.6.14.zip
RUN unzip v1.6.14.zip

RUN git clone https://github.com/Sage-Bionetworks/Genie.git
RUN git clone https://github.com/cBioPortal/cbioportal.git

WORKDIR /root/Genie/genie
