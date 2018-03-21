FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y build-essential wget unzip python2.7 python-dev git python-pip bats awscli curl \
					   libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add files
RUN mkdir /usr/metaspades
ADD requirements.txt /usr/metaspades

# Install python requirements
RUN pip install -r /usr/metaspades/requirements.txt && rm /usr/metaspades/requirements.txt


# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz && \
	tar xzvf sratoolkit.2.8.2-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.8.2-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.8.2-ubuntu64.tar.gz

# Install metaSPAdes
RUN cd /usr/metaspades/ && \
	wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz && \
	tar xzvf SPAdes-3.11.1-Linux.tar.gz && \
	rm SPAdes-3.11.1-Linux.tar.gz
ENV PATH="/usr/metaspades/SPAdes-3.11.1-Linux/bin:${PATH}"


# Install Prokka
RUN mkdir /usr/prokka && \
	cd /usr/prokka && \
	wget https://github.com/tseemann/prokka/archive/v1.12.tar.gz && \
	tar xzvf v1.12.tar.gz && \
	prokka-1.12/bin/prokka --setupdb
ENV PATH="/usr/prokka/prokka-1.12/bin:${PATH}"


# Install tbl2asn
RUN cd /usr/prokka && \
	wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
	gunzip linux64.tbl2asn.gz && \
	chmod +x linux64.tbl2asn && \
	mv /usr/prokka/linux64.tbl2asn /usr/local/bin/tbl2asn


# Add the run script to the PATH
ADD run.py /usr/metaspades
ADD helpers /usr/metaspades/helpers
RUN cd /usr/metaspades && \
	chmod +x run.py && \
	ln -s /usr/metaspades/run.py /usr/bin/


# Run tests and then remove the folder
ADD tests /usr/metaspades/tests
RUN bats /usr/metaspades/tests/ && rm -r /usr/metaspades/tests/
