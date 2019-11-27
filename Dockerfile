FROM ubuntu:16.04

LABEL Bio-Rad Support <support@bio-rad.com>

RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    perl \
    parallel \
    pigz \
    default-jdk \
    wget \
    samtools \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev
 

######### FastQC Setup ###############
RUN apt-get install -y \
    openjdk-8-jre-headless

ENV FASTQC_URL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV FASTQC_VERSION 0.11.7
ENV FASTQC_RELEASE fastqc_v${FASTQC_VERSION}.zip
ENV DEST_DIR /opt/
# Make destination directory
RUN mkdir -p $DEST_DIR
RUN echo "This is a test"
# Do this in one command to avoid caching the zip file and its removal in separate layers
RUN curl -SLO ${FASTQC_URL}/${FASTQC_RELEASE} && unzip ${FASTQC_RELEASE} -d ${DEST_DIR} && rm ${FASTQC_RELEASE}
# Make the wrapper script executable
RUN chmod a+x ${DEST_DIR}/FastQC/fastqc
# Include it in PATH
ENV PATH ${DEST_DIR}/FastQC:$PATH
######### End FastQC Setup ###########

######### Cutadapt Setup #############
RUN apt-get install -y \
    python3-pip
RUN pip3 install --upgrade pip setuptools

RUN pip3 install --upgrade pip setuptools
RUN pip3 install --user --upgrade cutadapt
RUN ln -s ~/.local/bin/cutadapt /usr/bin/
######### End Cutadapt Setup #########

######### UMI Tools Setup ############
RUN pip3 install --user --upgrade umi_tools
RUN ln -s ~/.local/bin/umi_tools /usr/bin/umi_tools
######### END UMI Tools Setup ########

######### STAR Setup #################
ENV STAR_VERSION 2.7.0f
# Same deal as above with FASTQC; using precompiled executable
RUN curl -SLO https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && tar -zxvf ${STAR_VERSION}.tar.gz --directory /opt/ && rm ${STAR_VERSION}.tar.gz
ENV PATH /opt/STAR-${STAR_VERSION}/bin/Linux_x86_64:$PATH
######### End STAR Setup #############

######### BEDTOOLS Setup #############
ENV BEDTOOLS_VERSION v2.28.0
RUN mkdir -p ${DEST_DIR}/bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/${BEDTOOLS_VERSION}/bedtools -P ${DEST_DIR}/bedtools/
RUN chmod a+x ${DEST_DIR}/bedtools/bedtools
ENV PATH ${DEST_DIR}/bedtools:$PATH
######### END BEDTOOLS Setup #########

######### PICARD Setup ###############
ENV PICARD_VERSION 2.20.0
RUN mkdir -p ${DEST_DIR}/picard
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -P ${DEST_DIR}/picard/
ENV PATH ${DEST_DIR}/picard:$PATH
######### End PICARD Setup ##########

######### Subread Setup #############
ENV SUBREAD_VERSION 1.6.4 
RUN mkdir -p ${DEST_DIR}/subread
RUN curl -SLO https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz \
    && tar -zxvf subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz --directory /opt/ \
    && rm subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz
ENV PATH ${DEST_DIR}/subread-${SUBREAD_VERSION}-Linux-x86_64/bin/:$PATH
######### End Subread Setup ##########

######### Sambamba Setup #############
ENV SAMBAMBA_VERSION 0.6.9
RUN mkdir -p ${DEST_DIR}/sambamba
RUN curl -SLO https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-static.gz \
    && unpigz sambamba-${SAMBAMBA_VERSION}-linux-static.gz && mv sambamba-${SAMBAMBA_VERSION}-linux-static ${DEST_DIR}/sambamba/sambamba
RUN chmod a+x ${DEST_DIR}/sambamba/sambamba
ENV PATH ${DEST_DIR}/sambamba/:$PATH
######### End Sambamba Setup #########

######### Pysam Setup ################
RUN pip3 install pysam
######### End Pysam Setup ############

######### RPKM/TPM Setup #############
RUN pip3 install pandas
######### End RPKM/TPM Setup #########

######### BigWig Setup ###############
RUN mkdir /opt/tools && wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -O /opt/tools/bedGraphToBigWig && chmod a+x /opt/tools/bedGraphToBigWig
ENV PATH /opt/tools/:$PATH
RUN mkdir -p /ref/hg38/sizes
RUN wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O /ref/hg38/sizes/hg38.chrom.sizes && chmod a+r /ref/hg38/sizes/hg38.chrom.sizes 
######### End BigWig Setup ###########

######### R Setup ###############
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update && apt-get install -y \
    apt-transport-https \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    fonts-freefont-ttf
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
RUN apt-get install -y r-base
RUN Rscript -e 'install.packages(c("dplyr", "knitr", "rmarkdown", "kableExtra", "ggplot2", "plotly", "fastqcr", "data.table", "tibble", "rlist", "tinytex", "webshot"), repos = "http://cran.r-project.org")'
RUN Rscript -e 'tinytex::install_tinytex()'
RUN Rscript -e 'webshot::install_phantomjs()'
RUN Rscript -e 'tinytex::tlmgr_install(pkgs = c("xcolor", "colortbl", "multirow", "wrapfig", "float", "tabu", "varwidth", "threeparttable", "threeparttablex", "environ", "trimspaces", "ulem", "makecell"))'
######### End R Setup ###########

WORKDIR /opt/biorad

COPY . .

# Pull in some ARGS for defining container name
ARG IMAGE_NAME
ARG SOURCE_BRANCH
ARG SOURCE_COMMIT
RUN printf "Container Name: ${IMAGE_NAME:-local}\n" > imageInfo.txt
RUN printf "Source Branch: ${SOURCE_BRANCH:-local}\n" >> imageInfo.txt
RUN printf "Source Commit: ${SOURCE_COMMIT:-local}" >> imageInfo.txt

