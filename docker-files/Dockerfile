#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM broadinstitute/gatk:4.1.8.0

# Set maintainer.
LABEL  description='Dockerfile for installing tools and copying files and scripts for GATK JointGenotypeGVCFs' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    bcftools && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar 

# by default /bin/bash is executed
CMD ["/bin/bash"]
