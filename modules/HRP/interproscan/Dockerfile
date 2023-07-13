# Created by: Lee Bergstrand 2018
# Modified by Niklas Schandry 2023
#
# Description:  A Dockerfile for building a container containing InterProScan5.63-95-0.

FROM openjdk:22


ARG IPR=5
ENV IPR $IPR
ENV IPRSCAN "5.63-95.0"

RUN mkdir -p /opt

# Download Interproscan

RUN curl -o /opt/interproscan-$IPRSCAN-64-bit.tar.gz ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/$IPR/$IPRSCAN/interproscan-$IPRSCAN-64-bit.tar.gz
RUN curl -o /opt/interproscan-$IPRSCAN-64-bit.tar.gz.md5 ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/$IPR/$IPRSCAN/interproscan-$IPRSCAN-64-bit.tar.gz.md5

WORKDIR /opt

# Install InterProScan5.
RUN md5sum -c interproscan-$IPRSCAN-64-bit.tar.gz.md5

RUN mkdir -p /opt/interproscan

RUN  tar -pxzf interproscan-$IPRSCAN-64-bit.tar.gz \
    -C /opt/interproscan

RUN rm -f interproscan-$IPRSCAN-64-bit.tar.gz interproscan-$IPRSCAN-64-bit.tar.gz.md5

RUN apt-get update && apt-get -y install libgomp1 python python3 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/interproscan
RUN python3 setup.py -f interproscan.properties

ENTRYPOINT ["/bin/bash", "interproscan/threaded_interproscan.sh"]
