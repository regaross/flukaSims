# Use a base image with a suitable Linux distribution
FROM ubuntu:23.10

# # Set the working directory
  WORKDIR /root

RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    git \
    libfuse-dev \
    libglib2.0-dev \
    libseccomp-dev \
    libtool  make \
    pkg-config \
    runc \
    squashfs-tools \
    squashfs-tools-ng \
    uidmap \
    wget \
    zlib1g-dev \
    debootstrap

RUN export VERSION=1.21.0 OS=linux ARCH=amd64 && \
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
rm go$VERSION.$OS-$ARCH.tar.gz

RUN echo 'export PATH=/usr/local/go/bin:$PATH' >> /root/.bashrc && \
    /bin/bash -c "source /root/.bashrc"

RUN export VERSION=4.1.0 && \
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
tar -xzf singularity-ce-${VERSION}.tar.gz 

COPY compile_singularity.sh /root/singularity-ce-4.1.0/

RUN cd singularity-ce-4.1.0 && ./compile_singularity.sh


# Update and install necessary packages
RUN apt-get update && \
    apt-get install -y \
        software-properties-common \
        python3.8 \
        python3-tk \
        xauth \
        x11-apps \
        x11-utils


# Install any other dependencies you need for your GUI application

# Install X11 and other necessary tools for GUI
RUN apt-get update && \
    apt-get install -y \

RUN wget -q -O - https://flair.web.cern.ch/flair/download/ubuntu/flair.gpg |  tee /etc/apt/trusted.gpg.d/flair.asc && \
 add-apt-repository 'deb [arch=all,amd64] https://cern.ch/flair/download/ubuntu/23.10 /'

RUN apt-get install -y flair
