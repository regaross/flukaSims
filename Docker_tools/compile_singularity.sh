#!/usr/bin/bash
PATH=/usr/local/go/bin:$PATH
source /root/.bashrc
cd /app/singularity-ce-4.1.0/
export VERSION=4.1.0 && /root/singularity-ce-${VERSION}/mconfig && \
make -C /root/singularity-ce-${VERSION}/builddir && \
make -C /root/singularity-ce-${VERSION}/builddir install
