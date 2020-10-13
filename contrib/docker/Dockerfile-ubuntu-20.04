FROM ubuntu:20.04 as base
LABEL Maintainer="Rob Egan<RSEgan@lbl.gov>"

WORKDIR /root

# This is necessary because the upgrade sometimes prompts for input
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/usr/local/bin:${PATH}

RUN apt-get update  && \
    apt-get install -y libgomp1 zlib1g zlib1g-dev libbz2-1.0 git curl wget python3 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 10 && \
    apt-get upgrade -y && \
    echo "Cleanup apt files" && \
    apt-get autoremove -y && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* && \
    find / -name '*.pyc' -exec rm {} \; && \
    echo "Done"

#    
# build env with build tools, compilers, etc
#
FROM base as buildenv
LABEL Maintainer="Rob Egan<RSEgan@lbl.gov>"

ARG MHM2_VER

RUN apt-get update  && \
    apt-get install -y build-essential cmake && \
    echo "Cleanup apt files" && \
    apt-get autoremove -y && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* && \
    find / -name '*.pyc' -exec rm {} \; && \
    echo "Done installing build-essential"

# Copy source tree to mhm2.  If not run from the mhm2 top level, then a git clone of ${MHM2_VER} will occur
COPY . mhm2

RUN echo "Build MHM2" && \
    cd mhm2 && \
    ( [ -z "${MHM2_VER}" -o -f src/mhm2.py ] || ( \
      git clone https://bitbucket.org/berkeleylab/mhm2.git . && \
      git checkout ${MHM2_VER} && \
      git submodule init && \
      git submodule update ) ) && \
    git describe && \
    echo "Build upcxx" && \
    mkdir build-upcxx && \
    cd build-upcxx && \
    TMPDIR=$(pwd) ../upcxx-utils/contrib/install_upcxx.sh /usr/local && \
    cd .. && \
    echo ${PATH} && \
    which upcxx && \
    mkdir build && \
    cd build && \
    cmake -DENABLE_CUDA=off -DMHM2_VECTORS=off .. && \
    make -j 8 all check install && \
    echo "Done building MHM2"

#
# run env with upcxx and mhm2 installs from buildenv
#
FROM base as runenv
LABEL Maintainer="Rob Egan<RSEgan@lbl.gov>"

COPY --from=buildenv /usr/local/ /usr/local/

ENV PATH=/usr/local/bin:${PATH}

CMD ["/bin/bash"]