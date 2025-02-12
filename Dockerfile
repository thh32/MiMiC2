FROM mambaorg/micromamba:git-99f55e2-cuda12.2.2-ubuntu22.04

USER root

RUN apt-get update && apt-get -y install git && rm -rf /var/lib/apt/lists/*

USER mambauser

# RUN git clone https://github.com/ebi-jlu8/MiMiC2.git --verbose

WORKDIR /MiMiC2
COPY . . 

RUN micromamba install -y --no-channel-priority -n base \
    -c bioconda -c conda-forge \
    "python=3.11" "numpy=1.24.3" "scipy=1.10.1" \
    "conda-forge::matplotlib-base" "seaborn=0.13.0" \
    "pandas=1.5.3" "tqdm" "glpk" "liblzma-devel" \
    "libzlib" "zlib"

RUN micromamba install -y --no-channel-priority -n base \
    -c bioconda -c conda-forge \
    "r-essentials" "r-base"

RUN micromamba run -n base Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
RUN micromamba run -n base Rscript -e 'remotes::install_github("SysBioChalmers/sybil")'
RUN micromamba run -n base Rscript -e 'remotes::install_github("euba/bacarena")'

ENV PATH=$PATH:/MiMiC2
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python"]

# Usage:
# docker build -t mimic2 .
# docker run mimic2 MiMiC2.py
# docker run mimic2 MiMiC2-BUTLER.py

# docker run -it --entrypoint /bin/bash mimic2

# Clear cache space
# docker system prune