FROM python:3.14-slim-trixie

SHELL ["/bin/bash", "-euxo", "pipefail", "-c"]


# Install R and necessary dependencies and Git
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates \
    cmake pkg-config libnlopt-dev libuv1-dev \
    git \
    r-base-core=4.5.0-3 \
    r-base-dev=4.5.0-3 \
    && rm -rf /var/lib/apt/lists/*

# Install uv
ADD https://astral.sh/uv/install.sh /uv-installer.sh
RUN sh /uv-installer.sh && rm /uv-installer.sh
ENV PATH="/root/.local/bin:$PATH"

# Install Python packages
RUN uv pip install --system pandas numpy scipy

# Install CRAN packages used by the analysis scripts. Do not install Suggests here:
# ggpubr's optional dependency tree pulls in many build-time extras.
RUN Rscript -e "\
packages <- c('ggplot2','reshape2','RColorBrewer','dplyr','ggpubr','ggsci','gridExtra','gplots','ggrepel','scales','RUnit','whereami'); \
install.packages(packages, repos='https://cloud.r-project.org', dependencies=c('Depends','Imports','LinkingTo')); \
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]; \
if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', ')); \
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
BiocManager::install('DESeq2') \
"