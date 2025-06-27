# Base image
FROM rocker/r-base:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    pandoc \
    jags \
    automake \
    autoconf \
    libtool \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libxml2-dev \
    libglpk-dev \
    libicu-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libpng-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \   
    make \
    g++ \
    && apt-get clean

# Set CRAN mirror and install R packages
COPY packages.R /tmp/
RUN Rscript -e "options(repos = c(CRAN='http://cran.us.r-project.org')); source('/tmp/packages.R')"

# Copy project files
COPY . /home/project/
WORKDIR /home/project/

# Optional command to run script automatically
# CMD ["Rscript", "analysis.R"]
