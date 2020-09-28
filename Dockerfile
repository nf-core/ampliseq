FROM nfcore/base:1.10.2
LABEL authors="Daniel Straub, Alexander Peltzer" \
      description="Docker image containing all software requirements for the nf-core/ampliseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-ampliseq-1.2.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-ampliseq-1.2.0dev > nf-core-ampliseq-1.2.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

## Required to build the container properly
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc
## Don't recache on each execution, do that once per build process
RUN qiime dev refresh-cache
