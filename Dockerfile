FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/ampliseq pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-ampliseq-1.0.0/bin:$PATH
## Required to build the container properly
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc
## Don't recache on each execution, do that once per build process
RUN qiime dev refresh-cache

