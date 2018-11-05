FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/rrna-ampliseq pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN mkdir -p ~/.config/matplotlib
RUN echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
ENV PATH /opt/conda/envs/nf-core-rrna-ampliseq-1.0dev/bin:$PATH
