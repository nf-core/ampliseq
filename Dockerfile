FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/rrna-ampliseq pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rrna-ampliseq-1.0dev/bin:$PATH
RUN qiime dev refresh-cache

