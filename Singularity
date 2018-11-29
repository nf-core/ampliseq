From:nfcore/base
Bootstrap:docker

%labels
    DESCRIPTION Singularity image containing all requirements for nf-core/ampliseq pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-ampliseq-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
    /opt/conda/bin/conda clean -a
