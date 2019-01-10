# nfcore/ampliseq: Local Configuration

If running the pipeline in a local environment, we highly recommend using Singularity.

### Pipeline versions
The public Docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r 1.0`. This uses pipeline code and Docker image from this tagged version.


## Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

> NB: The "tag" at the end of this command corresponds to the pipeline version.
> Here, we're pulling the docker image for version 1.0 of the nfcore/ampliseq pipeline
> Make sure that this tag corresponds to the version of the pipeline that you're using

```bash
singularity pull --name nfcore-ampliseq-1.0.img docker://nfcore/ampliseq:1.0
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/nfcore-ampliseq -with-singularity /path/to/nfcore-ampliseq-1.0.img
```
