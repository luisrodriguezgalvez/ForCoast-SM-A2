FROM continuumio/miniconda3

WORKDIR /usr/src/app

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate env" > /root/.bashrc
RUN apt-get update && apt-get -y install gcc
ENV PATH /usr/src/app:$PATH
ENV PATH /usr/src/app/preprocessing:$PATH
ENV PATH /usr/src/app/processing:$PATH
ENV PATH /usr/src/app/postprocessing:$PATH
ENV PATH /usr/src/app/bulletinscript:$PATH

COPY PreProcessing /usr/src/app/preprocessing
COPY Processing /usr/src/app/processing
COPY PostProcessing /usr/src/app/postprocessing
COPY BulletinScript /usr/src/app/bulletinscript
COPY parcels /opt/conda/envs/env/lib/python3.6/site-packages/parcels
COPY coastlines /usr/src/app/coastlines
COPY run.sh /usr/src/app/run.sh
COPY data /usr/src/app/data

RUN chmod a+x /usr/src/app/run.sh

# ENTRYPOINT ["bash", "-c"]
ENTRYPOINT ["/usr/src/app/run.sh"]

