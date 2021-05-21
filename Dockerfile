FROM continuumio/miniconda3
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate env" > /root/.bashrc
RUN apt-get update && apt-get -y install gcc
ENV PATH /opt/forcoast:$PATH
ENV PATH /opt/forcoast/preprocessing:$PATH
ENV PATH /opt/forcoast/processing:$PATH
ENV PATH /opt/forcoast/postprocessing:$PATH

COPY PreProcessing /opt/forcoast/preprocessing
COPY Processing /opt/forcoast/processing
COPY PostProcessing /opt/forcoast/postprocessing
COPY parcels /opt/conda/envs/env/lib/python3.6/site-packages/parcels
COPY coastlines /opt/forcoast/coastlines
COPY run.sh /opt/forcoast/run.sh

RUN chmod a+x /opt/forcoast/run.sh

ENTRYPOINT ["bash", "-c"]
# ENTRYPOINT ["python", "/opt/forcoast/preprocessing/forcoast_download.py"]
# ENTRYPOINT ["/opt/forcoast/run.sh"]
# CMD ["python /root/glossis.py"]
