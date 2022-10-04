FROM continuumio/miniconda3

WORKDIR /usr/src/app

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate forcoastA2" > /root/.bashrc
RUN apt-get update && apt-get -y install gcc && apt-get -y install imagemagick ffmpeg libsm6 libxext6
RUN pip install wget==3.2 pyyaml==6.0 telepot==12.7 argparse==1.4.0 requests==2.28.1

ENV PATH /usr/src/app:$PATH
ENV PATH /usr/src/app/PreProcessing:$PATH
ENV PATH /usr/src/app/Processing:$PATH
ENV PATH /usr/src/app/PostProcessing:$PATH
ENV PATH /usr/src/app/BulletinScript:$PATH

COPY PreProcessing /usr/src/app/PreProcessing
COPY Processing /usr/src/app/Processing
COPY PostProcessing /usr/src/app/PostProcessing
COPY BulletinScript /usr/src/app/BulletinScript
#COPY parcels /opt/conda/envs/env/lib/python3.6/site-packages/parcels
COPY coastlines /usr/src/app/coastlines
COPY run_py.sh /usr/src/app/run_py.sh
COPY data /usr/src/app/data
COPY usr /usr/src/app/usr
COPY Telegram /usr/src/app/Telegram

RUN chmod a+x /usr/src/app/run_py.sh
SHELL ["conda", "run", "-n", "forcoastA2", "/bin/bash", "-c"]

# ENTRYPOINT ["bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "forcoastA2", "/usr/src/app/run_py.sh"]

