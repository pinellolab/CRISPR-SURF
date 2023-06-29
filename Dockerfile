############################################################
# Dockerfile to build STREAM & webapp
############################################################

# Set the base image to anaconda python 2.7
FROM continuumio/anaconda2

# File Author / Maintainer
MAINTAINER Jonathan Y. Hsu

ENV SHELL bash

RUN conda install r-base
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

#Add build tools
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install build-essential zlib1g-dev -y

#add Python dependencies
RUN apt-get install libreadline-dev -y
RUN pip install CVXCanon==0.1.0 # ADDED 1/18/2023 - prevent incompatibility between cvxpy and CVXCanon
RUN pip install cvxpy==0.4.11
RUN apt-get install unzip libxml2 libxml2-dev -y

#website dependencies
RUN pip install Flask-Compress==1.4.0
RUN pip install dash==0.21.0  # The core dash backend
RUN pip install dash-renderer==0.12.1  # The dash front-end
RUN pip install dash-html-components==0.8.0  # HTML components
RUN pip install dash-core-components==0.22.1  # Supercharged components
RUN pip install dash-core-components==0.21.0rc1 # Tabs pre-release
RUN pip install plotly --upgrade  # Plotly graphing library used in examples
RUN pip install dash-table-experiments
RUN pip install gunicorn

# install zips
RUN apt-get update && apt-get install zip zlib1g liblzo2-dev -y

#new dependencies
RUN pip install bx-python
RUN git clone https://github.com/lucapinello/bioutilities.git && cd bioutilities/ && python setup.py install

# create environment
COPY SURF /SURF
#COPY /SURF/static/CRISPR-SURF.css /opt/conda/lib/python2.7/site-packages/dash_core_components/
#COPY /SURF/static/Loading-State.css /opt/conda/lib/python2.7/site-packages/dash_core_components/
#COPY /SURF/static/jquery-3.3.1.min.js /opt/conda/lib/python2.7/site-packages/dash_core_components/
RUN mkdir /tmp/UPLOADS_FOLDER
RUN mkdir /tmp/RESULTS_FOLDER

# Reroute to enable the STREAM CLI and STREAM webapp
WORKDIR /SURF
EXPOSE 9993
#CMD ["bash", "start_server_docker.sh"]
ENTRYPOINT ["/opt/conda/bin/python", "/SURF/crisprsurf_router.py"]
