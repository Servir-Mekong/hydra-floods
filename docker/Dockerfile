FROM continuumio/anaconda3

USER root

# install google cloud software development kit and support libraries
# Downloading gcloud package
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz

# Installing the package
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh

# Adding the package path to local
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin


# install supporting libraries
RUN apt-get update \
  && apt-get install -y \
            tar \
            git \
            curl \
            nano \
            wget \
            dialog \
            net-tools \
            build-essential\
            apt-transport-https \
            ca-certificates \
            gnupg \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# install Python dependencies
RUN conda install -c conda-forge -y\
            numpy \
            scipy \
            pandas \
            requests \
            yaml \
            xmltodict \
            gdal \
            shapely \
            pyproj \
            netCDF4 \
            xarray \
            pyresample \
            pillow \
            pandas \
            pyproj \
            tensorflow=2.1 \
            geopandas \
            oauth2client \
            earthengine-api \
            fire \
  && conda clean --all -f -y

# install the additional packages through pip/github
RUN pip install git+git://github.com/kmarkert/simple-cmr.git

# install the hydrafloods package
RUN pip install git+git://github.com/servir-mekong/hydra-floods.git

USER $NB_UID
