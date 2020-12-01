FROM continuumio/miniconda3:latest

USER root

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
            build-essential \
            libssl-dev \
            libffi-dev \
            vim \
            curl \
  && apt-get purge -y build-essential \
            libssl-dev \
            libffi-dev \
            dpkg-dev \
            fakeroot \
            libfakeroot:amd64 \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# install google cloud software development kit and support libraries
# Downloading gcloud package
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz

# Installing the package
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh

# remove the gcp-sdk install tarball
RUN rm /tmp/google-cloud-sdk.tar.gz

# Adding the package path to local
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin

# install Python dependencies
RUN conda install -c conda-forge -y\
            oauth2client \
            "earthengine-api<=0.1.232" \
            gcsfs \
            fire \
  && conda clean --all -f -y

# install the additional packages through pip/github
# install the hydrafloods package
RUN pip install --no-cache-dir hydrafloods

USER $NB_UID
