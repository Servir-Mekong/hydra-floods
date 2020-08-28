# Installation

`hydrafloods` itself is a pure Python package, but its dependencies are not. Furthermore, the package relies on Google Cloud and Google Earth Engine to stage and process data relying on specific software and even more importantly account authentication. There are two ways to install for use, one through a Docker Image and another via a manual installation.

## Using the Docker Image

The easiest way to get up and started using the `hydrafloods` packages is via a Docker Image. The Docker Image comes with pre-installed software and dependencies so you do not have to deal with mis-matching dependencies or sometimes difficult installations, such as GDAL.

To start you will need to have [Docker installed](https://docs.docker.com/get-docker/) on your system and running. You will need to pull the pre-built Docker Image for `hydrafloods` and start a new Container from the Image using the following command:

```sh
docker run -it \
  -v ~/<PROJECT-DIR>/:/mnt/<PROJECT-DIR> \
  --name hydrafloods_container kmarkert/hydrafloods
```

This command should be a one-time process to download the package and start the Container. Additionally, this command will mount a local directory (i.e. `~/<PROEJCT-DIR>`) for use within the Docker Container which allows you to edit files locally and use within the container. Be sure to change `<PROJECTD-DIR>` within the command to an exisiting local directory. Now the Docker Container is running for use!

Within the Docker Container the `hydrafloods` package and dependencies are pre-installed so all that is left is to [authenticate the cloud APIs](https://servir-mekong.github.io/hydra-floods/installation#cloud-authentication) then we will be ready to test and start processing.

If you have exited the Docker Container and want to start it again, use the following command:

```sh
docker start -ia hydrafloods_container
```

_This command to restart an existing Container is important especially after authenticating the cloud environment so that you do not have to go through the authentication process everytime you run the Docker container._ For more information on working with Docker Images/Containers using the CLI see the [Docker command line documentation](https://docs.docker.com/engine/reference/commandline/cli/).

## Manual installation

Another convient way to install the package and its dependencies is using [anaconda](https://www.anaconda.com/). It is recommend using the community maintained [conda-forge](https://conda-forge.github.io/) channel to handle dependencies. Furthermore, it is good practice to use a virtual environment within conda. To create a new environment, install dependencies, and activate the environment:

```sh
conda create -n hydra -c conda-forge python=3.7 \
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
  scikit-learn \
  pyresample \
  geopandas \
  earthengine-api \
  gcsfs \
  fire -y
conda activate hydra
```

Finally, we need to install the `hydrafloods` package and one last dependency via `pip`:

```sh
pip install simplecmr hydrafloods
```

You will now also need to install the [Google Cloud SDK](https://cloud.google.com/sdk/docs/downloads-versioned-archives) to interface to with the Google cloud. Follow the directions provided by the website.

Once all of the source code and dependencies has been installed successfully, you will need to [authenticate the cloud APIs](https://servir-mekong.github.io/hydra-floods/installation#cloud-authentication)

## Cloud authentication

After successful installation of the package and dependencies we will need to authenticate our local installation (or within the Docker Container) to interface with Google Cloud and Earth Engine. Running these command will prompt you through the authentication process using a web browser.

<span style="color:red">
**Warning:** *Make sure you initialize the `earthengine` and `gcloud` APIs with Google accounts that have permissions to read and write to Google Cloud Storage and Google Earth Engine assets.*
</span>

To intialize the Google Cloud environment and authenticate using your credentials, run the following command:

```sh
gcloud init
```

To authenticate the Earth Engine Python API with your credentials, run the following:

```sh
earthengine authenticate
```

Now we are ready to test our installation!

## Testing installation

🚧 Coming soon! 🚧
