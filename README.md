# hydra-floods

[![docs](https://github.com/Servir-Mekong/hydra-floods/workflows/docs/badge.svg)](https://servir-mekong.github.io/hydra-floods)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/hydrafloods.svg)](https://badge.fury.io/py/hydrafloods)
[![hydra-floods](https://circleci.com/gh/Servir-Mekong/hydra-floods/tree/master.svg?style=shield)](https://circleci.com/gh/Servir-Mekong/hydra-floods/tree/master)


## Introduction

The Hydrologic Remote Sensing Analysis for Floods (or HYDRAFloods) is an open source Python application for downloading, processing, and delivering surface water maps derived from remote sensing data. The bases behind the tool is to provide sensor agnostic approaches to produce surface water maps. Furthermore, there are workflows that leverage multiple remote sensing dataset in conjunction to provide daily surface water maps for flood application.

### Installation

The recommended way to get up and started using the `hydrafloods` packages is to install using `pip`:

```
pip install hydrafloods
```

`pip` should handle some of the basic dependencies such as the Earth Engine Python API that we need for the majority of the functionality. It is planned to add hydrafloods to the conda-forge channel but that is currently not completed.

To use the `hydrafloods` package successfully, Google Cloud and Earth Engine authentication is necessary. Tointialize the Google Cloud environment and authenticate using your credentials, run the following command:

```
gcloud init
```

To authenticate the Earth Engine Python API with your credentials, run the following:

```
earthengine authenticate
```

For more information on setup and installation of the `hydrafloods` package, please see the [Installation Docs](https://servir-mekong.github.io/hydra-floods/installation/).

## Example

To highlight a quick example of the `hydrafloods` API and simplicity to produce high-quality surface water maps we provide a quick example of mapping surface water using Sentinel-1 over the confluence of the Mekong and Tonle Sap rivers, which experiences frequent flooding.

```python
# import the hydrafloods and ee package
import hydrafloods as hf
import ee
ee.Initialize()

# specify start and end time as well as geographic region to process
start_time = "2019-10-05"
end_time =  "2019-10-06"
region = ee.Geometry.Rectangle([104, 11.5, 106, 12.5 ])

# get the Sentinel-1 collection
# the hf.dataset classes performs the spatial-temporal filtering for you
s1 = hf.datasets.Sentinel1(region, start_time, end_time)

# apply a water mapping function to the S1 dataset
# this applies the "Edge Otsu" algorithm from https://doi.org/10.3390/rs12152469
water_imgs = s1.apply_func(
    hf.thresholding.edge_otsu,
    initial_threshold=-14,
    edge_buffer=300
)

# take the mode from multiple images
# since this is just imagery from one day, it will simply mosaic the images
water_map = ee.Image(water_imgs.collection.mode())

# export the water map
hf.geeutils.export_image(
    water_map,
    region,
    "users/<YOUR_USERNAME>/water_map_example",
    scale=30,
)
```

_(This script is complete, it should run "as is")_

At the end of the script execution, there will be an [Earth Engine export task](https://developers.google.com/earth-engine/exporting#to-asset) running the process on the EE servers for use later in the EE platform. The resulting surface water image should look like the following figure. It should be noted that `hydrafloods` can scale quickly and easily by simply changing the start or end time and region to process, allowing for processing of surface water maps with minimal effort in terms of coding.

![Quick Start Results](docs/img/quick_start_results.png)

<!-- html code for figure caption -->
<span class="img_caption" style="display: block; text-align: center; font-size: 14px">
    __Figure 1.__ Sentinel-1 backscatter image (left) and resulting surface water map (right) from 2019-10-05 for a region in Cambodia as in the example.
</span>

Learn more about the package throughout the documentation such as [installation](https://servir-mekong.github.io/hydra-floods/installation/), the [algorithms available](https://servir-mekong.github.io/hydra-floods/algorithms/), or setting up the package to run operationally using the [CLI](https://servir-mekong.github.io/hydra-floods/cli/).

## Get in touch

- Report bugs, suggest features or view the source code on [GitHub](https://github.com/servir-mekong/hydra-floods).
- Contact us through a [Technical Assistance Request](https://servir.adpc.net/services/technical-assistance) and mention "hydrafloods"

## Contribute

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given. Please see the [Contributing Guidelines](https://github.com/servir-mekong/hydra-floods/blob/master/CONTRIBUTING.md) for details on where to contribute and how to get started.

## License

`hydrafloods` is available under the open source [GNU General Public License v3.0](https://github.com/Servir-Mekong/hydra-floods/blob/master/LICENSE).
