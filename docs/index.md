# Welcome to the HYDRAFloods Documentation

The Hydrologic Remote Sensing Analysis for Floods (or HYDRAFloods) is an open source Python application for downloading, processing, and delivering surface water maps derived from remote sensing data. The bases behind the tool is to provide sensor agnostic approaches to produce surface water maps. Furthermore, there are applications that leverage multiple remote sensing dataset in conjunction to provide daily surface water maps for flood application. 
<!-- Mention something about EE. -->

The HYDRAFloods application is built using [Google Earth Engine](https://earthengine.google.com/) and [Google Cloud Platform](https://cloud.google.com/) to leverage cloud computing for large-scale computations and handling high data volume outputs. 


## Quick Start

To highlight an example of the `hydrafloods` API and simplicity to produce high-quality surface water maps we provide a quick example of mapping surface water using Sentinel-1 over the confluence of the Mekong and Tonle Sap rivers, which experiences frequent flooding.

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
# this dataset classes does the spatial-temporal filtering for you 
s1 = hf.datasets.Sentinel1(region, start_time end_time)

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
    "users/<YOUR_USERNAME>/water_map_example,
    scale=30,
)
```

At the end of processing, there will be an [Earth Engine export task](https://developers.google.com/earth-engine/exporting#to-asset) that running the process on the EE servers for use later in the EE platform.

![Quick Start Results](img/quick_start_results.png)
Figure 1. Sentinel-1 backscatter image (left) and resulting surface water map (right) from the example.

Learn more throughout the documentation such as [Installation](https://servir-mekong.github.io/hydra-floods/installation/), the [Algorithms available](https://servir-mekong.github.io/hydra-floods/algorithms/), or setting up the package to run operationally using the [CLI](https://servir-mekong.github.io/hydra-floods/cli/).


## Get in touch
* Report bugs, suggest features or view the source code on [GitHub](https://github.com/servir-mekong/hydra-floods).
* Contact us through a [Technical Assistance Request](https://servir.adpc.net/services/technical-assistance) and mention "hydrafloods"

## Contribute

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given. Please see the [Contributing Guidelines](https://github.com/servir-mekong/hydra-floods/blob/master/CONTRIBUTING.md) for details on where to contribute and how to get started.


## License
`hydrafloods` is available under the open source [GNU General Public License v3.0](https://github.com/Servir-Mekong/hydra-floods/blob/master/LICENSE).
