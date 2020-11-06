import ee
import functools


def carry_metadata(func):
    """Decorator function to set the properties of an image from computations to that of the input
    Function to decorate should take an ee.Image object and return an ee.Image object

    args:
        func (object): function object to wrap. Expects that an element within args is of type 
            ee.Image and will use first ee.Image to carry metadata

    Example:
        ```python
        @decorators.carry_metadata
        def ndvi(img):
            return img.normalizedDifference([b1,b2])
        ```
        Returned image(s) will have all of the same metadata properties as the input including `system:time_start`
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # expects an element within args is img
        # will set the metadata to first ee.Image instance
        # this assumption is true for 99% of fuctions used for ee.ImageCollection.map()
        result = ee.Image(func(*args, **kwargs))
        img = [i for i in args if isinstance(i, ee.Image)][0]
        return ee.Image(
            result.copyProperties(img).set(
                "system:time_start", img.get("system:time_start")
            )
        )

    return wrapper
