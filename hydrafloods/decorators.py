import ee
import functools


def carry_metadata(func):
    """Decorator function to set the properties of an image from computations to that of the input
        
    Function to decorate should take an ee.Image object and return an ee.Image object

    Example:
        @decorators.carry_metadata
        def ndvi(img):
            return img.normalizedDifference([b1,b2])
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
