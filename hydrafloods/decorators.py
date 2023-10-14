import ee
import functools


def keep_attrs(func):
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
        # this assumption is true for 99% of functions used for ee.ImageCollection.map()
        result = ee.Image(func(*args, **kwargs))
        img = [i for i in args if isinstance(i, ee.Image)][0]
        return img.select().addBands(result)

    return wrapper


def keep_names(func):
    """Decorator function to set the band names of an image from computations to that of the input
    Function to decorate should take an ee.Image object and return an ee.Image object

    args:
        func (object): function object to wrap. Expects that an element within args is of type
            ee.Image and will use first ee.Image to carry metadata

    Example:
        ```python
        @decorators.keep_names
        def my_computation(img):
            const = ee.Image.constant(2)
            return const.multiply(img)
        ```
        Returned image(s) will have the same band names as the input.
        This should only be used if the input and output images will have the same number of bands!
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # expects an element within args is img
        # will set the metadata to first ee.Image instance
        # this assumption is true for 99% of functions used for ee.ImageCollection.map()
        result = ee.Image(func(*args, **kwargs))
        img = [i for i in args if isinstance(i, ee.Image)][0]
        return result.rename(img.bandNames())

    return wrapper
