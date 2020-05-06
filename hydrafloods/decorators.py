import ee
import functools


def carryMetadata(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # expects first element of args is img
        # this assumption is true for 99.9% of fuctions used for ee.ImageCollection.map()
        result = ee.Image(func(*args,**kwargs))
        return ee.Image(result.copyProperties(args[0])\
            .set('system:time_start',args[0].get('system:time_start')))

    return wrapper
