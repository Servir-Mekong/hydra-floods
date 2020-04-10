import ee
import functools


def carryMetadata(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # expects first element of args is img
        result = func(*args,**kwargs)
        return ee.Image(result.copyProperties(args[0])\
            .set('system:time_start',args[0].get('system:time_start')))

    return wrapper
