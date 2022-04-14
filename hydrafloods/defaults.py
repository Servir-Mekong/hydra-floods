s1 = {
    "edge_otsu": {
        "initial_threshold": -16,
        "edge_buffer": 300,
        "scale": 150,
        "band": "VV",
    },
    "bmax_otsu": {
        "initial_threshold": -16,
        "grid_size": 0.25,
        "scale": 150,
        "band": "VV",
    },
}

lc8 = {
    "edge_otsu": {
        "initial_threshold": 0.0,
        "edge_buffer": 300,
        "scale": 150,
        "invert": True,
        "band": "mndwi",
    },
    "bmax_otsu": {"initial_threshold": 0.0, "scale": 150, "invert": True},
}

fusion = {
    "edge_otsu": {
        "initial_threshold": 0.0,
        "edge_buffer": 300,
        "scale": 150,
        "invert": True,
    },
    "bmax_otsu": {"initial_threshold": 0.0, "scale": 150, "invert": True},
}
