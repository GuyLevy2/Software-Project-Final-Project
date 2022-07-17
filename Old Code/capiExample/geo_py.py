import numpy as np

def geo(z:float, n: int) -> float:
    return sum(z+n)

def geo_np(z: float, n: int):
    return np.sum(n**2+z**2)