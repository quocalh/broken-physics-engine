import numpy as np

WIDTH = 400
HEIGHT = 400

FPS = 100000
FPS = 60

GRAVITATIONAL_ACCELERATION = np.array([0, 780], dtype = np.float64)

def CrossProductVector(vector1: np.array, vector2: np.array) -> float:
    """
    Determinant references
    """
    return (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])

def CrossProductScalar(vector: np.array, z_component: float):
    """
    the cross prod of (vector[0], vector[1], 0) x (0, 0, z_component)
    the z always appear to be 0 so we cut it, making it a 2d vec | it is a 2d engine anyway
    """
    return np.array([vector[1] * z_component, -vector[0] * z_component], dtype = np.float64)
