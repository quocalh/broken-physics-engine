import numpy as np
import pygame as pg
import math

pg.init()

from components.debug_kit import *
from components.settings import *

default_array  = np.zeros(2, dtype = np.float64)

# def CrossProductFloat(a: np.ndarray, b: np.ndarray) -> float:
#     return a[0] * b[1] - a[1] * b[0]

class PolygonsAttributes:
    def __init__(self, vertices: np.ndarray):
        self.local_vertices: list[np.ndarray] = vertices 

        
        self.local_faces: list[np.array] = [np.zeros(2, dtype = np.float64) for _ in range(len(self.local_vertices))]
        self.GetLocalFaces()

        self.is_clockwise: bool = None
        self.SetClockwise()
        self.sign = 1 - 2 * int(not self.is_clockwise)

        self.local_normals: list[np.ndarray] = [np.zeros(2, dtype = np.float64) for _ in range(len(self.local_vertices))]
        self.GetLocalNormals()

        self.radian : float = 0
        self.previous_radian : float = 0 # not sure if i need it
        self.delta_radian : float = 0
        
        self.world_vertices : list[np.ndarray] = [np.zeros(2, dtype = np.float64) for _ in range(len(self.local_vertices))]
        self.world_normals : list[np.ndarray] = [np.zeros(2, dtype = np.float64) for _ in range(len(self.local_vertices))]
        self.world_faces : list[np.ndarray] = [np.zeros(2, dtype = np.float64) for _ in range(len(self.local_vertices))]
        
        # self.world_prime_axis_index: int = -1 # id of axis of failure


    def GetLocalFaces(self):
        len_local_vertices = len(self.local_vertices)
        for i in range(len_local_vertices):
            i_vertex: int = i
            next_i_vertex: int =  0 if i + 1 == len_local_vertices else i + 1
            faces = self.local_vertices[next_i_vertex] - self.local_vertices[i_vertex]
            self.local_faces[i][0] = faces[0]
            self.local_faces[i][1] = faces[1]

    def GetLocalNormals(self):
        for i in range(len(self.local_vertices)):
            length = np.linalg.norm(self.local_faces[i])
            self.local_normals[i][0] = self.local_faces[i][1]
            self.local_normals[i][1] = -self.local_faces[i][0]
            self.local_normals[i] *= self.sign / length
    
    def DebugNormalVector(self, length: int = 50):
        len_world_vertices = len(self.world_vertices)
        for i in range(len_world_vertices):
            vertex_i = i
            next_vertex_i = 0 if i + 1 == len_world_vertices else i + 1
            normal = self.world_normals[i] * length
            mid_position = self.world_vertices[vertex_i] / 2 + self.world_vertices[next_vertex_i] / 2
            inspecting_vectors.Add(normal, mid_position)
    
            

    def GetAABB(self):
        """
        THE FUNC IS ALREADY WRITTEN IN THE BODY CLASS
        """
        pass
    
    def isShape(self):
        return len(self.local_vertices) > 2 # would only call once every init

    @staticmethod
    def RotateMatrixMultiplication2D(vector: np.array, radian: float):
        cos = np.cos(radian)
        sin = np.sin(radian)
        x = vector[0] * cos - vector[1] * sin
        y = vector[0] * sin + vector[1] * cos
        return np.array([x, y], dtype = np.float64)

    def Rotating(self):
        for i in range(len(self.local_vertices)):
            world_vertices = self.RotateMatrixMultiplication2D(self.local_vertices[i], self.radian)
            world_normals = self.RotateMatrixMultiplication2D(self.local_normals[i], self.radian)
            world_faces = self.RotateMatrixMultiplication2D(self.local_faces[i], self.radian)

            self.world_vertices[i][0] = world_vertices[0]
            self.world_vertices[i][1] = world_vertices[1]

            self.world_normals[i][0] = world_normals[0]
            self.world_normals[i][1] = world_normals[1]

            self.world_faces[i][0] = world_faces[0]
            self.world_faces[i][1] = world_faces[1]

    


    def Translating(self, translation_value: np.ndarray):
        for i in range(len(self.world_vertices)):
            self.world_vertices[i][0] = self.world_vertices[i][0] + translation_value[0]
            self.world_vertices[i][1] = self.world_vertices[i][1] + translation_value[1]
        

    def SetRadian(self, new_radian: float):
        self.radian = new_radian

    def SetClockwise(self):
        a = - self.local_faces[0]
        b = self.local_faces[1]
        self.is_clockwise = CrossProductVector(a, b) < 0




class Transform:
    def __init__(self):
        self.orientational_matrix = np.array([[1, 0], [0, 1]], dtype = np.float64)
            
    def Rotate(self, rad: float):
        self.orientational_matrix[0][0] = math.cos(rad)
        self.orientational_matrix[0][1] = math.sin(rad)
        self.orientational_matrix[1][0] = -math.sin(rad)
        self.orientational_matrix[1][1] = math.cos(rad)

class MassData:
    def __init__(self, iMass: float, inverse_inertia: float =  1/10):
        self.iMass: float = iMass
        if iMass != 0:
            self.mass: float = 1 / iMass
        else: 
            self.mass = None
        self.inverse_inertia: float = inverse_inertia
        if inverse_inertia != 0:
            self.inertia: float = 1 / inverse_inertia
    
class Materials:
    def __init__(self, restitution: float, denstiy: float = 1):
        self.restitution: float = restitution 
        self.density: float = denstiy

class FrictionCoefficient:
    def __init__(self, dynamic_friction: float = 0.05, static_friction: float = 0.01):
        self.static_friction: float = static_friction
        self.dynamic_friction: float = dynamic_friction