"""
This suppose to give a perspective to contact classes
Testing collsion, giving our contact normal | penetration | and other trivial things ig
"""

import pygame as pg
from pygame.locals import *

import numpy as np


default_font = pg.font.SysFont("Arial", 20)
# default_font = pg.font.SysFont("Consolas", 15)
# default_font = pg.font.SysFont("Courier New", 10)

class InspectingTexts:
    def __init__(self, color: tuple[int, int, int] = (255, 255, 255), universal_font: pg.font.Font = default_font, anti_alias = True):
        self.font = universal_font
        self.color = color
        self.anti_alias = anti_alias

        self.img_queue: list[pg.surface.Surface] = []
        self.pos_queue: list[np.array] = []
        self.center_queue: list[int] = []
    
    @staticmethod
    def ManualLog(message: str, surface: pg.surface.Surface, position, color: tuple[int, int, int] = (255, 255, 255), font: pg.font.Font = default_font, center = False):
        if type(message) != str:
            message = str(message)
        font_image = font.render(message, True, color).convert_alpha()
        x, y = position
        if center == True:
            x -= font_image.get_width()  / 2
            y -= font_image.get_height() / 2 
        surface.blit(font_image, (x, y))

    def Add(self, message: str, position: np.array, is_center = False):
        if type(message) != str:
            message = str(message)
        self.img_queue.append(default_font.render(message, self.anti_alias, self.color))
        self.pos_queue.append(position)
        self.center_queue.append(is_center)
        
    def Clear(self):
        self.img_queue = []
        self.pos_queue = []
        self.center_queue = []
    
    def Draw(self, surface: pg.surface.Surface):
        length = len(self.img_queue)
        for i in range(length):
            image = self.img_queue[i]
            x, y = self.pos_queue[i]
            if self.center_queue[i]:
                x -= image.get_width() / 2
                y -= image.get_height() / 2
            surface.blit(image, (x, y))
            


class InspectingPoints:
    """
    Book a point, go through the engine interval, storing the hired points. 
    Eventually, would be queued to show up on the "bug screen",
    The function is created wholeheartly, exclusively for debug purposes in position, maybe vectors in the future
    """
    def __init__(self, color: tuple[int, int, int] = (255, 0, 0), radius: int = 5):
        self.color = color
        self.radius = radius
        self.pos_queue: list[np.array] = []
        self.color_queue: list[tuple[int, int, int]] = []
    
    def Clear(self):
        self.pos_queue.clear()
        self.color_queue.clear()
    
    def Add(self, point: np.array, color: tuple[int, int, int] = None):
        if color == None: 
            color = self.color
        self.color_queue.append(color)
        self.pos_queue.append(point)
    
    def Draw(self, surface):
        for i in range(len(self.pos_queue)):
            pg.draw.circle(surface, self.color_queue[i], self.pos_queue[i].tolist(), self.radius)

sqrt_3 = np.sqrt(3)

class InspectingLines:
    def __init__(self, width: int = 1):
        self.queue_p1: list[np.ndarray] = []
        self.queue_p2: list[np.ndarray] = []
        self.color_queue: list[np.array] = []
        self.width = width
    
    def Add(self, p1: np.ndarray, p2: np.ndarray, is_vector, color: tuple[int, int, int] = (255, 255, 255)):
        # if is_vector:
        #     p2[0] = p1[0] + p2[0]
        #     p2[1] = p1[1] + p2[1]
        self.queue_p1.append(p1)
        if is_vector:
            self.queue_p2.append(p1 + p2)
        else: 
            self.queue_p2.append(p2)
        self.color_queue.append(color)
    
    def Clear(self):
        self.queue_p1.clear()
        self.queue_p2.clear()
        self.color_queue.clear()

    def Draw(self, surface: pg.surface.Surface):
        length = len(self.queue_p1)
        for i in range(length):
            pg.draw.line(surface, self.color_queue[i], self.queue_p1[i], self.queue_p2[i], self.width)
        

class InspectingVectors:
    """
    THIS FUNCTION IS RATHER EXPENSIVE TO PULL OUT, TRIPLE THE COST COMPARED TO THE OTHER TWO METHODS.
    USE IT WISELY.
    """
    default_pointer_vertices = [
        np.array([-sqrt_3 / 6, -0.5], dtype=np.float64),
        np.array([-sqrt_3 / 6,  0.5], dtype=np.float64),
        np.array([sqrt_3 / 3,  0], dtype=np.float64),
    ]

    def __init__(self, scale, thickness: int = 2, color: tuple[int, int, int] = (118, 185, 0)):
        self.scale = scale
        self.color = color
        self.thickness = thickness
        self.position_queue = []
        self.vector_queue = []

    def Add(self, vector: np.array, position: np.array):
        if vector @ vector == 0:
            return
        self.vector_queue.append(vector)
        self.position_queue.append(position)
    
    def CreatePointerVertices(self, vector: np.array, position: np.array) -> np.array:
        vector_mag = np.linalg.norm(vector)
        tangent = vector / vector_mag  # Normalize the vector

        pnt_vertices = np.zeros((3, 2), dtype = np.float64)
        for i in range(len(pnt_vertices)):
            pnt_vertices[i] = self.default_pointer_vertices[i] * self.scale
            pnt_vertices[i][0] += vector_mag
        
        rotated_translated_vertices = np.array([self.Rotate(v, tangent[0], tangent[1]) + position for v in pnt_vertices])

        return rotated_translated_vertices

    @staticmethod
    def Rotate(vertex, cos: float, sin: float) -> np.array:
        """
        Rotate a vertex using the angle's cosine and sine values.
        NOTE: Return np.array([cos * x - sin * y, sin * x + cos * y], dtype=np.float64)
        """
        x, y = vertex
        vertex[0] = cos * x - sin * y
        vertex[1] = sin * x + cos * y
        return vertex

    def Draw(self, surface: pg.surface.Surface):
        for position, vector in zip(self.position_queue, self.vector_queue):
            polygon_vertices = self.CreatePointerVertices(vector, position)
            destination = position + vector
            pg.draw.line(surface, self.color, position.tolist(), destination.tolist(), self.thickness)
            pg.draw.polygon(surface, self.color, polygon_vertices)

    def Clear(self):
        self.position_queue.clear()
        self.vector_queue.clear()

class WarningLogs:
    """
    supposed to work the same as window caution log windows. 
    USAGE: It prints the warning log once in the start of the loop, then stays in the function as an O(1) check statement
    
    EXAMPLE: 
            warning_logs = WarningLogs()
            while (true)
            {
                ...
                warning_logs.Warning("Taiwan is an indepent country, why not :)")
                ...
                warning_logs.reset()
            }

    not coreect ! (rework needed)
    """

    def __init__(self):
        self.i: int = 0 # should be reset every loop
        self.total_i_passed: int = 0
    
    def Warning(self, message: str = "Taiwan is an independent country"):
        self.i += 1
        # if the log have not been found before, in the previous loops
        if self.total_i_passed < self.i: 
            print(message)
            self.total_i_passed += 1
    
    def Reset(self):
        self.i = 0
        

    

    
    

inspecting_points = InspectingPoints(radius = 4)
inspecting_texts = InspectingTexts()
inspecting_vectors = InspectingVectors(scale = 5, thickness = 2, color = (255, 255, 255))
inspecting_lines = InspectingLines(width = 1)
warning_logs = WarningLogs()

