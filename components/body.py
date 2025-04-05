import numpy as np
from numpy.core.multiarray import array as array
import pygame as pg
import math

from components.body_components import *
from components.body_components import MassData, Materials

pg.init()

    

    
"""
██████╗░██╗░██████╗░██╗██████╗░  ██████╗░░█████╗░██████╗░██╗░░░██╗
██╔══██╗██║██╔════╝░██║██╔══██╗  ██╔══██╗██╔══██╗██╔══██╗╚██╗░██╔╝
██████╔╝██║██║░░██╗░██║██║░░██║  ██████╦╝██║░░██║██║░░██║░╚████╔╝░
██╔══██╗██║██║░░╚██╗██║██║░░██║  ██╔══██╗██║░░██║██║░░██║░░╚██╔╝░░
██║░░██║██║╚██████╔╝██║██████╔╝  ██████╦╝╚█████╔╝██████╔╝░░░██║░░░
╚═╝░░╚═╝╚═╝░╚═════╝░╚═╝╚═════╝░  ╚═════╝░░╚════╝░╚═════╝░░░░╚═╝░░░
"""

Materials_Density_StringDict = {
    "Rock":         0.6,
    "Wood":         0.3,
    "Metal":        1.2,
    "BouncyBall":   0.3,
    "SuperBall":    0.3,
    "Pillow":       0.1,
    "Static":       0.0,
} 

Materials_Restitution_StringDict = {
    "Rock":         0.1,
    "Wood":         0.2,
    "Metal":        0.05,
    "BouncyBall":   0.8,
    "SuperBall":    0.95,
    "Pillow":       0.2,
    "Static":       0.4,
} 

class AxisAlignedBoundingBox:
    def __init__(self, x: float, y: float, w: int, h: int):
        self.x = x
        self.y = y
        self.w = w
        self.h = h

class Body:
    def __init__(self, position: np.array, 
                 mass_data: MassData = MassData(1/200), 
                 materials: Materials = Materials(0.1), 
                 color: tuple[int, int, int] = (0, 0, 255),
                ):
        self.mass_data = mass_data    
        self.color = color

        self.position = position
        self.velocity = np.copy(default_array)
        self.acceleration = np.copy(default_array)
        self.force = np.copy(default_array)

        self.material: Materials =  materials
        # self.friction_static: float = friction_static
        
        # they only exist for aabb debugs
        self.AABBwidth: int = None
        self.AABBheight: int = None
        self.AABBposition: np.array = None
        self.rect: pg.Rect = pg.Rect(0, 0, 0, 0)

        self.AABB: AxisAlignedBoundingBox = None 

        # is rest
        self.is_rest = False

        # ultimate forces
        self.orientation: float
        self.angular_velocity: float
        self.torque: float

    def ResetVelocity(self):
        self.velocity[0] = 0
        self.velocity[1] = 0
    def ResetAcceleration(self):
        self.acceleration[0] = 0
        self.acceleration[1] = 0
    def ResetForce(self):
        self.force[0] = 0
        self.force[1] = 0

    def AccumulatingForce(self, value: np.array):
        self.force += value # D'Larmbert Principal   
    def AccumulatingAcceleration(self, value: np.array):
        self.acceleration += value
    def AccumulatingVelocity(self, value: np.array):
        self.velocity += value
    def AccumulatingPosition(self, value: np.array):
        self.position += value      
    
    def UpdateAcceleration(self):
        accelerationGenerated = self.force * self.mass_data.iMass
        self.AccumulatingAcceleration(accelerationGenerated)
    def UpdateVelocity(self, dt: float):
        self.velocity += self.acceleration * dt
    def UpdatePosition(self, dt: float):
        self.position += self.velocity * dt + 0.5 * self.acceleration * dt * dt

    def ClearKinetic(self):
        self.ResetAcceleration()
        self.ResetForce()

    def Integrate(self, dt: float):
        self.UpdateAcceleration()
        self.UpdateVelocity(dt)
        self.UpdatePosition(dt)
    
    def Draw(self, surface: pg.surface.Surface):
        print("drawing in process ... ")
    
    def DrawAABB(self, surface: pg.surface.Surface):
        assert False, f"you have not created the self.DrawwAABB function for the children class: {type(self)}"

    def ComputeAABB(self):
        assert False, f"you have not created the self.ComputeAABB function for the children class!{type(self)} "
        print("creating an AABB broadphase collision check")
        raise "you havent not created the AABB broad check for the ObjectType"
    def CreateAABB(self):
        assert False, f"you have not created the self.CreateAABB function for the children class!{type(self)} "

    def get_AABBheight(self):
        print("hello height")
    def get_AABBwidth(self):
        print("hello width")
    def get_AABBposition(self):
        print("hello position")

def is_collided_AABB_vs_AABB_body(object1: Body, object2: Body):
        body1pos: np.array = object1.get_AABBposition()
        body1w: int = object1.get_AABBwidth()
        body1h: int = object1.get_AABBheight()
        
        body2pos: np.array = object2.get_AABBposition()
        body2w: int = object2.get_AABBwidth()
        body2h: int = object2.get_AABBheight()
        # the if statement is more intuitive by isolating the checking axis, only try to illustrate in your mind 
        if (body1pos[0]                                         < body2pos[0] + body2w  and
            body1pos[0] + body1w                                > body2pos[0]           and
            body1pos[1]                                         < body2pos[1] + body2h  and
            body1pos[1] + body1h                                > body2pos[1]
            ):
            return True
        return False

def is_collided_AABB_vs_AABB_AABB(aabb1: AxisAlignedBoundingBox, aabb2: AxisAlignedBoundingBox):
        aabb1x: np.array = aabb1.x
        aabb1y: np.array = aabb1.y
        aabb1w: int = aabb1.w
        aabb1h: int = aabb1.h

        aabb2x: np.array = aabb2.x
        aabb2y: np.array = aabb2.y
        aabb2w: int = aabb2.w
        aabb2h: int = aabb2.h
        # the if statement is more intuitive by isolating the checking axis, only try to illustrate in your mind 
        if (aabb1x                                          < aabb2x + aabb2w  and
            aabb1x + aabb1w                                 > aabb2x           and
            aabb1y                                          < aabb2y + aabb2h  and
            aabb1y + aabb1h                                 > aabb2y
            ):
            return True
        return False


class BodyRect(Body):
    def __init__(self, position: np.array, width: int, height: int,
                 materials: Materials = Materials(0.1),
                 mass_data: MassData = MassData(1 / 20),
                 color: tuple = (0, 0, 255)):
        super().__init__(position, mass_data, materials, color)
        self.width = width
        self.height = height
        self.rect = pg.Rect(int(position[0] - width / 2), int(position[1] - height / 2), width, height)

    def get_AABBheight(self):
        return self.height

    def get_AABBwidth(self):
        return self.width

    def get_AABBposition(self):
        return self.position - np.array([self.width / 2, self.height / 2], dtype=np.float64)

    def ComputeAABB(self):
        self.AABBheight = self.get_AABBheight()
        self.AABBwidth = self.get_AABBwidth()
        self.AABBposition = self.get_AABBposition()
    
    def CreateAABB(self):
        x = self.position[0] - self.width / 2
        y = self.position[1] - self.height / 2
        self.AABB = AxisAlignedBoundingBox(x, y, self.width, self.height)

    def Draw(self, surface: pg.surface.Surface):
        self.rect.left = int(self.position[0] - self.width / 2)
        self.rect.top = int(self.position[1] - self.height / 2)
        pg.draw.rect(surface, self.color, self.rect)
    
    def DrawAABB(self, surface: pg.surface.Surface):
        self.rect.left = int(self.position[0] - self.width / 2)
        self.rect.top = int(self.position[1] - self.height / 2)
        pg.draw.rect(surface, self.color, self.rect, 10)
        


class BodyCircle(Body):
    def __init__(self, position: np.array, 
                 radius: int, 
                 materials: Materials = Materials(0.1), 
                 mass_data: MassData = MassData(1 / 200), 
                 color: tuple = (0, 0, 255)):
        super().__init__(position, mass_data, materials, color)
        self.radius = radius

    def get_AABBwidth(self):
        return self.radius * 2
    def get_AABBheight(self):
        return self.radius * 2
    def get_AABBposition(self):
        return self.position - np.array([self.radius, self.radius], dtype = np.float64)

    def ComputeAABB(self):
        self.AABBheight = self.get_AABBheight()
        self.AABBwidth = self.get_AABBwidth()
        self.AABBposition = self.get_AABBposition()
    
    def CreateAABB(self):
        x = self.position[0] - self.radius
        y = self.position[1] - self.radius
        self.AABB = AxisAlignedBoundingBox(x, y, self.radius * 2, self.radius * 2)
        
    def Draw(self, surface: pg.Surface):
        pg.draw.circle(surface, self.color, self.position.tolist(), self.radius)

class BodyPolygon(Body):
    def __init__(self, position: np.array, 
                 polygon: PolygonsAttributes, 
                 mass_data: MassData = MassData(1 / 200), 
                 materials: Materials = Materials(0.1), 
                 color: tuple[int, int, int] = (0, 0, 255),
                 is_fixed_world_pivot = False
                 ):
        """THE ONLY CONSTRAINT FOR THE VERTICES: EACH NEXT VERTEX HAS TO BE ADJANCED TO EACH OTHER, THE LAST ELEMENT WOULD BE USED TO LINK BACK TO THE START IF NEEDED"""
        super().__init__(position, mass_data, materials, color)

        self.polygon: PolygonsAttributes = polygon
        self.angular_velocity: float = 0
        self.polygon.Rotating()
        self.polygon.Translating(self.position)

        self.world_pivot: np.ndarray = np.zeros(2, dtype = np.float64)

        # reset the pivot to the center of mass 
        self.world_pivot[0] = self.position[0]
        self.world_pivot[1] = self.position[1]
        self.is_fixed_world_pivot: bool = is_fixed_world_pivot


    def Integrate(self, dt: float):
        super().Integrate(dt)

        # update the orientation using the angular velocity
        self.UpdateRadian(dt)

        # all you need to know is that it worth two hours of debug
        self.polygon.delta_radian = self.polygon.radian - self.polygon.previous_radian; # print(self.polygon.delta_radian)  # CHECKED | SEEMS TO WORK PERFECTLY.    
        self.position = self.world_pivot + self.polygon.RotateMatrixMultiplication2D(self.position - self.world_pivot, self.polygon.delta_radian)# update the postion after having done rotating the vertices

        self.polygon.Rotating()
        self.polygon.Translating(self.position)

        self.polygon.previous_radian = self.polygon.radian
        
        # reset it every frame
        self.world_pivot[0] = self.position[0]
        self.world_pivot[1] = self.position[1]
    
    def LocalizePoint(self, value: np.ndarray): 
        value = self.polygon.RotateMatrixMultiplication2D(value - self.position, - self.polygon.radian)
        return value
    
    def ComputeAABB(self):
        assert False, "have not create it yet, and never have to! (do the CreateAABB instead!)"

    def get_AABBheight(self):
        return super().get_AABBheight()
    def get_AABBwidth(self):
        return super().get_AABBwidth()
    def get_AABBposition(self):
        return super().get_AABBposition()

    def CreateAABB(self):
        # assert False, f"do it right now ({type(self)})"
        vertex = self.polygon.world_vertices[0]
        min_x = vertex[0]
        max_x = vertex[0]

        min_y = vertex[1]
        max_y = vertex[1]

        for i in range(1, len(self.polygon.world_vertices)):
            vertex = self.polygon.world_vertices[i]
            x = vertex[0]
            y = vertex[1]

            if x < min_x:
                min_x = x
            elif x > max_x:
                max_x = x
            
            if y < min_y:
                min_y = y
            elif y > max_y:
                max_y = y
        self.AABB = AxisAlignedBoundingBox(min_x, min_y, max_x - min_x, max_y - min_y)
        self.AABBwidth = max_x - min_x
        self.AABBheight = max_y - min_y
        self.AABBposition = np.array([min_x, min_y], dtype = np.float64)


    def UpdateRadian(self, dt: float):
        self.polygon.radian += self.angular_velocity * dt

    def Draw(self, surface: pg.surface.Surface):
        pg.draw.polygon(surface, self.color, self.polygon.world_vertices)
    
    def DrawAABB(self, surface: pg.surface.Surface):
        self.rect.x = self.AABBposition[0]
        self.rect.y = self.AABBposition[1]
        self.rect.width = self.AABBwidth
        self.rect.height = self.AABBheight
        pg.draw.rect(surface, (0, 255, 255), self.rect, 1)
    
    def DebugXYOrientation(self):
        cos = np.cos(self.polygon.radian)
        sin = np.sin(self.polygon.radian)
        x_orientation = np.array((cos, sin), dtype = np.float64)
        y_orientation = np.array((-sin, cos), dtype = np.float64)
        inspecting_vectors.Add(x_orientation * 50, self.position)
        inspecting_vectors.Add(y_orientation * 50, self.position)
        

    
    

    
    
    

    

    
    