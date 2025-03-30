import pygame as pg
from pygame.locals import *
import numpy as np
import time
"""
THIS FILE IS CREATED TO GENERATE POLYGONS | I HAVE NOT DETERMINED HOW IT WILL EXPORT THE VERTICES FILE THOUGH
WHAT I KNOW FOR SURE IS THAT IT CAN BE VARIED
"""

pg.init()
from components.debug_kit import *
inspecting_texts = InspectingTexts((255, 0, 0))

class Points:
    def __init__(self, position: np.array, 
                 color: tuple[int, int, int] = (164,30,34),
                 ):
        self.position = position
        self.color = color
        self.radius = 5

    def Draw(self, surface: pg.surface.Surface):
        pg.draw.circle(surface, self.color, self.position, self.radius)

class PolygonVertices:
    def __init__(self,
                 center_of_mass: np.array = np.zeros(2, dtype = np.float32)
                 ):
        assert type(center_of_mass) != np.array, f"[POLYGON_VERTIECS][INIT]: center of mass suppose to be a np.array type not {type(center_of_mass)}!"
        self.center_of_mass: np.array = center_of_mass
        self.center_of_mass_point: Points = Points(self.center_of_mass)
        self.coordinate_points: list[Points] = [] # the vertex data purelty retrived from pygame application, have not been through translation correction
        self.localized_vertices: list[np.array] = [] # the final result, the exporting file | would be given with a copy afterwards
        self.is_convex_shape = True
        self.imaginary_mass_center = np.zeros(2, dtype = np.float32)
        self.is_clockwise: bool = None
        # self.SetClockwise()
        # self.sign: int = 1 - 2 * (int(not self.is_clockwise)) 
        # print(self.sign)

    def SetClockwise(self, is_shape_check: bool = False):
        if is_shape_check:
            if self.IsShape(): return
        len_coordinate_points = len(self.coordinate_points)
        tail_vertex  = self.coordinate_points[len_coordinate_points - 1].position
        pretail_vertex = self.coordinate_points[len_coordinate_points - 2].position
        prepretail_vertex = self.coordinate_points[len_coordinate_points - 3].position
        pre_vector = tail_vertex - pretail_vertex
        prepre_vector = prepretail_vertex - pretail_vertex
        cross_product = self.CrossProduct(pre_vector, prepre_vector)
        self.is_clockwise  = cross_product > 0

    def AddPoint(self, point: Points):
        assert type(point) == Points, f"[POLYGON_VERTIECS][ADDING POINTS]: center of mass suppose to be a np.array type not {type(point)}!"
        self.imaginary_mass_center = self.imaginary_mass_center * len(self.coordinate_points) + point.position
        self.imaginary_mass_center /= len(self.coordinate_points) + 1
        self.coordinate_points.append(point)

    def IsShape(self):
        return len(self.coordinate_points) > 2
    
    @staticmethod
    def CrossProduct(a: np.ndarray, b: np.ndarray) -> float:
        z_component = a[0] * b[1] - a[1] * b[0]
        return z_component
    
    def ThresholdConvexCheck(self):
        if not self.IsShape():
            return
        len_coordinate_vertex = len(self.coordinate_points)

        previous_position = self.coordinate_points[len_coordinate_vertex - 2].position
        current_position = self.coordinate_points[len_coordinate_vertex - 1].position
        new_position = self.coordinate_points[0].position

        self.SetClockwise()
        print(self.is_clockwise)

        a = new_position - current_position
        b = previous_position - current_position
        z_component = self.CrossProduct(a, b)
        if z_component <= 0: # if angle between a and b is 
            self.is_convex_shape = False    

    def ThresholdConvexCheck(self):
        """
        WOULD ABSTRACT IT SOON, NOT NOW
        """
        if not self.IsShape():
            return
        len_coordinate_vertex = len(self.coordinate_points)
        previous_position = self.coordinate_points[len_coordinate_vertex - 2].position
        current_position = self.coordinate_points[len_coordinate_vertex - 1].position
        new_position = self.coordinate_points[0].position
        next_position = self.coordinate_points[1].position

        self.SetClockwise()
        clock_wise_coef = 1 - 2 * int(not self.is_clockwise)

        a = new_position - current_position
        b = previous_position - current_position
        z_component = self.CrossProduct(a, b * clock_wise_coef)
        if z_component <= 0: 
            self.is_convex_shape = False    

        a = next_position - new_position
        b = current_position - new_position
        z_component = self.CrossProduct(a, b * clock_wise_coef)
        if z_component <= 0: 
            self.is_convex_shape = False    
         

            
    def SetCenterOfMass(self, position: np.array):
        assert type(position) == np.ndarray, "Nigga, that is not a numpy array!"
        self.center_of_mass = position
        self.center_of_mass_point.position = position

    def Draw(self, surface: pg.surface.Surface):
        if not self.is_convex_shape:
            inspecting_texts.Add("i think that is not a convex shape (maybe yes, this checking algorithm suck) :sob:", (10, 10))
            inspecting_texts.Add("our engine succumbed to those concave creatures", (10, 40))

        pg.draw.circle(surface, self.center_of_mass_point.color, self.center_of_mass_point.position.tolist(), self.center_of_mass_point.radius)
        len_coordinate_vertex = len(self.coordinate_points)

        if not self.IsShape():
            for vertex in self.coordinate_points:
                pg.draw.circle(screen, vertex.color, vertex.position.tolist(), vertex.radius + 2)
            return

        for i in range(len_coordinate_vertex):
            next_i = 0 if (i + 1 == len_coordinate_vertex) else i + 1 # sealing up the vertex coordinate
            pg.draw.line(screen, (85,134,164), self.coordinate_points[i].position.tolist(), self.coordinate_points[next_i].position.tolist(), 2)
            vertex: Points = self.coordinate_points[i]
            pg.draw.circle(screen, vertex.color, vertex.position.tolist(), vertex.radius)

    def DrawImaginaryMassCenter(self, surface: pg.Surface):
        pg.draw.circle(surface, (0, 255, 0), self.imaginary_mass_center.tolist(), 2)

    def TranslatingCenterOfMass(self, position: np.array):
        assert type(position) == np.ndarray, "Nigga, that is not a numpy array!"
        coordinated_vertices = []
        for vertex in self.localized_vertices:
            vertex: np.array = vertex
            vertex[0] += position
            vertex[1] += position
            coordinated_vertices.append(vertex)
        return coordinated_vertices
        
    def LocalizingCoordinate_pylist(self):
        """
        We translate the general coordinate system to the object's local coordinate system.
        py_list flag makes sure that the list only contains pythonic lists.
        """
        for vertex in self.coordinate_points:
            vertex : Points = vertex
            localized_vertex = - self.center_of_mass + vertex.position
            self.localized_vertices.append(localized_vertex.tolist())
    def LocalizingCoordinate_numpy(self):
        """
        We translate the general coordinate system to the object's local coordinate system.
        numpy flag tells that the element of the list is numpy list (type -> np.ndarray).
        """
        for vertex in self.coordinate_points:
            vertex : Points = vertex
            localized_vertex = - self.center_of_mass + vertex.position
            self.localized_vertices.append(localized_vertex)

    def ExportAsFile(self, file_name: str = "polygon_vertices"):
        pass
    
    def ExportAsLog(self, var_name: str = "vertices_array", str_type = "np.float64"):
        """
        vertices = [array([-154., -108.]), array([ 91., -90.]), array([ 52., 150.])]
        """
        self.LocalizingCoordinate_numpy()
        sum_string = f"{var_name} = "
        sum_string += "[\n"
        for index, vertex in enumerate(self.localized_vertices):
            sum_string += "\t"
            sum_string += "np.array("
            sum_string += f"[{vertex[0]},{vertex[1]}], dtype = {str_type}"
            sum_string += ")"
            sum_string += ", " 
            if index != len(self.localized_vertices) - 1:
                sum_string += "\n"
        sum_string += "\n]"
        print(sum_string)
        return

    def Export(self):
        self.LocalizingCoordinate()
        print(self.localized_vertices)

    def Import(self, vertices_list: list[np.array]):
        self.localized_vertices = vertices_list        

WIDTH, HEIGHT = 400, 400
FPS = 200
clock = pg.time.Clock()

screen = pg.display.set_mode((WIDTH, HEIGHT))

polygon_vertices = PolygonVertices(np.array([WIDTH / 2, HEIGHT / 2], dtype = np.float64))

running = True
while running:
    screen.fill((0, 0, 0))
    for event in pg.event.get():
        if event.type == QUIT:
            running = False
        if event.type == KEYDOWN:
            if event.key == K_ESCAPE:
                running = False
            if event.key == K_q:
                polygon_vertices.center_of_mass[0] = polygon_vertices.imaginary_mass_center[0]
                polygon_vertices.center_of_mass[1] = polygon_vertices.imaginary_mass_center[1]
            if event.key == K_m:
                polygon_vertices.ExportAsLog("vertices")
        if event.type == MOUSEBUTTONDOWN:
            if event.button == BUTTON_LEFT:
                vertex = Points(np.array(pg.mouse.get_pos(), dtype = np.float32))
                polygon_vertices.AddPoint(vertex)
                polygon_vertices.ThresholdConvexCheck()
    
    mouse = pg.mouse.get_pressed()
    if mouse[2]:
        polygon_vertices.SetCenterOfMass(np.array(pg.mouse.get_pos(), dtype = np.float32))
    
    polygon_vertices.Draw(screen)
    polygon_vertices.DrawImaginaryMassCenter(screen)
    inspecting_texts.Add("Right click to offset the center of mass", (100, HEIGHT - 30))
    inspecting_texts.Add("Press Q to make the mass center matchs the green dot!", (100, HEIGHT - 60))
    inspecting_texts.Add("Press M to export the local vertices", (100, HEIGHT - 90))

    inspecting_texts.Draw(screen)
    inspecting_points.Draw(screen)

    inspecting_texts.Clear()
    inspecting_points.Clear()

    clock.tick(FPS)
    pg.display.set_caption(f"{clock.get_fps() // 1}")
    pg.display.flip()

