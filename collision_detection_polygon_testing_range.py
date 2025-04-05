import numpy as np
from pygame.locals import *
import pygame as pg
import time

from components.body import *
from components.debug_kit import *
from components.scene import *

from components.settings import *

WIDTH, HEIGHT = 600, 600
FPS = 10000
clock = pg.time.Clock()

screen = pg.display.set_mode((WIDTH, HEIGHT))


SIMULATION = Scene(GRAVITATIONAL_ACCELERATION)
SIMULATION.AddForceGenerator(GRAVITY_FORCE_GENERATOR)
SIMULATION.AddForceGenerator(DRAG_FORCE_GENERATOR)

FORCE_REGISTORY = ForceRegistory()
BROADPHASEQUEUE = BroadPhaseQueue()
NARROWPHASEQUEUE = NarrowPhaseQueue()
CONTACT_SOLVER = ContactSolver()

SIMULATION = Scene(GRAVITATIONAL_ACCELERATION)
SIMULATION.AddForceGenerator(GRAVITY_FORCE_GENERATOR)
SIMULATION.AddForceGenerator(DRAG_FORCE_GENERATOR)

SIMULATION.InsertForceRegistry(FORCE_REGISTORY)
SIMULATION.InsertNarrowPhaseQueue(NARROWPHASEQUEUE)
SIMULATION.InsertBroadPhaseQueue(BROADPHASEQUEUE)
SIMULATION.InsertContactSolver(CONTACT_SOLVER)

polygon_attributes = PolygonsAttributes(vertices = [
        np.array([-128.3333282470703,-42.83332824707031], dtype = np.float64), 
        np.array([9.666671752929688,-114.83332824707031], dtype = np.float64), 
        np.array([89.66667175292969,-35.83332824707031], dtype = np.float64), 
        np.array([101.66667175292969,39.16667175292969], dtype = np.float64), 
        np.array([17.666671752929688,91.16667175292969], dtype = np.float64), 
        np.array([-90.33332824707031,63.16667175292969], dtype = np.float64), 
])

polygon_attributes1 = PolygonsAttributes(vertices = [
        np.array([-128.3333282470703,-42.83332824707031], dtype = np.float64), 
        np.array([9.666671752929688,-114.83332824707031], dtype = np.float64), 
        np.array([89.66667175292969,-35.83332824707031], dtype = np.float64),
        np.array([101.66667175292969,39.16667175292969], dtype = np.float64), 
        np.array([17.666671752929688,91.16667175292969], dtype = np.float64), 
        np.array([-90.33332824707031,63.16667175292969], dtype = np.float64), 
])

polygon_attributes2 = PolygonsAttributes(vertices = [
        np.array([-140.75,-19.5], dtype = np.float64), 
        np.array([140.75,-19.5], dtype = np.float64), 
        np.array([140.75,16.5], dtype = np.float64), 
        np.array([-140.75,16.5], dtype = np.float64), 
]
)


polygon_attributes3 = PolygonsAttributes(vertices = [
        np.array([-340.75,-19.5], dtype = np.float64), 
        np.array([340.75,-19.5], dtype = np.float64), 
        np.array([340.75,16.5], dtype = np.float64), 
        np.array([-340.75,16.5], dtype = np.float64), 
]
)
# polygon_attributes3.radian = - np.pi / 18 


polygon = BodyPolygon(
    np.array([200, 300], dtype = np.float64),
    polygon_attributes2,
    mass_data = MassData(0, 1/1_00_000),
    color = (0, 255, 0),
    is_fixed_world_pivot = True,
)
SIMULATION.AddBody(polygon)

polygon1 = BodyPolygon(
    np.array([200, 0], dtype = np.float64),
    polygon_attributes1,
    mass_data = MassData(1/20, 1/100_000),
)
SIMULATION.AddBody(polygon1)


polygon2 = BodyPolygon(
    np.array([400, 500], dtype = np.float64),
    polygon_attributes3,
    mass_data = MassData(0, 0),
)
SIMULATION.AddBody(polygon2)



SIMULATION.SetupForceRegistry()

default_velocity = 100
velocity = default_velocity

# polygon.angular_velocity = 0
# polygon.polygon.SetRadian(np.pi / 2)

running = True
previous_time = time.time()
while running:

    dt = time.time() - previous_time
    previous_time = time.time()

    # dt *= 0.55

    screen.fill((0, 0, 0))

    mouse = pg.mouse.get_pressed()
    
    for event in pg.event.get():
        if event.type == QUIT:
            running = False
        if event.type == KEYDOWN:
            if event.key == K_ESCAPE:
                running = False

    keys = pg.key.get_pressed()
    
    if keys[K_q]:
        polygon.polygon.radian -= 2 * dt
    if keys[K_e]:
        polygon.polygon.radian += 2 * dt

    if keys[K_LSHIFT]:
        velocity *= 2.7

    if keys[K_w]:
        polygon.position[1] -= velocity * dt
    if keys[K_s]:
        polygon.position[1] += velocity * dt
    if keys[K_a]:
        polygon.position[0] -= velocity * dt
    if keys[K_d]:
        polygon.position[0] += velocity * dt
    
    velocity = default_velocity

    if len(SIMULATION.FORCE_REGISTRY.queue) != len(SIMULATION.FORCE_GENERATOR_POOL) * len(SIMULATION.body_list):
        SIMULATION.FORCE_REGISTRY.Clear()
        SIMULATION.SetupForceRegistry()
    SIMULATION.ForceRegistryExecute(dt)

    for _ in range(3):
        SIMULATION.BroadPhase()
        SIMULATION.NarrowPhase(dt)

    # inspecting_points.Add(polygon.world_pivot, "purple")
    # inspecting_points.Add(polygon1.world_pivot, "purple")
    
    SIMULATION.Integrate(dt)
    SIMULATION.ClearKinetic()
    SIMULATION.Draw(screen)

    polygon.CreateAABB()
    


    # polygon.DrawAABB(screen)
    # polygon.polygon.DebugNormalVector()
    
    # polygon1.DrawAABB(screen)
    # polygon1.polygon.DebugNormalVector()

    inspecting_texts.Draw(screen)
    inspecting_points.Draw(screen)
    inspecting_vectors.Draw(screen)
    inspecting_lines.Draw(screen)


    inspecting_texts.Clear()
    inspecting_points.Clear()
    inspecting_vectors.Clear()
    inspecting_lines.Clear()
    warning_logs.Reset()

    BROADPHASEQUEUE.Clear()

    clock.tick(FPS)
    pg.display.set_caption(f"{clock.get_fps() // 1}")
    pg.display.flip()

print("YOU RE CURRENTLY AT THE LAST VERSION I THINK")
"""
UPDATE LOG:
 30-3-2025:
  - THE POLYGON COLLISION LOOOKS CLEAN. HOWEVER ...
  - THE FRICTION ONE IS BROKEN (STIL WORKS WITH FIX-ORIENTED SHAPES).
  - THE IMPULSE CALCULATION IS SILL NOT GOOD (IF THE OBJECT MASS IS LARGE, THE CONTACT SOLVER MAKES IT VIBRATING LIKE CRAZY).
  QUESTIONS REMAIN:
  - SINCE IT HAVE TO ABILITY TO ROTATE, WHAT IF THEY GOT INTO SITUATION WHERE THEY BEING STUCK BETWEEN FIXED OBJECTS (IT WILL PROBABLY EXPLODE VIOLENTLY).
  - THE IMPULSE SYSTEM IS WRONG, (AS YOU CAN SEE) | YET, I DO NOT KNOW WHY IT JUST HAPPENS. 
AS THE MATTER OF FACTS, THIS PROJECT STILL SATISFIED THE OBJECTIVE DRAWN OUT FROM THE START - LEARN HOW OBJECT INTEGRATES IN GAMES.
I THINK GO DEEPER DOWN WOULD JEOPARDISE MY ACADEMEIC PERFORMANCE FROM SCHOOL (GOTTA LOCK IN THIS SEMESTER). THUS, I WOULD ONLY STRIVE FORM SIMPLE MATH PROJECTS OR READING BOOKS.

IN CONCLUSION, I WILL JUST LOOKING FOR PROOFS, PUTTING ALL MY THOUGHT INTO DOCUMENTS THIS TIME | LOOKING FORWARD TO DO IT AGAIN AS SOON AS I HAS DONE THE COLLEGE ENTRANCE TEST. 
SOL:
 0- ISOLATION

"""
