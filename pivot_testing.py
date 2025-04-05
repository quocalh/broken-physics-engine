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

polygon = BodyPolygon(
    np.array([200, 300], dtype = np.float64),
    polygon_attributes2,
    # mass_data = MassData(0, 0),
    mass_data = MassData(0, 1/100_000_000),
    # mass_data = MassData(1/1_000_000, 0),
    color = (0, 255, 0) 
)
SIMULATION.AddBody(polygon)

polygon1 = BodyPolygon(
    np.array([200, 0], dtype = np.float64),
    polygon_attributes1,
    # mass_data = MassData(1/200, 1/100_000),
    mass_data = MassData(1/200, 0),
    # mass_data = MassData(0, 1/200),
)
SIMULATION.AddBody(polygon1)

# polygon.mass_data.inverse_inertia = 0
# polygon1.mass_data.inverse_inertia = 0

SIMULATION.SetupForceRegistry()

default_velocity = 100
velocity = default_velocity

# ============================ TESTING VALUE ===============================
world_pivot = np.zeros(2, dtype = np.float64)
world_pivot[0] = polygon.position[0]; world_pivot[1] = polygon.position[1]

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
            if event.key == K_f:
                polygon.polygon.radian = 0
        

    keys = pg.key.get_pressed()
    
    if mouse[0]:
        world_pivot[0] = pg.mouse.get_pos()[0]
        world_pivot[1] = pg.mouse.get_pos()[1]
        world_pivot[0] = polygon.position[0] ; world_pivot[1] = polygon.position[1]
        

    if keys[K_q]:
        polygon.polygon.radian -= 5 * dt
    if keys[K_e]:
        polygon.polygon.radian += 5 * dt
    
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

    # usually, we wont need to manually turn it on, it will reveal itself in the contact_solver (contact solver deal with this thing)
    
    polygon.world_pivot = world_pivot
    SIMULATION.Integrate(dt)
    SIMULATION.ClearKinetic()
    SIMULATION.Draw(screen)

    inspecting_points.Add(polygon.position)
    inspecting_points.Add(polygon.world_pivot, "blue")
    inspecting_points.Add(polygon.LocalizePoint(polygon.world_pivot), "green")

    polygon.CreateAABB()
    
    polygon.DrawAABB(screen)
    # polygon.polygon.DebugNormalVector()
    # polygon.DebugXYOrientation()
    

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

print()
print("NOTE: THIS IS SUPPOSED TO BE THE FINAL CHECKPOINT, ONLY STRIVE TO SOLVE FOR ROTATION BUGS")
print()
# print("THIS ATTEMPT REPRESENTS THE ANGULAR VELOCITY IN CHARGE, BUT STILL YIELD UNREALISTIC RESULT IF THE PENERTRATION OBJECT' INERTIA IS 0 (THESIS: SUFFER FROM WRONG ROTATION PIVOT)")


"""
    UPDATE LOG 20/3/25:
    BUG DETECTED: 
    - LOCAL VERTICES ARE NOT SUPPOSED TO MOVE
    - AFTER DEBUGGING SESSION, SHIT GOT WORSE, NOW IT EVEN DOES NOT ROTATE
    - RADIAN: 2 0 2 0 2 0 2 0 # I DONT KNOW HOW THE FUCK IT RETURNS A 0 EVERY FUCKING FRAME :SOB:
    shit is incorrect after added tangents| before rotating function  (likely to be wrong though)
    UPDATE LOG 21/3/2025:
     - THE EXAM IS DUE TOMMOROW :SOB:, 
     - WHATEVER, THE BUG IS FORTUNATELY FIXED (DUDE WAS REALLY THAT STUPID, HE FORGOT AN ADDITION IN THE TRANSLATION).
        10:33 AM
     - SHIT IS RUNNING GOOD -> PROCEED TO NEXT STEP
        9:37PM
     - JUST CONSIDER MOVE ON TO THE NEXT CHECKPOINT
     - ALMOST DONE! THE PROBLEM: THE PENETRATION VALUE SOMEHOW KIND OF REVERSED, THE FORMUA MUST HAVE BEEN WRONG SOMEWHERE
     - AFTER DONE WITH POLYGON COLLISION, I WOULD UPLOAD THE FILE TO GITHUB :d

    UPDATE LOG 26/3/2025:
     - THE ROTATION IS NOT LOOKING GOOD (FORMULA)
     - FIX THE COUPLE OF BUGS | STILL SUFFER FROM EXCESSIVE ROTATION
     - SOME SCENARIOS DOES HIDE THE BUG (ROTATION IMPULSE)
        1:16 NEXT DAY
     - ALMOST LOOKING GOOD, THE ROTATION MAKE THE OBJECT GAINS ENERY MID WAY, TOO BAD
     - EXCESS AMOUNT OF ROT EXERTED ON THE OBJECT
     - ALMOST PERFECT
    UPDATE LOG 28/3/2025:
     - ERM SOMETIMES IN COLLSION RESPONES, WE DONT ROTATE THE DTHING IN COM
     - IF IT HAVE A PEIEC OF PERICE = > ROT AXIS IS A CONTACT POINT
      - IF IT IS A PASIVELY HIT BYT A CORNER - ROT IN CENTER
"""
