import numpy as np
import pygame as pg
from pygame.locals import *
import time

from components.settings import *
from components.body import *
from components.force import *
from components.narrow_phase import *

from components.debug_kit import *
from components.broad_phase import *
from components.scene import *

running = True

FPS = 100000

size = WIDTH, HEIGHT
screen = pg.display.set_mode(size)
clock = pg.time.Clock()
# print(np.random.uniform(1, 10))

DRAG_FORCE_GENERATOR.SetCoefficients(0.01, 0.03)
GRAVITY_FORCE_GENERATOR.SetGravitationalAcceleration(GRAVITATIONAL_ACCELERATION)

FORCE_REGISTORY = ForceRegistory()
BROADPHASEQUEUE = BroadPhaseQueue()
NARROWPHASEQUEUE = NarrowPhaseQueue()
CONTACT_SOLVER = ContactSolver()


SIMULATION = Scene(GRAVITATIONAL_ACCELERATION)
SIMULATION.AddForceGenerator(GRAVITY_FORCE_GENERATOR)
SIMULATION.AddForceGenerator(DRAG_FORCE_GENERATOR)

SIMULATION.InsertForceRegistory(FORCE_REGISTORY)
SIMULATION.InsertNarrowPhaseQueue(NARROWPHASEQUEUE)
SIMULATION.InsertBroadPhaseQueue(BROADPHASEQUEUE)
SIMULATION.InsertContactSolver(CONTACT_SOLVER)

platform = BodyRect(np.array([WIDTH / 2, 350], dtype = np.float64),
                    WIDTH, 40, Materials(0.0),  MassData(0, 10),
                    color = (0, 255, 134),
                    )
SIMULATION.AddBody(platform)
dirt = BodyRect(np.array([WIDTH / 2, 390], dtype = np.float64),
                    WIDTH, 40, mass_data = MassData(0, 10),
                    color = (150, 75, 0),
                    )
SIMULATION.AddBody(dirt)
box = BodyRect(np.array([120, 255], dtype = np.float64),
                    50, 155, Materials(0.0),
                    MassData(0), 
                    color = (85,144,194),
                    )
SIMULATION.AddBody(box)
box = BodyRect(np.array([350, 255], dtype = np.float64),
                    50, 155, Materials(0.0),
                    MassData(0), 
                    color = (85,144,194),
                    )
SIMULATION.AddBody(box)
wheel = BodyCircle(np.array([180, 50], dtype = np.float64),
                   15, Materials(0.0),
                   MassData(1/20, 10), 
                   color = (165,30,34),
                   )
SIMULATION.AddBody(wheel)
wheel_1 = BodyCircle(np.array([190, 150], dtype = np.float64),
                   25, Materials(0.0),
                   MassData(1/10, 10), 
                   color = (165,30,34),
                   )
SIMULATION.AddBody(wheel_1)
wheel_2 = BodyCircle(np.array([180, 150], dtype = np.float64),
                   45, Materials(0.0),
                   MassData(1/10, 10), 
                   color = (165,30,34),
                   )
SIMULATION.AddBody(wheel_2)

# SIMULATION.body_list = [platform, box, wheel, wheel_1, wheel_2]

SIMULATION.SetupForceRegistory()
SIMULATION.SetupBugCheck()

default_slowmo_coef = 1
default_slowmo_divisor = 1/5
slowmo_coef = default_slowmo_coef # 1 to 0.000000000001 (use to slow time down)
slowmo_toggle = False

PreviousTime = time.time()
while running:
    dt = time.time() - PreviousTime
    PreviousTime = time.time()
 
    dt = dt * slowmo_coef

    screen.fill((0, 0, 0))
    
    keys = pg.key.get_pressed()
    mouse = pg.mouse.get_pressed()
    
    for event in pg.event.get():
        if event.type == QUIT:
            running = False
        if event.type == KEYDOWN:
            if event.key == K_ESCAPE:
                running = False
        if event.type == pg.MOUSEBUTTONDOWN:
            NewBody = BodyCircle(np.array(pg.mouse.get_pos(), dtype = np.float64),
                   15, Materials(0.0),
                   MassData(1/1, 10), 
                   color = (165,30,34),
                   )
            SIMULATION.AddBody(NewBody)
    
    """
    As long as the self.bodylist is not changed (being added or deleted with new ones), the force registor can stay the same
    """
    if len(SIMULATION.FORCE_REGISTORY.queue) != len(SIMULATION.FORCE_GENERATOR_POOL) * len(SIMULATION.body_list):
        SIMULATION.FORCE_REGISTORY.Clear()
        SIMULATION.SetupForceRegistory()

    slowmo_toggle = False
    if keys[K_SPACE]:
        slowmo_toggle = True
        inspecting_texts.ManualLog(f"slowed by {default_slowmo_divisor}x (be cautious with float errors!)", screen, (20, 10), center = False)
    if slowmo_toggle:
        slowmo_coef = default_slowmo_coef / default_slowmo_divisor           
    else:
        slowmo_coef = default_slowmo_coef

    
    SIMULATION.ForceRegistoryExecute(dt)
    # FORCE_REGISTORY.Add(box            , GRAVITY_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(box            , DRAG_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel          , GRAVITY_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel          , DRAG_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel_1        , GRAVITY_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel_1        , DRAG_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel_2        , DRAG_FORCE_GENERATOR)
    # FORCE_REGISTORY.Add(wheel_2        , GRAVITY_FORCE_GENERATOR )
    # FORCE_REGISTORY.Execute(dt)
    # FORCE_REGISTORY.Clear()

    for _ in range(4):
        SIMULATION.BroadPhase()
        SIMULATION.NarrowPhase(dt)
    # BROADPHASEQUEUE.GeneratingPairs([box, wheel, wheel_1, platform, wheel_2])
    
    # NARROWPHASEQUEUE.ResolveContactPair(BROADPHASEQUEUE.pair_queue, CONTACT_SOLVER, dt)

    SIMULATION.Integrate(dt)
    # box.                        Integrate(dt)
    # wheel.                      Integrate(dt)
    # platform.                   Integrate(dt)
    # wheel_1.                    Integrate(dt)
    # wheel_2.                    Integrate(dt)

    SIMULATION.ClearKinetic()
    # box.                        ClearKinetic()
    # wheel.                      ClearKinetic()
    # platform.                   ClearKinetic()
    # wheel_1.                    ClearKinetic()
    # wheel_2.                    ClearKinetic()
    
    SIMULATION.Draw(screen)
    # platform.                   Draw(screen)
    # dirt.                       Draw(screen)
    # wheel.                      Draw(screen)
    # box.                        Draw(screen)
    # wheel_1.                    Draw(screen) 
    # wheel_2.                    Draw(screen)
    
    # debug
    # inspecting_vectors.Add(wheel.velocity / 5, wheel.position)     
    # inspecting_vectors.Add(box.velocity / 5, box.position)
    # inspecting_vectors.Add(wheel_1.velocity / 5, wheel_1.position)
    # inspecting_vectors.Add(wheel_2.velocity / 5, wheel_2.position)
    
    inspecting_points.Draw(screen)
    inspecting_texts.Draw(screen)
    inspecting_vectors.Draw(screen)
    inspecting_lines.Draw(screen)

    inspecting_points.Clear()
    inspecting_texts.Clear()
    inspecting_vectors.Clear()
    inspecting_lines.Clear()

    BROADPHASEQUEUE.Clear()

    pg.display.flip()
    pg.display.set_caption(f"{clock.get_fps() // 1}")
    clock.tick(FPS)

SIMULATION.ClearForceGeneratorPool()
