import numpy as np
import pygame as pg

from components.body import *
from components.settings import *
from components.force import *
from components.broad_phase import *
from components.narrow_phase import *

class IntegrationRegistory:
    def __init__(self):
        self.body_queue = []

    def Add(self, body):
        self.body_queue.append(body)
    
    def Integrate(self, dt: float):
        for body in self.body_queue:
            body.Integrate(dt)

    def ClearKinetic(self):
        for body in self.body_queue:
            body.ClearKinetic()

    def Clear(self):
        self.body_queue.clear()



class Scene:
    def __init__(self, 
                 
                gravitational_acceleration: np.array = GRAVITATIONAL_ACCELERATION,
                force_registory: ForceRegistory = None,
                broad_phase_queue: BroadPhaseQueue = None, 
                narrow_phase_queue: NarrowPhaseQueue = None,
                contact_solver: ContactSolver = None, 
                 ):
        self.gravitational_acceleration: np.array = gravitational_acceleration
        self.dt: float = None
        self.body_list: list[Body] = []

        self.FORCE_REGISTRY: ForceRegistory = force_registory
        self.BROAD_PHASE_QUEUE: BroadPhaseQueue = broad_phase_queue
        self.NARROW_PHASE : NarrowPhaseQueue = narrow_phase_queue
        self.CONTACT_SOLVER: ContactSolver = contact_solver
        
        self.FORCE_GENERATOR_POOL: list[ForceGenerator] = []
        self.HasSetUp = False

    def InsertForceRegistry(self, force_registory: ForceRegistory):
        if type(force_registory) != ForceRegistory:
            raise "[WARNING][SCENE][INSERT FORCE REGISTORY]: the insert object is not a force_registory!"
        self.FORCE_REGISTRY = force_registory
    def InsertBroadPhaseQueue(self, broad_phase_queue: BroadPhaseQueue):
        if type(broad_phase_queue) != BroadPhaseQueue:
            raise "[WARNING][SCENE][INSERT BROAD PHASE QUEUE]: The insert object is not a BroadPhaseQueue!"
        self.BROAD_PHASE_QUEUE = broad_phase_queue

    def InsertNarrowPhaseQueue(self, narrow_phase_queue: NarrowPhaseQueue):
        if type(narrow_phase_queue) != NarrowPhaseQueue:
            raise "[WARNING][SCENE][INSERT NARROW PHASE QUEUE]: The insert object is not a NarrowPhaseQueue!"
        self.NARROW_PHASE = narrow_phase_queue

    def InsertContactSolver(self, contact_solver: ContactSolver):
        if type(contact_solver) != ContactSolver:
            raise "[WARNING][SCENE][INSERT CONTACT SOLVER]: The insert object is not a ContactSolver!"
        self.CONTACT_SOLVER = contact_solver
    
    def InsertForceRegistry(self, force_registory: ForceRegistory):
        if not isinstance(force_registory, ForceRegistory):
            raise TypeError("[WARNING][SCENE][INSERT FORCE REGISTORY]: The insert object is not a ForceRegistory!")
        self.FORCE_REGISTRY = force_registory
    
    def InsertBroadPhaseQueue(self, broad_phase_queue: BroadPhaseQueue):
        if not isinstance(broad_phase_queue, BroadPhaseQueue):
            raise TypeError("[WARNING][SCENE][INSERT BROAD PHASE QUEUE]: The insert object is not a BroadPhaseQueue!")
        self.BROAD_PHASE_QUEUE = broad_phase_queue
    
    def InsertNarrowPhaseQueue(self, narrow_phase_queue: NarrowPhaseQueue):
        if not isinstance(narrow_phase_queue, NarrowPhaseQueue):
            raise TypeError("[WARNING][SCENE][INSERT NARROW PHASE QUEUE]: The insert object is not a NarrowPhaseQueue!")
        self.NARROW_PHASE = narrow_phase_queue
    
    def InsertContactSolver(self, contact_solver: ContactSolver):
        if not isinstance(contact_solver, ContactSolver):
            raise TypeError("[WARNING][SCENE][INSERT CONTACT SOLVER]: The insert object is not a ContactSolver!")
        self.CONTACT_SOLVER = contact_solver
    
    
    def AddForceGenerator(self, force_generator: ForceGenerator):
        if not isinstance(force_generator, ForceGenerator):
            raise "[WARNING][SCENE][INSERT FORCE GENERATOR]: The inserted object is not a ForceGenerator"
        self.FORCE_GENERATOR_POOL.append(force_generator)
        # self.FORCE_REGISTORY
        
    def ClearForceGeneratorPool(self):
        self.FORCE_GENERATOR_POOL.clear()

    def Run(self, surface: pg.surface.Surface, dt: float):
        print("This should be exclusively call as the standard run | otherwise, the functions need to be spreaded? out, because that would be more easy to debug")
        if self.HasSetUp == False:
            raise ("[Warning][Scene][Creat_A_Simulation]: Should have setup the simulation before running")
        self.FORCE_REGISTRY.Execute(dt)
        self.BroadPhase()
        self.NarrowPhase(dt)
        self.Integrate(dt)
        self.ClearKinetic()
        self.Draw(surface)
        
    def SetupBugCheck(self):
        if len(self.body_list) == 0:
            print("[WARNING][SCENE][SET UP]: Bro, no rigid body is added in the simulation ...")
        if len(self.FORCE_GENERATOR_POOL) == 0:
            print("[WARNING][SCENE][SET UP]: There is no force_generator in the simulation, please watch out!")
        if self.FORCE_REGISTRY == None:
            print("[WARNING][SCENE][SET UP]: Force registory is missed")
        if self.BROAD_PHASE_QUEUE == None:
            print("[WARNING][SCENE][SET UP]: Broad phase queue is not found")
        if self.NARROW_PHASE == None:
            print("[WARNING][SCENE][SET UP]: Narrow phase queue is not found")

    def Set_dt(self, value: float):
        self.dt = value
    def AddBody(self, body: Body):
        self.body_list.append(body)
    def RemovePreviousBody(self):
        self.body_list.pop()
    def RemoveSpecificBody(self):
        self.body_list

    def SetupForceRegistry(self):
        for body in self.body_list:
            body: Body = body
            for force_generator in self.FORCE_GENERATOR_POOL:
                self.FORCE_REGISTRY.Add(body, force_generator)
    def ForceRegistryExecute(self, dt: float):
        self.FORCE_REGISTRY.Execute(dt)

    def BroadPhase(self):
        self.BROAD_PHASE_QUEUE.GeneratingPairs(self.body_list)
    def NarrowPhase(self, dt: float):
        self.NARROW_PHASE.ResolveContactPair(self.BROAD_PHASE_QUEUE.pair_queue, self.CONTACT_SOLVER, dt)

    def Integrate(self, dt: float):
        for body in self.body_list:
            body.Integrate(dt)
    def ClearKinetic(self):
        for body in self.body_list:
            body.ClearKinetic()
    
    

    def Draw(self, surface: pg.surface.Surface):
        for body in self.body_list:
            body.Draw(surface)

        
        
    