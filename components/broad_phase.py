import numpy as np

from components.narrow_phase import *
from components.body import *
from components.debug_kit import *
"""

██████╗░██████╗░░█████╗░░█████╗░██████╗░  ██████╗░██╗░░██╗░█████╗░░██████╗███████╗
██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔══██╗  ██╔══██╗██║░░██║██╔══██╗██╔════╝██╔════╝
██████╦╝██████╔╝██║░░██║███████║██║░░██║  ██████╔╝███████║███████║╚█████╗░█████╗░░
██╔══██╗██╔══██╗██║░░██║██╔══██║██║░░██║  ██╔═══╝░██╔══██║██╔══██║░╚═══██╗██╔══╝░░
██████╦╝██║░░██║╚█████╔╝██║░░██║██████╔╝  ██║░░░░░██║░░██║██║░░██║██████╔╝███████╗
╚═════╝░╚═╝░░╚═╝░╚════╝░╚═╝░░╚═╝╚═════╝░  ╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═════╝░╚══════╝
"""

pg.init()

class Pair:
    def __init__(self, body1, body2):
        self.object1 = body1
        self.object2 = body2

class BroadPhaseQueue:
    def __init__(self):
        self.pair_queue: list = []

    def Clear(self):
        self.pair_queue = []

    def GeneratingPairs(self, array_of_bodies: list[Body]):
        for i in range(len(array_of_bodies)):
            for j in range(i + 1, len(array_of_bodies)):
                body1 = array_of_bodies[i]
                body2 = array_of_bodies[j]
                
                # layer check
                "if a.layer != b.layer => skip"

                body1.CreateAABB()
                body2.CreateAABB()
                
                if is_collided_AABB_vs_AABB_AABB(body1.AABB, body2.AABB):
                    self.pair_queue.append(Pair(body1, body2))
                

class NarrowPhaseQueue:
    CollisionTypeCallBack_dict = {
        BodyCircle: {
            BodyCircle: COLLSION_CIRCLE_VS_CIRCLE,
            BodyRect: COLLISION_CIRCLE_VS_AABB, 
            # polygons are comming soon
        },
        BodyRect: {
            BodyCircle: COLLISION_AABB_VS_CIRCLE,
            BodyRect: COLLSION_AABB_VS_AABB, 
            # polygons are comming soon
        },
        BodyPolygon: {
            BodyPolygon: COLLISION_POLYGON_VS_POLYGON,
        }
    }
    def __init__(self):
        pass

    def ResolveContactPair(self, broad_phase_pairs: list[Pair], contact_solver: ContactSolver, dt: float):
        for pair in broad_phase_pairs:
            body1 = pair.object1
            body2 = pair.object2
            collision: Collision = self.CollisionTypeCallBack_dict[type(body1)][type(body2)]
            
            collision.FillingObjects(body1, body2)
            
            collision.SetIndependenceRun(contact_solver, dt)
        return

            