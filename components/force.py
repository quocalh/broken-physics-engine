from components.settings import GRAVITATIONAL_ACCELERATION
from components.body import *


class ForceGenerator:
    def __init__(self):
        pass
    def ApplyUpdate(self, body: Body):
        print("Applying Force ... ")
    
class GravityForceGenerator(ForceGenerator):   
    def __init__(self, gravitational_acceleration: np.array = None):
        super().__init__()
        if gravitational_acceleration == None:
            self.gravitational_acceleration = GRAVITATIONAL_ACCELERATION
        else:
            self.gravitational_acceleration = gravitational_acceleration
    
    def SetGravitationalAcceleration(self, value = np.array):
        self.gravitational_acceleration = value
        
    def ApplyUpdate(self, body: Body):
        if body.mass_data.iMass == 0:
            return
        GravityForce = self.gravitational_acceleration * body.mass_data.mass
        body.AccumulatingForce(GravityForce)

class DragForceGenerator(ForceGenerator):
    def __init__(self, k1: float = 0.01, k2: float = 0.03):
        self.k1: float = k1
        self.k2: float = k2
    
    def SetCoefficients(self, k1_value, k2_value):
        self.k1: float = k1_value
        self.k2: float = k2_value

    def ApplyUpdate(self, body: Body):
        Speed = np.linalg.norm(body.velocity)
        if Speed == 0:
            return
        DragForceMagnitude = self.k1 * Speed + self.k2 * Speed * Speed
        DragForce = -(body.velocity / Speed) * DragForceMagnitude
        body.AccumulatingForce(DragForce)
    """
    FRICITON IS VELOCITY-DEPENDENT | LOL
    """

class ForceRegistoryPair:
    def __init__(self, body: Body, force_generator: ForceGenerator):
        self.body: Body = body
        self.force_generator = force_generator

class ForceRegistory:
    def __init__(self):
        self.queue = []
    
    def Add(self, body: Body, force_generator: ForceGenerator):
        self.queue.append(ForceRegistoryPair(body, force_generator))
    
    def Execute(self, dt: float):
        # have not use dt parameter yet 
        for pair in self.queue:
            pair: ForceRegistoryPair = pair
            pair.force_generator.ApplyUpdate(pair.body)

    def Clear(self):
        self.queue = []
        
GRAVITY_FORCE_GENERATOR = GravityForceGenerator()
DRAG_FORCE_GENERATOR = DragForceGenerator()