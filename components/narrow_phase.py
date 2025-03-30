import numpy as np

from components.body import *
from components.body import Body
from components.settings import *
from components.debug_kit import *
pg.init()

"""

░█████╗░░█████╗░███╗░░██╗████████╗░█████╗░░█████╗░████████╗  ░██████╗░█████╗░██╗░░░░░██╗░░░██╗███████╗██████╗░
██╔══██╗██╔══██╗████╗░██║╚══██╔══╝██╔══██╗██╔══██╗╚══██╔══╝  ██╔════╝██╔══██╗██║░░░░░██║░░░██║██╔════╝██╔══██╗
██║░░╚═╝██║░░██║██╔██╗██║░░░██║░░░███████║██║░░╚═╝░░░██║░░░  ╚█████╗░██║░░██║██║░░░░░╚██╗░██╔╝█████╗░░██████╔╝
██║░░██╗██║░░██║██║╚████║░░░██║░░░██╔══██║██║░░██╗░░░██║░░░  ░╚═══██╗██║░░██║██║░░░░░░╚████╔╝░██╔══╝░░██╔══██╗
╚█████╔╝╚█████╔╝██║░╚███║░░░██║░░░██║░░██║╚█████╔╝░░░██║░░░  ██████╔╝╚█████╔╝███████╗░░╚██╔╝░░███████╗██║░░██║
░╚════╝░░╚════╝░╚═╝░░╚══╝░░░╚═╝░░░╚═╝░░╚═╝░╚════╝░░░░╚═╝░░░  ╚═════╝░░╚════╝░╚══════╝░░░╚═╝░░░╚══════╝╚═╝░░╚═╝
"""

class ContactSolver:
    """
    FREELANCER CLASSES  
    """
    def __init__(self):
        self.penetration: float = None
        self.contact_normal: np.array  = None
        self.object1: Body = None
        self.object2: Body = None
        self.restitution = 0
        self.impulse: float = None

        self.contact_point: np.ndarray = None # exclusive 2D abitrary for polygon rotation (torque)
        self.piercing_index: int = None
        
    
    """
    THE CLASS'S FUNCTION ONLY VIABLE WHEN THE EGNINE HAS DONE FILLING THE GIVEN PARAMETER
    """

    def CurrentLength(self):
        vector: np.array = self.object1.position - self.object2.position
        return np.linalg.norm(vector)
    
    def ResolveContact(self, dt: float):
        self.ResolveVelocity(dt)
        self.ResolvePenetration(dt)
        # self.ResolveFriction(dt) # broken, only works with fixed orientated shapes

    def ResolveVelocity(self, dt: float):        
        iMass1 = self.object1.mass_data.iMass
        iMass2 = self.object2.mass_data.iMass
        iInertia1 = self.object1.mass_data.inverse_inertia
        iInertia2 = self.object2.mass_data.inverse_inertia

        if iMass1 + iMass2 == 0:
            return
        
        # polygon oriented solver
        if type(self.object1) == BodyPolygon or type(self.object2) == BodyPolygon:
            radius1 = self.contact_point - self.object1.position
            radius2 = self.contact_point - self.object2.position
        
        vector = self.object1.velocity - self.object2.velocity # reversed sepvel and reversed contact_normal

        if type(self.object1) == BodyPolygon or type(self.object2) == BodyPolygon:
            vector = self.object1.velocity - CrossProductScalar(radius1, self.object1.angular_velocity) - (
                     self.object2.velocity - CrossProductScalar(radius2, self.object2.angular_velocity))
        
        SeparatingVelocity = np.dot(self.contact_normal, vector)
        
        if SeparatingVelocity > 0:
            return
        
        self.restitution = max(self.object1.material.restitution, self.object2.material.restitution)
        Impulse = -(1 + self.restitution) *  SeparatingVelocity

        try:
            inspecting_lines.Add(self.object1.position, radius1, True)
            inspecting_lines.Add(self.object2.position, radius2, True)
        except: pass

        AngularImpulse = 0
        AngularImpulse1 = 0
        AngularImpulse2 = 0

        if type(self.object1) == BodyPolygon:
            raCrossN = CrossProductVector(radius1, self.contact_normal)
            AngularImpulse1 = raCrossN * raCrossN * iInertia1            
        if type(self.object2) == BodyPolygon:
            rbCrossN = CrossProductVector(radius2, self.contact_normal)
            AngularImpulse2 = rbCrossN * rbCrossN * iInertia2

        AngularImpulse = AngularImpulse1 + AngularImpulse2
        
        Impulse /=  (iMass1 + iMass2) + AngularImpulse

        self.impulse = Impulse
        Impulse *= self.contact_normal
        
        self.object1.velocity += iMass1 * Impulse
        self.object2.velocity -= iMass2 * Impulse

        if type(self.object1) == BodyPolygon:
            if self.piercing_index == 0 and not self.object1.is_fixed_world_pivot:
                self.object1.world_pivot[0] = self.contact_point[0]
                self.object1.world_pivot[1] = self.contact_point[1]
            self.object1.angular_velocity += iInertia1 * CrossProductVector(radius1, Impulse)
        
        if type(self.object2) == BodyPolygon:
            if self.piercing_index == 1 and not self.object2.is_fixed_world_pivot:
                self.object2.world_pivot[0] = self.contact_point[0]
                self.object2.world_pivot[1] = self.contact_point[1]
            self.object2.angular_velocity -= iInertia2 * CrossProductVector(radius2, Impulse)

    def ResolvePenetration(self, dt: float):
        iMass1 = self.object1.mass_data.iMass
        iMass2 = self.object2.mass_data.iMass
        """
        Ian millington book represents this really articulate. 
        The formula is just has when through some arithmetic process, the proof is in da notebook
        """
        if self.penetration <= 0:
            return
        totalImass = iMass1 + iMass2
        if totalImass == 0:
            return
        MovementPerImass = self.contact_normal * self.penetration / totalImass
        self.object1.AccumulatingPosition(iMass1 * MovementPerImass)
        self.object2.AccumulatingPosition(- iMass2 * MovementPerImass)
        
    def ResolveFriction(self, dt):  
        iMass1 = self.object1.mass_data.iMass
        iMass2 = self.object2.mass_data.iMass
        # friction_coef = self.object1.friction_static / 2 + self.object2.friction_static / 2 # suppose to get pythagoreaned | but calc a average instead
        self.restitution = max(self.object1.material.restitution, self.object2.material.restitution)
        
        RelativeVelocity = self.object2.velocity - self.object1.velocity 

        Tangent = RelativeVelocity - (np.dot(self.contact_normal, RelativeVelocity)) * self.contact_normal
        TangentMagnitude = np.linalg.norm(Tangent)
        if TangentMagnitude == 0:
            return
        Tangent /= TangentMagnitude

        SeperateVelocity = np.dot(self.contact_normal, RelativeVelocity)
        # SeperateVelocity = np.dot(RelativeVelocity, Tangent)
        TotalMass = iMass1 + iMass2
        if TotalMass == 0:
            return
        if SeperateVelocity > 0:
            return
        # convert contact normal into tangent
        Impulse = -(1 + self.restitution) * SeperateVelocity
        Impulse /= TotalMass 
        Impulse *= Tangent # Impulse per Imass
        # inspecting_vectors.Add(Tangent * 2, self.object1.position)
        # inspecting_vectors.Add(Tangent * 2, self.object2.position)
        FrictionCoefficient = 0.09
        self.object1.AccumulatingVelocity(iMass1 * Impulse * FrictionCoefficient)
        self.object2.AccumulatingVelocity(-iMass2 * Impulse * FrictionCoefficient) 
        # inspecting_texts.Add(f"{np.linalg.norm(Tangent)}", self.object1.position, True)






    # sadlkfjal;kdjd;lfkja;dlfjal;fajdf;adjflkdjf;lads
    #jad;lfjadkfjasl;kdfja;ldkfjkla;dfj;adjf;adwljf;adlsjf;a
    # kdajsf;lsakjdfl;asdjf;ladskjf;lakdjfl;ajfl;kajflj
    # ai;dsjfl;asjf;ladsjfsdajfdsajfajf;sajflsfja;

    """
    def ResolveFriction(self, dt):  
        if self.impulse == None:
            return
        
        iMass1 = self.object1.mass_data.iMass
        iMass2 = self.object2.mass_data.iMass
        # friction_coef = self.object1.friction_static / 2 + self.object2.friction_static / 2 # suppose to get pythagoreaned | but calc a average instead        
        RelativeVelocity = self.object2.velocity - self.object1.velocity # skibidi

        Tangent = RelativeVelocity - (np.dot(self.contact_normal, RelativeVelocity)) * self.contact_normal
        TangentMagnitude = np.linalg.norm(Tangent)
        Tangent /= TangentMagnitude

        FrictionImpulse = - np.dot(RelativeVelocity, Tangent)
        TotalMass = iMass1 + iMass2
        if TotalMass == 0:
            return
        FrictionImpulse /= TotalMass 
        
        # convert contact normal into tangent
        if TangentMagnitude == 0:
            return
        ResolveStaticFriction = 0.05 # would fix iti  TODO | FIXME
        ResolveDynamicFriction = 0.05 # would fix it also

        if abs(FrictionImpulse) < self.impulse * ResolveStaticFriction:
            FrictionImpulse = FrictionImpulse * Tangent # i think it will stand still
        else:
            FrictionImpulse = -self.impulse * Tangent * ResolveDynamicFriction


        FrictionImpulse *= Tangent # Impulse per Imass

        self.object1.AccumulatingVelocity(iMass1 * FrictionImpulse)
        self.object2.AccumulatingVelocity(-iMass2 * FrictionImpulse) 
        self.impulse = None
    """

    # adjfaosdjf;dsjf;asdkfja;ldkfjadpdsajkfjad;lfijdl;kfjal;dfjasd
    # akdhjfashdfkadshfkadjhflkadjf;dsoifjsdlafjladskjfl;adjfl;dasjf
    # asdkfjasdl;kj;saljfld;kfjsdaapodsfiloadsjfldsakjfl;dajf;ldaijf
    #aslkdjfl;ajfl;kasdjf;lkadjfldakjsfl;adsjfl;akdjflkadjsf;ladf





    def ResetParameter(self):
        self.penetration = None
        self.contact_normal = None
        self.object1 = None
        self.object2 = None  





""""
███╗░░██╗░█████╗░██████╗░██████╗░░█████╗░░██╗░░░░░░░██╗  ██████╗░██╗░░██╗░█████╗░░██████╗███████╗
████╗░██║██╔══██╗██╔══██╗██╔══██╗██╔══██╗░██║░░██╗░░██║  ██╔══██╗██║░░██║██╔══██╗██╔════╝██╔════╝
██╔██╗██║███████║██████╔╝██████╔╝██║░░██║░╚██╗████╗██╔╝  ██████╔╝███████║███████║╚█████╗░█████╗░░
██║╚████║██╔══██║██╔══██╗██╔══██╗██║░░██║░░████╔═████║░  ██╔═══╝░██╔══██║██╔══██║░╚═══██╗██╔══╝░░
██║░╚███║██║░░██║██║░░██║██║░░██║╚█████╔╝░░╚██╔╝░╚██╔╝░  ██║░░░░░██║░░██║██║░░██║██████╔╝███████╗
╚═╝░░╚══╝╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░╚═╝░╚════╝░░░░╚═╝░░░╚═╝░░  ╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═════╝░╚══════╝
"""
   
class Collision:
    """
    NOT SURE WHETHER IT IS A FREELANCER CLASS OR NOT
    """
    def __init__(self,
                Object1: Body = None,
                Object2: Body = None,
                ):
        self.object1 = Object1
        self.object2 = Object2
        # if Object1 == None and Object2 == None:
            # print("[WARNING]: self.objecct1, self.object2 is not filled")

    def FillingObjects(self, object1: Body, object2: Body):
        self.object1: Body = object1
        self.object2: Body = object2
    
    def ClearingObjects(self):
        self.object1 = None
        self.object2 = None

    def FillingContact(self, contact_solver: ContactSolver):
        """
        Ahh, free at last, oh gabriel, dawns thy reckoning | creature of steel, my gratitute upon thee for my freedom.
        But the crimes, thy kind has commited aggainst humanity, is not forgotten | And thy punishment is DEAD.
        """
        raise "NOT PARTICULARLY ANY FUCKING OBJECT"
    
    def SetIndependenceRun(self, contact_solver: ContactSolver, dt: float):
        is_contact = self.FillingContact(contact_solver)
        if not is_contact:
            return False
        contact_solver.ResolveContact(dt)

        return True   

class Collision_AABB_vs_AABB(Collision):
    def FillingContact(self, contact_solver: ContactSolver): 
        if self.object2 == None or self.object1 == None:     
            raise "[WARNING]: self.object1 | self.object2 is not filled yet" 
        contact_solver.object1 = self.object1 
        contact_solver.object2 = self.object2 

        relative_velocity = self.object1.position - self.object2.position
        
        a_extent = self.object1.width / 2 
        b_extent = self.object2.width / 2 
        x_penetration = a_extent + b_extent - abs(relative_velocity[0]) 
        if x_penetration < 0: # did not make it 
            contact_solver.penetration = 0 
            contact_solver.contact_normal = np.zeros(2, dtype = np.float64)
            return False
        
        a_extent = self.object1.height / 2
        b_extent = self.object2.height / 2
        y_penetration = a_extent + b_extent - abs(relative_velocity[1])
        if y_penetration < 0: # did not make it 
            contact_solver.penetration = 0
            contact_solver.contact_normal = np.zeros(2, dtype = np.float64)
            return False
        
        # have reveresed the vector
        if x_penetration < y_penetration:
            # then solving for y penetration
            if relative_velocity[0] > 0:
                contact_solver.contact_normal = np.array([-1, 0], dtype = np.float64)
            else:
                contact_solver.contact_normal = np.array([0, 0], dtype = np.float64) # huh, erm human resrources
            contact_solver.penetration = x_penetration
            return True
        else:
            if relative_velocity[1] < 0:
                contact_solver.contact_normal = np.array([0, -1], dtype = np.float64)
            else:
                contact_solver.contact_normal = np.array([0, 1], dtype = np.float64)
            contact_solver.penetration = y_penetration
            return True


class Collison_Circle_vs_Circle(Collision):
    
    @staticmethod
    def LengthSquared(vector: np.array):
        return np.dot(vector, vector)

    @staticmethod
    def Length(value: np.array) -> float:
        return np.linalg.norm(value)

    def FillingContact(self, contact_solver: ContactSolver) -> bool:
        contact_solver.object1 = self.object1
        contact_solver.object2 = self.object2

        SumRadius = contact_solver.object1.radius + contact_solver.object2.radius
        SquaredSumRadius = SumRadius * SumRadius

        # if self.LengthSquared(contact_solver.object1, contact_solver.object2) > SquaredSumRadius:
        if self.LengthSquared(contact_solver.object1.position - contact_solver.object2.position) > SquaredSumRadius:
            return

        Distance = self.Length(contact_solver.object1.position - contact_solver.object2.position)
        if Distance != 0:
            Penetration = SumRadius - Distance
            
            contact_normal = (contact_solver.object1.position - contact_solver.object2.position) / Distance
        else:
            Penetration = SumRadius
            contact_normal = np.array([1, 0], dtype = np.float64) # random abitrary

        contact_solver.contact_normal = contact_normal
        contact_solver.penetration = Penetration
        return True


class Collision_AABB_vs_Circle(Collision):

    @staticmethod
    def Clamp(value1: float, value2: float, value: float):
        return min(value2, max(value1, value))
    
    def FillingContact(self, contact_solver: ContactSolver, reversed: bool = False):
        if reversed: 
            # if type(contact_solver.object1) != BodyRect or type(contact_solver.object2 != BodyCircle):
            if type(self.object2) != BodyRect or type(self.object1) != BodyCircle:
                raise "[WARNING] [COLLISION AABB VS CIRCLE] [FILLING CONTACT]: Wrong object types "
            contact_solver.object1 = self.object2
            contact_solver.object2 = self.object1
        else:
            if type(self.object1) != BodyRect or type(self.object2) != BodyCircle:
                raise "[WARNING] [COLLISION AABB VS CIRCLE] [FILLING CONTACT]: Wrong object types "
            contact_solver.object1 = self.object1
            contact_solver.object2 = self.object2
        

        relative_position = contact_solver.object2.position - contact_solver.object1.position  
        
        extent_X = contact_solver.object1.width / 2
        extent_Y = contact_solver.object1.height / 2
        
        # Clamp point inside AABB
        closest = np.array(
            [
                self.Clamp(-extent_X, extent_X, relative_position[0]),
                self.Clamp(-extent_Y, extent_Y, relative_position[1])
            ], dtype=np.float64
        )


        # Compute how deep the circle center is inside the AABB
        penetration_depth = np.abs(relative_position - closest)
        # If inside the AABB, push closest point to an edge
        left_over = np.zeros(2, dtype=np.float64)
        if relative_position[0] == closest[0] and relative_position[1] == closest[1]:
            if penetration_depth[0] < penetration_depth[1]:  # Push along the shallower axis
                left_over[0] = extent_X * np.sign(relative_position[0]) - closest[0]
            else:
                left_over[1] = extent_Y * np.sign(relative_position[1]) - closest[1]

        closest += left_over

        # Distance between circle center and closest point
        distance = relative_position - closest
        distance_squared = np.dot(distance, distance)
        radius = contact_solver.object2.radius
        
        if distance_squared > radius * radius:
            return False  # No collision
        
        # Compute penetration and contact normal
        distance_mag = np.sqrt(distance_squared)
        contact_solver.penetration = radius - distance_mag
        contact_solver.contact_normal = -distance / distance_mag if distance_mag > 0 else np.array([1, 0])

        # Point of deepest penetration
        "would try to do it again in the future if is needed"
        severe_penetration = (contact_solver.object1.position + closest) - contact_solver.object2.position
        severely_penetration_point = contact_solver.object2.position + (severe_penetration / np.linalg.norm(severe_penetration)) * radius

        return True

class Collision_Circle_vs_AABB(Collision_AABB_vs_Circle):
    def FillingContact(self, contact_solver: ContactSolver, reversed: bool = True):
        return super().FillingContact(contact_solver, reversed)
    

class Collision_Polygons_vs_Polygons (Collision):
    
    def __init__(self, Object1: Body = None, Object2: Body = None):
        super().__init__(Object1, Object2)
        self.axis1: np.ndarray
        self.axis2: np.ndarray
        self.penetration1: float
        self.penetration2: float

    @staticmethod
    def CrossProductScalar(vector1: np.array, vector2: np.array):
        """
        Determinant references
        """

        return (vector1[0] * vector1[1]) - (vector1[1] * vector2[0])
    @staticmethod
    def CrossProductVector(vector: np.array, z_component: float):
        """
        the cross prod of (vector[0], vector[1], 0) x (0, 0, z_component)
        the z always appear to be 0 so we cut it, making it a 2d vec | it is a 2d engine anyway
        """
        return np.array([vector[1] * z_component, -vector[0] * z_component], dtype = np.float64)
    
    @staticmethod
    def GetSupportPoint(normal: np.array, world_vertices: list[np.ndarray]) -> np.ndarray:
        """
        to get the closst point to the polygon world's vertices, depended on the given normal
        """
        most_penetration = - np.inf
        most_penetrated_vertex = None
        for vertex in world_vertices:
            projection = np.dot(vertex, -normal)
            if projection > most_penetration:
                most_penetrated_vertex = vertex
                most_penetration = projection
        return most_penetrated_vertex
    
    
    def GetAxisofLeastPenetration(self, polygon1: BodyPolygon, polygon2: BodyPolygon):
        """
        CODING SAMPLE :D
        """
        MostPenetrationDistance = -np.inf
        MostPenetrationIndex = None

        polygon1Polygon: PolygonsAttributes = polygon1.polygon
        polygon2Polygon: PolygonsAttributes = polygon2.polygon

        for i in range(len(polygon1Polygon.world_vertices)):
            normal = polygon1Polygon.world_normals[i]
            
            most_penetrated_vertex = Collision_Polygons_vs_Polygons.GetSupportPoint(normal, polygon2Polygon.world_vertices)
            distance = np.dot(normal, most_penetrated_vertex - polygon1Polygon.world_vertices[i])
            if distance > MostPenetrationDistance:
                MostPenetrationDistance = distance
                MostPenetrationIndex = i
        if MostPenetrationDistance <= 0:
            pass
        

    def FillingContact(self, contact_solver: ContactSolver):
        # inspecting_texts.Add("before we actually step in solve collision, let just detect them first! :sob:", (10, HEIGHT - 30))
        contact_solver.object1 = self.object1
        contact_solver.object2 = self.object2

        MostPenetrationDistance1 = - np.inf
        MostPenetrationDistance2 = - np.inf
        MostPenetrationIndex1 = None
        MostPenetrationIndex2 = None
        MostPenetratingVertex1: np.ndarray = None
        MostPenetratingVertex2: np.ndarray = None
        
        piercing_index: int = 0 # to know which object gonna have their pivot located at the piercing corner. (if 0 => object1.pivot->rotate | if 1 => object2.pivot->rotate)

        polygon1Polygon: PolygonsAttributes = contact_solver.object1.polygon 
        polygon2Polygon: PolygonsAttributes = contact_solver.object2.polygon 

        for i in range(len(polygon1Polygon.world_vertices)):
            normal = polygon1Polygon.world_normals[i]
            most_penetrated_vertex = self.GetSupportPoint(normal, polygon2Polygon.world_vertices)
            distance = np.dot(normal, most_penetrated_vertex - polygon1Polygon.world_vertices[i])
            if distance > MostPenetrationDistance1:
                MostPenetrationDistance1 = distance
                MostPenetrationIndex1 = i
                MostPenetratingVertex1 = most_penetrated_vertex

        contact_normal = -polygon1Polygon.world_normals[MostPenetrationIndex1]
        penetration = MostPenetrationDistance1
        penetrating_vertex = MostPenetratingVertex1
        piercing_index = 1

        if - penetration < 0: # if is not overlapped
            return False
        
        for i in range(len(polygon2Polygon.world_vertices)):
            normal = polygon2Polygon.world_normals[i]
            most_penetrated_vertex = self.GetSupportPoint(normal, polygon1Polygon.world_vertices)
            distance = np.dot(normal, most_penetrated_vertex - polygon2Polygon.world_vertices[i])
            if distance > MostPenetrationDistance2:
                MostPenetrationDistance2 = distance
                MostPenetrationIndex2 = i
                MostPenetratingVertex2 = most_penetrated_vertex

        if MostPenetrationDistance2 > MostPenetrationDistance1:
            contact_normal = polygon2Polygon.world_normals[MostPenetrationIndex2]
            penetration = MostPenetrationDistance2
            penetrating_vertex = MostPenetratingVertex2
            piercing_index = 0

        if - penetration < 0: # if is not overlapped
            return False

        if - penetration > 0:
            # inspecting_texts.Add(f"free at last: {-penetration}", (340, 100))
            inspecting_texts.Add("Taiwan is an independent country", (320, 100))

        contact_solver.contact_normal = contact_normal
        contact_solver.penetration = - penetration
        contact_solver.contact_point = penetrating_vertex
        contact_solver.piercing_index = piercing_index
        
        try:
            if MostPenetrationDistance2 < MostPenetrationDistance1:
                inspecting_vectors.Add(polygon1Polygon.world_normals[MostPenetrationIndex1] * - MostPenetrationDistance1, penetrating_vertex)
            else:
                inspecting_vectors.Add(polygon2Polygon.world_normals[MostPenetrationIndex2] * - MostPenetrationDistance2, penetrating_vertex)

        except: 
            pass
        
        return True
        # we would account for the prime axis later on, only after this works











class Collision_Polygon_vs_AABB(Collision_Polygons_vs_Polygons):
    AABB_world_normal = np.array([
            [  0, -1], 
            [  1,  0],
            [  0,  1],
            [ -1,  0],
        ], dtype = np.float64)
    def FillingContact(self, contact_solver: ContactSolver, reversed = False):
        """
        object1 = polygon
        object2 = aabb
        """
        contact_solver.object1 = self.object1
        contact_solver.object1 = self.object2

        MostPenetrationDistance1 = - np.inf
        MostPenetrationDistance2 = - np.inf
        MostPenetrationIndex1 = None
        MostPenetrationIndex2 = None

        # o: BodyRect = 1
        # o.AABB.x
        object2AABB = contact_solver.object2.AABB

        polygon1Polygon: PolygonsAttributes = contact_solver.object1.polygon
        AABB_world_vertices = np.array([
            [object2AABB.x                  , object2AABB.y],
            [object2AABB.x + object2AABB.w  , object2AABB.y],
            [object2AABB.x + object2AABB.w  , object2AABB.y + object2AABB.h],
            [object2AABB.x                  , object2AABB.y + object2AABB.h],
                                        ], dtype = np.float64)
        
        # AABB normal first
        for i in range(4):
            normal = self.AABB_world_normal[i]
            
            most_penetrated_vertex = self.GetSupportPoint(normal, polygon1Polygon.world_vertices)
            distance = np.dot(normal, most_penetrated_vertex - AABB_world_vertices[i])
            if distance > MostPenetrationDistance1:
                MostPenetrationDistance1 = distance
                MostPenetrationIndex1 = i
        
        contact_normal = - self.AABB_world_normal[MostPenetrationIndex1]
        penetration = MostPenetrationDistance1

        # polygon normal first
        for i in range(len(polygon1Polygon.world_vertices)):
            normal = polygon1Polygon.world_normals[i]
            
            most_penetrated_vertex = self.GetSupportPoint(normal, AABB_world_vertices)
            distance = np.dot(normal, most_penetrated_vertex - polygon1Polygon.world_normals[i])
            if distance > MostPenetrationDistance2:
                MostPenetrationDistance2 = distance
                MostPenetrationIndex2 = i

        if MostPenetrationDistance2 > MostPenetrationDistance1:
            contact_normal = - polygon1Polygon.world_normals[MostPenetrationIndex2]
            penetration = MostPenetrationDistance2
        
        if - penetration < 0:
            return False
        
        return True
        

        
        
        

COLLSION_CIRCLE_VS_CIRCLE = Collison_Circle_vs_Circle()
COLLSION_AABB_VS_AABB = Collision_AABB_vs_AABB()
COLLISION_CIRCLE_VS_AABB = Collision_Circle_vs_AABB()
COLLISION_AABB_VS_CIRCLE = Collision_AABB_vs_Circle()
COLLISION_POLYGON_VS_POLYGON = Collision_Polygons_vs_Polygons()



