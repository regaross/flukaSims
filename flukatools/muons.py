#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'
## Last Edited Sept 27, 2022

'''
muons.py

Contact:
Regan Ross
regan.ross@mail.mcgill.ca

A module for simulating muons underneath SNOLAB's overburden for nEXO's Outer Detector.

'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################
from . import *
#################################################
#                   CONSTANTS                   }
#             Referenced Throughout             }
#################################################

# Miscellaneous Information
ior_water = 1.333       # water index of refraction
c = 299792458           # [m/s]
alpha = 7.297353e-3     # Fine structure constant
R5912_min = 300e-9      # m     Min detectable wavelength for PMT
R5912_max = 650e-9      # m     Max detectable wavelength for PMT
mu_rest_mass_kg = 1.8835e-28   # [kg]
mu_rest_mass_MeV = 105.66      # [MeV]
elementary_charge = 1.60218e-19        # [C]
#lXe_refractive_index = 

# nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 12.800      # m
OD_CENTER = (0,0,0)     # m         defines coordinate system with respect to literal centre of OD
OC_RADIUS = 2.270       # m
OC_POSITION = (0,0,0.40) # m         positions OC with respect to OD centre
TPC_RADIUS = 0.575      # m         from the pre-conceptual design report
TPC_HEIGHT = 0.625      # m

# MUON FLUX PARAMETERS AT SNOLAB
SNOLAB_MU_FLUX = 3.31e-10       # \pm (0.01 (stat) \pm 0.09 (sys))e-10 mu/cm^2/s # arXiv:0902.2776v1
SNOLAB_MU_E_AVG = 363.0         # \pm 1.2 GeV # arXiv:1909.11728v1
SNOLAB_DEPTH = 5.890    #km.w.e        # \pm 94 km.w.e.  # arXiv:1909.11728v1


#################################################
#                    CLASSES                     }
#                                                }
#################################################

class TPC:
    '''
    A class for the attributes of the cylindrical TPC. Very much like the OuterDetector, but filled with liquid Xenon.
    '''

    def __init__(self, radius=TPC_RADIUS, height=TPC_HEIGHT, vertical_depth=SNOLAB_DEPTH) -> 'TPC':
        ''' A basic constructor; see class docstring'''
        self.radius = radius
        self.height = height
        self.vertical_depth = vertical_depth


class OuterCryostat:
    '''
    A class for the spherical outer cryostat within nEXO's outer detector
        - radius                Sphere radius [meters]
        - center      *         Sphere center [meters, meters, meters]

        * Center positions the outer cryostat w.r.t the outer detector center.

        '''

    def __init__(self, radius=OC_RADIUS, position=OC_POSITION) -> 'OuterCryostat':
        '''A basic constructor; see class docstring'''
        self.radius = radius
        self.position = position


class OuterDetector:
    '''
    A class for nEXO's *cylindrical* outer detector.
        - radius                Cylinder radius [meters]
        - height                Cylinder height [meters]
        - center      *         Cylinder center [meters, meters, meters]
        - vertical depth        Attenuation by overburden [km.w.e]
        - fill_height   **      Vertical level to which tank is filled with water [meters]
        - outer_cryo            The outer cryostat within the detector - a sphere for path lengths

    * Subordinate components will be positioned w.r.t the center which is the literal center of the cylinder,
    NOT the center of the bottom of the cylinder.

    ** fill_height specifies to which vertical level the cylinder is filled with liquid. This number is used to
        define the COVER GAS region.

    '''

    def __init__(self, radius=OD_RADIUS, height=OD_HEIGHT, vertical_depth=SNOLAB_DEPTH,\
         fill_height=OD_HEIGHT-0.2, outer_cryo = OuterCryostat()) -> 'OuterDetector':
        ''' A basic constructor; see class docstring'''
        self.radius = radius
        self.height = height
        self.vertical_depth = vertical_depth
        self.fill_height = fill_height
        self.outer_cryo = outer_cryo

class Muon:
    '''
    A class for muons parameterized by empirical functions.

        - zenith                [rad]
        - azimuth               [rad]
        - energy                [GeV]
        - initial_position      [meters, meters, meters]
        - speed
        - path_length
        - impact_param
    '''
    
    # Attributes of the class:
    rest_mass_kg = 1.8835e-28   # [kg]
    rest_mass_MeV = 105.66      # [MeV]
    charge = 1.60218e-19        # [C]
    

    def __init__(self, zenith=0, azimuth=0, energy=SNOLAB_MU_E_AVG, initial=(0,0,0)) -> 'Muon':
        ''' A constructor for the muon. Defaults to vertical muon at average SNOLAB energy'''

        self.zenith = zenith
        self.azimuth = azimuth
        self.energy = energy
        self.initial = initial
        self.speed = c*np.sqrt(1-(self.rest_mass_MeV/(1000*self.energy + self.rest_mass_MeV))**2) # Relativistic Kinetic Energy
        self.path_length = 0
        self.impact_param = self.closest_approach((0,0,0))
        self.hits_cryostat = False
        self.pos_neg = np.random.random() < 0.72

        if self.pos_neg:
            self.fluka_number = 10 # Positive muon
        else:
            self.fluka_number = 11 # Negative muon


    ### Instance Functions

    def get_unit_vec(self) -> tuple:
        ''' Returns a unit vector for the muon direction (x,y,z)'''

        x = np.sin(self.zenith)*np.cos(self.azimuth)
        y = np.sin(self.zenith)*np.sin(self.azimuth)
        z = -np.cos(self.zenith)

        return (x,y,z)

    def get_cartesian_track(self) -> tuple:
        ''' Returns the muon track (x, y, z, x0, y0, z0)'''

        return self.get_unit_vec()+self.initial

    def closest_approach(self, point) -> float:
        ''' Finds the muon's closest distance to the provided point in 3-space (cartesian coords)'''

        unit_vector = self.get_unit_vec()

        # Parameterize for simplicity

        i, j, k = point[0], point[1], point[2]
        mx, my, mz = unit_vector[0], unit_vector[1], unit_vector[2]
        x0, y0, z0 = self.initial[0], self.initial[1], self.initial[2]

        t = ((mx*i + my*j + mz*k) - (mx*x0 + my*y0 + mz*z0))/(mx**2 + my**2 + mz**2)

        muon_point = (mx*t + x0, my*t + y0, mz*t + z0)

        x_dist = muon_point[0] - i
        y_dist = muon_point[1] - j
        z_dist = muon_point[2] - k

        distance = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)

        return distance

    def time_at_distance(self, dist):
        '''Returns the amount of time since the muon was instantiated based on the distance of the muon along its track. Assumes constant speed.'''
        return (dist/self.speed)

    def direction_cosines(self):
        '''returns a tuple of the direction cosines of the muon'''
        cos_z = np.cos(self.zenith)
        cos_x = np.sin(self.zenith)*np.cos(self.azimuth)
        cos_y = np.sin(self.zenith)*np.sin(self.azimuth)

        return (cos_x, cos_y, cos_z)

    def __str__(self):
        '''Returns a string of the particle fit for a phase space file for a FLUKA Source'''

        cos_x, cos_y, cos_z = self.direction_cosines()

        energy = self.energy
        initial = self.initial
        weight = 1
        fnumber = self.fluka_number

        return '{} {} {} {} {} {} {} {} {}'.format(fnumber, energy, initial[0], initial[1], initial[2], cos_x, cos_y, -cos_z, weight)




#################################################
#                   FUNCTIONS                    }
#                                                }
#################################################

def set_seed(seed):
    '''Sets the numpy random generation seed with a given seed'''
    np.random.seed(seed)

def mei_hime_intensity(zenith_angles, vert_depth = SNOLAB_DEPTH)-> np.ndarray:
    ''' Function from Mei & Hime's zenith angle intensity relation- exactly the same as Eqn. 3 from the paper.'''
    # Parameters same as in paper (Equation No. 3)
    I1 = 8.60e-6 #  /sec/cm^2/sr
    I2 = 0.44e-6 #  /sec/cm^2/sr
    lam1 = 0.45 #  km.w.e.
    lam2 = 0.87 #  km.w.e.

    intensity = 2*np.pi*(I1*np.exp(-vert_depth/(lam1*np.cos(zenith_angles)))+I2*np.exp(-vert_depth/\
        (lam2*np.cos(zenith_angles))))/np.cos(zenith_angles)

    return intensity

def mei_hime_normed_discrete(zenith_angles, vert_depth = SNOLAB_DEPTH) -> np.ndarray:
    ''' Returns the normalized through-going muon flux through horizontal surface for an array of angles 
        "angles" at a specified vertical depth vert_depth in km.w.e. NORMALIZED for the given array of angles
        by making array sum to 1. Differential solid angle factor included.'''

    normed_flux = mei_hime_intensity(zenith_angles)* np.cos(zenith_angles)*np.sin(zenith_angles)  # Horizonal projection * differential solid angle

    #Normalizes the distribution function such that array sums to 1 for any number of elements
    normed_flux = normed_flux / np.sqrt(normed_flux.sum()**2)

    return normed_flux

def mei_hime_normed_continuous(zenith_angles, vert_depth = SNOLAB_DEPTH):
    ''' Returns the normalized through-going muon flux through horizontal surface 
        for an array of angles "angles" at a specified vertical depth vert_depth in km.w.e
        NORMALIZED based on the integral, not the discrete array'''

    norm_const = 6.006989507403272e-11 #Integrated from 0 to pi/2

    #Normalizes the distribution function
    normed_flux = mei_hime_intensity(zenith_angles)*np.cos(zenith_angles)*np.sin(zenith_angles)

    return normed_flux(zenith_angles)/norm_const


def mh_energy_probs(energies, zenith = 0, sample = True):
    from scipy.integrate import quad
    
    b = 0.4 #km.w.e^{-1}
    epsilon = 693 #GeV
    gamma = 3.77
    h = (SNOLAB_DEPTH/np.cos(zenith))
    
    prob = lambda E : np.exp(-b*h*(gamma - 1))*(E + epsilon*(1 - np.exp(-b*h)))**(-gamma)
    norm_const = quad(prob, 0, np.pi/2)[0]
    if not sample:
        array = prob(energies)/norm_const
    else:
        array = prob(energies)/np.sum(prob(energies))

    return array


def gaisser_normed_discrete(energies, zenith):
    ''' Surface muon energy distribution based on Gaisser's formalism '''

    # Making Gaisser's distribution as a lambda function [cm^{-2}s^{-1}sr^{-1}GeV^{-1}]
    dNdE = lambda E, theta : 0.14*(E**(-2.7))*((1/(1+(1.1*E*np.cos(theta)/115)))+(0.054/(1+(1.1*E*np.cos(theta)/850))))

    energy_array = dNdE(energies, zenith)

    # Normalize it for use as PDF
    return energy_array/np.sum(energy_array)

def get_disk_radius(detector, gen_offset = 0):
    if gen_offset == 0:
        gen_offset = detector.height

    gen_radius = np.tan(1)*(detector.height + gen_offset) + detector.radius
    return gen_radius


def generate_muons(how_many, detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
    ''' Generates an array of Muons (instances of class) given an OuterDetector cylinder
        - how_many          The number of muons to generate
        - outer_detector    The OuterDetector volume to be used in simulations
        - gen_radius        The radius of the concentric disk on which muons will be instantiated
        - gen_offset        The distance from the top of the outer_detector to the concentric disk
        '''

    muons = np.empty(how_many, dtype=object)

    sampling_size = int(10*how_many)

    # If the user has not set the generator offset and or radius
    if gen_offset == 0:
        gen_offset = detector.height

    if gen_radius == 0:
            gen_radius = np.tan(1)*(detector.height + gen_offset) + detector.radius

    # Define muon initial positions
    rhos = np.random.random(size = how_many)*(gen_radius**2)
    gen_angles = np.random.random(size = how_many)*np.pi*2

    # Keeping in mind the coordinate transformations 
    initial_x = np.sqrt(rhos)*np.cos(gen_angles)
    initial_y = np.sqrt(rhos)*np.sin(gen_angles)
    initial_z = np.ones(how_many)*(gen_offset + detector.height/2)

    # Determining zenith and azimuthal angles
    theta_radians = np.linspace(0, np.pi/2, sampling_size)
    zenith_probabilities = mei_hime_normed_discrete(theta_radians, detector.vertical_depth)
    zeniths = np.random.choice(theta_radians, p = zenith_probabilities, size = how_many)
    
    azimuths = np.random.random(size = how_many)*np.pi*2

    # Selecting Energies
    surface_energy_range = np.linspace(5000, 100000, sampling_size)
    surface_energies = np.random.choice(surface_energy_range, p = gaisser_normed_discrete(surface_energy_range, 0), size = how_many)

    # Attenuate based on zenith angles
    b = 0.4 #km.w.e.^-1
    energies_underground = (surface_energies - SNOLAB_MU_E_AVG*(np.exp(b*np.cos(zeniths)*SNOLAB_DEPTH)-1))/np.exp(b*np.cos(zeniths)*SNOLAB_DEPTH)
    for e in range(len(energies_underground)): 
        if energies_underground[e] < 0: energies_underground[e] = 0
    
    
    # Instantiate muons
    for i in range(how_many):
        muons[i] = Muon(zeniths[i], azimuths[i], energies_underground[i], (initial_x[i], initial_y[i], initial_z[i]))

    return muons

def muons_per_square_meter(rate, detector, z_offset = 0, gen_radius = 0):
    '''Returns an array of intersecting muons evenly distributed over the generator area at the specified concentration'''
    
    gen_offset = z_offset #outer_detector.height
    
    if gen_radius == 0:
        gen_radius = np.tan(1)*((detector.height/2) + gen_offset + detector.height/2)\
            + detector.radius
        
    gen_area = np.pi*(gen_radius)**2 

    n_muons = int(gen_area*rate)

    muons = generate_muons(n_muons, detector, gen_radius, gen_offset)

    remove_indices = []

    for i in range(len(muons)):
        if not hits_detector(muons[i], detector):
            remove_indices.append(i)
    
    return np.delete(muons, remove_indices, axis = 0)

def intersection_points(muon, detector = OuterDetector(), labels = True, tolerance = 0.001):
    ''' A function for analytically determining the intersection points of a muon with an outer detector cylinder.  '''
    
    entryPoint, exitPoint = False, False
    entryLabel, exitLabel = '',''

    detRadius = detector.radius
    detHeight = detector.height
    #det_z_translation = detector.position[2]

    #We can parametrize the muon for simplicity:
    zenith, azimuth = muon.zenith, muon.azimuth
    mx = np.sin(zenith)*np.cos(azimuth)
    my = np.sin(zenith)*np.sin(azimuth)
    mz = -np.cos(zenith)        # Does this need to be negative?

    x0, y0, z0 = muon.initial[0], muon.initial[1], muon.initial[2] 

    # Quadratic equation parameters
    a = (mx**2 + my**2)
    b = 2*(mx*x0 + my*y0)
    c = (x0**2 + y0**2 - detRadius**2)
    det_squared = (b**2 - 4*a*c)
    
    # Initial setting
    qhighCheck = detRadius + 1
    qlowCheck = qhighCheck
    
    if det_squared > 0:
        qlow = float(-b/(2*a) - np.sqrt(det_squared)/(2*a))
        zlow = mz*qlow + z0
        qlowCheck = a*qlow**2 + b*qlow + (x0**2 + y0**2)
        qhigh = float(-b/(2*a) + np.sqrt(det_squared)/(2*a))
        zhigh = mz*qhigh + z0
        qhighCheck = a*qhigh**2 + b*qhigh + (x0**2 + y0**2)


    # Find the first point
    ptop = (detHeight/2 - z0)/mz
    xtop = mx*ptop + x0
    ytop = my*ptop + y0

    # ENTRY POINT POSSIBILITIES
    entryLabel = ''
    if xtop**2 + ytop**2 <= detRadius**2:
        # Hits top
        entryPoint = (xtop, ytop, detHeight/2)
        entryLabel = 'TOP'
    elif qhighCheck < (detRadius + tolerance)**2 and qhighCheck > (detRadius - tolerance)**2 and zhigh**2 < (detHeight/2)**2:
        # Hits the side at the higher z point
        entryPoint = (mx*qhigh + x0, my*qhigh + y0, zhigh)
        entryLabel = 'SIDE'
        qhigh = False
    elif qlowCheck < (detRadius + tolerance)**2 and qlowCheck > (detRadius - tolerance)**2 and zlow**2 < (detHeight/2)**2:
        # Hits the side at the lower z point
        entryPoint = (mx*qlow + x0, my*qlow + y0, zlow)
        entryLabel = 'SIDE'
        qlow = False

    if type(entryPoint) is tuple: #If it does actually enter the cylinder
        # Bottom Point Parameters
        pbottom = (-detHeight/2 - z0)/mz
        xbottom = mx*pbottom + x0
        ybottom = my*pbottom + y0

        exitLabel = ''
        # EXIT POINT POSSIBILITIES
        if xbottom**2 + ybottom**2 <= detRadius**2:
            # Hits the bottom of the cylinder
            exitPoint = (xbottom, ybottom, -detHeight/2)
            exitLabel = 'BOT'

        elif qhighCheck < (detRadius + tolerance)**2 and qhighCheck > (detRadius - tolerance)**2 and zhigh**2 < (detHeight/2)**2 and type(qhigh) is float:
            exitPoint = (mx*qhigh + x0, my*qhigh + y0, zhigh)
            exitLabel = 'SIDE'

        elif qlowCheck < (detRadius + tolerance)**2 and qlowCheck > (detRadius - tolerance)**2 and zlow**2 < (detHeight/2)**2 and type(qlow) is float:
            exitPoint = (mx*qlow + x0, my*qlow + y0, zlow)
            exitLabel = 'SIDE'
    
    if type(entryPoint) is bool:
        return False

    elif labels:
        return (entryPoint, entryLabel, exitPoint, exitLabel)
    else:
        return (entryPoint, exitPoint)


# def path_length(muon, outer_detector = OuterDetector(), cryostat = False):
#     ''' Returns the path length of a muon through the Outer Detector. False if it doesn't hit.'''

#     points = intersection_points(muon, outer_detector, labels = False)
#     path_length = 0

#     if type(points) is not bool:
#         x = points[1][0] - points[0][0]
#         y = points[1][1] - points[0][1]
#         z = points[1][2] - points[0][2]

#         path_length = np.sqrt(x**2 + y**2 + z**2)

#         return path_length

#     else:
#         return False
    
    
def intersecting_muons(how_many, detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
    ''' Does the same as generate_muons, but returns only muons that intersect the provided outer detector'''

    at_a_time = int(how_many/2)+1
    muon_list = []

    while len(muon_list) < how_many:
        temp_muons = generate_muons(at_a_time, detector, gen_radius, gen_offset)
        for mu in temp_muons:
            if hits_detector(mu, detector):
                muon_list.append(mu)
                mu.path_length = path_length(mu, detector, labels=False)

    return np.array(muon_list)[:how_many]

def intersecting_muons_with_time(how_many, detector = OuterDetector(), gen_radius=0, gen_offset=0) -> tuple:
    ''' Returns a tuple of an array with the intersecting muons and also the count of muons that had to be tested to produce the array'''

    at_a_time = int(np.sqrt(how_many))
    no_hits_count = 0
    hits_count = 0

    muon_list = []

    while len(muon_list) < how_many:
        temp_muons = generate_muons(at_a_time, detector, gen_radius, gen_offset)
        for mu in temp_muons:
            if hits_detector(mu, detector):
                if hits_count == how_many:
                    break
                hits_count += 1
                muon_list.append(mu)
                mu.path_length = path_length(mu, detector, labels=False)
            else:
                no_hits_count += 1

    total_muons = no_hits_count + hits_count
    area = (10**4)*np.pi*get_disk_radius(detector, detector.height)**2 ## cm squared

    seconds = (total_muons/area)/SNOLAB_MU_FLUX ## seconds
    hours = seconds / 3600

    return (np.array(muon_list)[:how_many], hours)

def get_gen_area(detector = OuterDetector(), cm2 = False) -> float:
    gen_radius = get_disk_radius(detector)
    area = np.pi*gen_radius**2
    if not cm2:
        return area # m^2
    else:
        return area*1e4

def muons_from_time(hours, detector = OuterDetector(), intersecting = True, gen_radius=0, gen_offset=0) -> tuple:
    ''' Returns a tuple of an array with the intersecting muons and also the count of muons that had to be tested to produce the array'''

    number_of_muons = int(hours*3600*SNOLAB_MU_FLUX*get_gen_area(detector, cm2=True)) + 1 ## plus one for conservative estimation
    muons = generate_muons(number_of_muons, detector)
    if intersecting:
        muon_list = [muon for muon in muons if hits_detector(muon, detector)]
        return muon_list
    else:
        return muons
    
def non_intersecting_muons(how_many, detector = OuterDetector(), gen_radius=0, gen_offset=0) -> np.ndarray:
    ''' Does the same as generate_muons, but returns only muons that do not intersect the provided outer detector'''

    at_a_time = int(how_many/2)+1
    muon_list = []

    while len(muon_list) < how_many:
        temp_muons = generate_muons(at_a_time, detector, gen_radius, gen_offset)
        for mu in temp_muons:
            if not hits_detector(mu, detector):
                muon_list.append(mu)
                mu.path_length = path_length(mu, detector, labels=False)

    return np.array(muon_list)[:how_many]

def hits_detector(muon, detector) -> bool:
    ''' Returns True if the muon intersects the detector, False otherwise'''
    return (type(intersection_points(muon, detector, labels= False)) is not bool)


def path_length(muon, detector = OuterDetector(), labels=True, ignore_cover_gas=True, ignore_cryostat=True):
    ''' Returns the path length of a muon through the Outer Detector. False if it doesn't hit.'''

    points = intersection_points(muon, detector)
    path_length = 0

    if type(points) is not bool:
        x = points[2][0] - points[0][0]
        y = points[2][1] - points[0][1]
        z = points[2][2] - points[0][2]

        path_length = np.sqrt(x**2 + y**2 + z**2)


        if ignore_cover_gas:
            z_fill_line = (detector.fill_height - detector.height/2)
            if points[0][2] > z_fill_line:
                # If we are to subtract the length of the path due to cover gas:
                path_length = path_length - path_through_covergas(muon, detector)
            
        if ignore_cryostat: path_length = path_length - path_through_cryostat(muon, detector)

        if path_length < 0: path_length = 0

        if not labels: return path_length
        else: return (path_length, points[1], points[3])

    else:
        return False

def path_through_covergas(muon, outer_detector) -> float:
    ''' Returns the pathlength of the muon through the outer_detector cover gas layer'''
    
    if not hits_detector(muon, outer_detector): return False
    
    mx, my, mz = muon.get_unit_vec()
    x0, y0, z0 = muon.initial

    z_fill_line = (outer_detector.fill_height - outer_detector.height/2)
    
    points = intersection_points(muon, outer_detector)

    if points[0][2] > z_fill_line:
        # If it does actually hit the cover gas region:
        if points[2][2] > z_fill_line: return 0 # If the muon enters and exits through the cover gas

        else:
            z_path = points[0][2] - z_fill_line
            t = (z_fill_line - z0)/mz
            x_exit = mx*t + x0
            y_exit = my*t + y0

            x_path = points[0][0] - x_exit
            y_path = points[0][1] - y_exit

            path_through_covergas = np.sqrt(x_path**2 + y_path**2 + z_path**2)

    return path_through_covergas

def path_through_cryostat(muon, outer_detector) -> float:
    ''' Returns the pathlength through the Outer Detector without the contribution of the cryostat'''

    if not hits_detector(muon, outer_detector): return False

    oc = outer_detector.outer_cryo
    oc_radius = oc.radius

    mx, my, mz = muon.get_unit_vec()
    x0, y0, z0 = muon.initial
    ocx, ocy, ocz = oc.position

    A_x, B_x, C_x = ((mx**2),(2*mx*(x0 - ocx)),(x0 - ocx)**2)
    A_y, B_y, C_y = ((my**2),(2*my*(y0 - ocy)),(y0 - ocy)**2)
    A_z, B_z, C_z = ((mz**2),(2*mz*(z0 - ocz)),(z0 - ocz)**2)

    A = A_x + A_y + A_z
    B = B_x + B_y + B_z
    C = C_x + C_y + C_z - oc_radius**2

    det_squared = B**2 - 4*A*C

    if det_squared < 0:
        # Muon doesn't go through Cryostat
        return 0

    # If the muon does, in fact, hit the OC...
    muon.hits_cryostat = True

    t_one = (-B + np.sqrt(det_squared))/2
    t_two = (-B - np.sqrt(det_squared))/2

    point_one = np.array([mx*t_one + x0, my*t_one + y0, mx*t_one + z0])
    point_two = np.array([mx*t_two + x0, my*t_two + y0, mx*t_two + z0])
    diff = point_one - point_two

    dist = np.sqrt(np.sum(diff**2))

    return dist


def get_cherenkov(muon, detector = OuterDetector(), photons_per_meter = False) -> np.ndarray:
    ''' Returns a cherenkov light cone for the provided muon through the provided detector'''
    speed = muon.speed # Relativistic Kinetic Energy
    beta = speed/c
    light_angle = np.arccos(1/(beta*ior_water))
    N = 2*alpha*np.pi*((1/R5912_min)-(1/R5912_max))*(1-(1/(beta**2*ior_water**2))) # Photons per meter

    if not photons_per_meter:
        total_photons = int(N*path_length(muon, detector, labels = False, ignore_cover_gas = True, ignore_cryostat = True))
    else:
        total_photons = N

    return np.array([int(total_photons), light_angle])

def make_phase_space_file() -> tuple:

    '''A function to make a phase space file for a number of muons. 
        This file is used by the larger simulation as the source for each particle.
        NOTE: if the filename is changed, it must also be changed in the .inp file for the sim.
        
        Muon initial units are to be METERS. This is converted to FLUKA native cm in the read_phase_space_file routine'''
    
    filename = PATHS['workdir'] + 'muons' + str(SEED) + '.txt'
    FLUKA_JOB_FILES['muons'] = filename

    num_muons = YAML_PARAMS['num_muons']
    roi_radius = YAML_PARAMS['roi_radius']
    roi_height = YAML_PARAMS['roi_height']
    intersecting = YAML_PARAMS['intersecting']

    file_stream = open(filename, 'w')

    if roi_radius <= 0:
            roi_radius = OD_RADIUS
            roi_height = OD_HEIGHT
    else:
        roi_radius = roi_radius + OD_RADIUS
        roi_height = roi_height + OD_HEIGHT

    roi = OuterDetector(roi_radius, roi_height)
    
    if intersecting:
        muarray, hours = intersecting_muons_with_time(num_muons, roi)
    else:
        muarray = non_intersecting_muons(num_muons, roi)


    for muon in muarray:
        file_stream.write(str(muon) + '\n')

    file_stream.close()

    return muarray, hours







# #################################################
# #                   PLOTTING                     }
# #                  FUNCTIONS                     }
# #################################################

# def plot_path_lengths(points, savefile = False):
#     ''' A function designed to plot a histogram of the points from the previous intersection_points function. Pass this function an array
#     argument equivalent to the output from the previous.
#     '''
#     import matplotlib.pyplot as plt
#     plt.rcParams.update({
#     "text.usetex": True,})

#     pathLengths = []
#     top_bottom = [] #Muons that pass through both top and bottom
#     top_side = []
#     side_side = []
#     ssCount = 0
#     side_bottom = []

#     ##(entryPoint, entryLabel, exitPoint, exitLabel)
#     for point in points:
#         if not type(point) is bool: #If there is an entry point
#             if type(point[1]) is str:
#                 pathLength = np.sqrt((point[0][0]-point[2][0])**2 + (point[0][1] - point[2][1])**2 + (point[0][2] - point[2][2])**2)
#                 #If points are both top and bottom, append there
#                 if point[1] == 'TOP' and point[3] == 'BOT':
#                     top_bottom.append(pathLength)
#                 #If points are top side, append there
#                 elif point[1] == 'TOP' and point[3] == 'SIDE':
#                     top_side.append(pathLength)
#                 #If points are side side, append there
#                 elif point[1] == 'SIDE' and point[3] == 'SIDE':
#                     side_side.append(pathLength)
#                     ssCount += 1
#                 #If points are side bottom, append there
#                 elif point[1] == 'SIDE' and point[3] == 'BOT':
#                     side_bottom.append(pathLength)
#                 pathLengths.append(pathLength)
#             else:
#                 print(' Points input for path_lengths function requires labels being switched on ')

#     bins = 100
#     counts_tb, bins_tb = np.histogram(top_bottom, bins = bins, density=False)
#     counts_ts, bins_ts = np.histogram(top_side, bins = bins, density=False)
#     counts_ss, bins_ss = np.histogram(side_side, bins = bins, density=False)
#     counts_sb, bins_sb = np.histogram(side_bottom, bins = bins, density=False)

#     counts, bins = np.histogram(pathLengths, bins = 50, density=False)
#     plt.figure()
#     plt.hist(bins[:-1], bins, weights=counts, histtype='stepfilled', alpha=0.4, color='orange', label = 'Total')
#     plt.xlabel('Path Length [m]', size = 'large'); plt.ylabel('Count', size = 'large')
#     #plt.text(1, counts[0]*1.5, 'Mean = ' + str(np.average(pathLengths)), size = 12)
#     #plt.title('Path Length Distribution: '+ str(len(pathLengths)) + ' paths', size = 'x-large')

#     #Plotting other subordinate hists
#     top_bot = str(len(top_bottom)*100/len(pathLengths))[0:4]
#     plt.hist(bins_tb[:-1], bins, weights=counts_tb, histtype='step', alpha=1.0, color='blue', label = r'Top $\rightarrow$ Bottom: ' + top_bot + '\%')
#     top_sides = str(len(top_side)*100/len(pathLengths))[0:4]
#     plt.hist(bins_ts[:-1], bins, weights=counts_ts, histtype='step', alpha=1.0, color='green', label = r'Top $\rightarrow$ Side: ' + top_sides + '\%')
#     sides = str(len(side_side)*100/len(pathLengths))[0:4]
#     plt.hist(bins_ss[:-1], bins, weights=counts_ss, histtype='step', alpha=1.0, color='red', label = r'Side $\rightarrow$ Side: ' + sides + '\%')
#     side_bot = str(len(side_bottom)*100/len(pathLengths))[0:4]
#     plt.hist(bins_sb[:-1], bins, weights=counts_sb, histtype='step', alpha=1.0, color='purple', label = r'Side $\rightarrow$ Bottom: ' + side_bot + '\%')

#     plt.yscale('log')
#     plt.legend(loc = 8, fontsize = 'large')
#     #plt.grid()
#     #print('Mean:', np.mean(pathLengths))
#     if savefile:
#         plt.savefig('pathLengths.png', facecolor = 'white')

#     plt.show()

# def muons_to_array(muons):
#     '''Converts an array of muon objects to a 2D array of the muon attributes'''

#     # [Zenith, Azimuth, Energy, Initial_x, Initial_y, Initial_z, impact_parameter]
#     rows = len(muons)
#     cols = 7

#     new_array = np.ndarray((rows, cols), dtype=float)

#     iter = 0
#     for muon in muons:
#         temp_array = np.array([muon.zenith, muon.azimuth, muon.energy, muon.initial[0], muon.initial[1], muon.initial[2], muon.impact_param])
#         new_array[iter] = temp_array
#         iter += 1


#     return new_array

# def plot_muon_counts_hist(days, od = OuterDetector(), with_oc = True, oc = OuterCryostat(), with_tpc = True, tpc = OuterDetector(0.6, 1.35)):
#     ''' Plots overlapping histograms of counts of muons hitting the OD, Outer Cryostat and TPC.'''
#     import matplotlib.pyplot as plt

#     od_counts = []
#     oc_counts = []
#     tpc_counts = []
    
#     for day in range(days):
#         muons = muons_from_time(24)
#         od_count = 0
#         oc_count = 0
#         tpc_count = 0

#         for muon in muons:
#             if hits_detector(muon, od):
#                 od_count += 1
            
#             if path_through_cryostat(muon, od) > 0:
#                 oc_count += 1

#             if hits_detector(muon, tpc):
#                 tpc_count += 1

#         od_counts.append(od_count)
#         oc_counts.append(oc_count)
#         tpc_counts.append(tpc_count)

#     od_counts = np.array(od_counts)
#     oc_counts = np.array(oc_counts)
#     tpc_counts = np.array(tpc_counts)

#     od_hist = plt.hist(od_counts, histtype = 'stepfilled', bins = (od_counts.max() - od_counts.min()), label = r'OD $\mu = $ ' + str(np.mean(od_counts))[:4], density=True, alpha = 0.7)
#     oc_hist = plt.hist(oc_counts, histtype = 'stepfilled', bins = (oc_counts.max() - oc_counts.min()), label = r'OC $\mu = $ ' + str(np.mean(oc_counts))[:4], density=True, alpha = 0.7)
#     tpc_hist = plt.hist(tpc_counts, histtype = 'stepfilled', bins = (tpc_counts.max() - tpc_counts.min()), label = r'TPC $\mu = $ ' + str(np.mean(tpc_counts))[:4], density=True, alpha = 0.7)

#     years = days/365

#     plt.xlabel('Muons counted per day')
#     plt.ylabel('Fraction of Runtime')
#     plt.title('Muon counts in detector for ' + str(years) + ' years')
#     plt.legend()
