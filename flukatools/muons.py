#!/usr/bin/python

__version__ = 2.2
__author__ = 'Regan Ross'
## Last Edited Sept 27, 2022

'''
muons.py

Contact:
Regan Ross
regan.ross@mail.mcgill.ca

Tools for producing a phase space file for FLUKA simulations

'''
#################################################
#                    IMPORTS                    }
#                                               }
#################################################
from . import *
from scipy.optimize import minimize_scalar
#################################################
#                   CONSTANTS                   }
#             Referenced Throughout             }
#################################################

# Miscellaneous Information
# ior_water = 1.333       # water index of refraction
# c = 299792458           # [m/s]
# alpha = 7.297353e-3     # Fine structure constant
# R5912_min = 300e-9      # m     Min detectable wavelength for PMT
# R5912_max = 650e-9      # m     Max detectable wavelength for PMT
# mu_rest_mass_kg = 1.8835e-28   # [kg]
# mu_rest_mass_MeV = 105.66      # [MeV]
# elementary_charge = 1.60218e-19        # [C]
# #lXe_refractive_index = 

# # nEXO OUTER DETECTOR PARAMETERS
OD_RADIUS = 6.1722      # m
OD_HEIGHT = 12.800      # m
OD_CENTER = (0,0,0)     # m         defines coordinate system with respect to literal centre of OD
OC_RADIUS = 2.270       # m
ROI_RADIUS = OD_RADIUS + 2 # m
ROI_HEIGHT = OD_HEIGHT + 4 # m
GEN_OFFSET = OD_HEIGHT
GEN_RADIUS = np.tan(1)*(OD_HEIGHT + GEN_OFFSET) + OD_RADIUS

# OC_POSITION = (0,0,0.40) # m         positions OC with respect to OD centre
# TPC_RADIUS = 0.575      # m         from the pre-conceptual design report
# TPC_HEIGHT = 0.625      # m

# MUON FLUX PARAMETERS AT SNOLAB
SNOLAB_MU_FLUX = 3.31e-10       # \pm (0.01 (stat) \pm 0.09 (sys))e-10 mu/cm^2/s # arXiv:0902.2776v1
SNOLAB_MU_E_AVG = 363.0         # \pm 1.2 GeV # arXiv:1909.11728v1
SNOLAB_DEPTH = 5.890    #km.w.e        # \pm 94 km.w.e.  # arXiv:1909.11728v1

MIN_ENERGY = 0      #GeV
MAX_ENERGY = 25e3   #GeV



def mei_and_hime_eq3(zenith, vert_depth = SNOLAB_DEPTH, return_bounds = False):
    '''Literally just Mei and Hime's equation three.'''

    # Save the parameters with their statistical uncertainty

    I_1 = (8.60e-6, 5.3e-7)
    I_2 = (4.4e-7, 6e-8)
    lam_1 = (4.5e-1, 1e-2)
    lam_2 = (8.7e-1, 2e-2)

    mh_intensity = (I_1[0]*np.exp(-vert_depth/(lam_1[0] * np.cos(zenith))) + I_2[0]*np.exp(-vert_depth/(lam_2[0] * np.cos(zenith))))/np.cos(zenith)
    if not return_bounds:
        return mh_intensity
    
    else:
        # Determine the uncertainties in the equation
        # Use the largest values of I, and maximum values of lam
        maximum = ((I_1[0] + I_1[1])*np.exp(-vert_depth/((lam_1[0] + lam_1[1]) * np.cos(zenith))) + \
                (I_2[0] + I_2[1])*np.exp(-vert_depth/((lam_2[0] + lam_2[1]) * np.cos(zenith))))/np.cos(zenith)
        
        # Use the smallest values of I, smallest values of lam
        minimum = ((I_1[0] - I_1[1])*np.exp(-vert_depth/((lam_1[0] - lam_1[1]) * np.cos(zenith))) + \
                (I_2[0] - I_2[1])*np.exp(-vert_depth/((lam_2[0] - lam_2[1]) * np.cos(zenith))))/np.cos(zenith)

        return mh_intensity, maximum, minimum
    
def mei_and_hime_eq3_integrated(zenith, vert_depth = SNOLAB_DEPTH, return_bounds = False):
    '''Mei and Hime's equation three but multiplied by the differential solid angle!'''

    # Save the parameters with their statistical uncertainty

    I_1 = (8.60e-6, 5.3e-7)
    I_2 = (4.4e-7, 6e-8)
    lam_1 = (4.5e-1, 1e-2)
    lam_2 = (8.7e-1, 2e-2)

    mh_intensity = 2*np.pi*(I_1[0]*np.exp(-vert_depth/(lam_1[0] * np.cos(zenith))) + I_2[0]*np.exp(-vert_depth/(lam_2[0] * np.cos(zenith))))*np.tan(zenith)
    if not return_bounds:
        return mh_intensity
    
    else:
        # Determine the uncertainties in the equation
        # Use the largest values of I, and maximum values of lam
        maximum = 2*np.pi*((I_1[0] + I_1[1])*np.exp(-vert_depth/((lam_1[0] + lam_1[1]) * np.cos(zenith))) + \
                (I_2[0] + I_2[1])*np.exp(-vert_depth/((lam_2[0] + lam_2[1]) * np.cos(zenith))))*np.tan(zenith)
        
        # Use the smallest values of I, smallest values of lam
        minimum = 2*np.pi*((I_1[0] - I_1[1])*np.exp(-vert_depth/((lam_1[0] - lam_1[1]) * np.cos(zenith))) + \
                (I_2[0] - I_2[1])*np.exp(-vert_depth/((lam_2[0] - lam_2[1]) * np.cos(zenith))))*np.tan(zenith)

        return mh_intensity, maximum, minimum
    
def mei_hime_energy_intensity(zenith, energies):
    '''The energy intensity equation from Mei and Hime's equation '''
    
    b = 0.4 #km.w.e^{-1}
    epsilon = 693 #GeV
    gamma = 3.77
    h = (SNOLAB_DEPTH/np.cos(zenith))
    
    intensity = lambda E : np.exp(-b*h*(gamma - 1))*(E + epsilon*(1 - np.exp(-b*h)))**(-gamma)

    return intensity(energies)

def intersection_points(vec, labels = False, tolerance = 0.001):
    ''' A function for analytically determining the intersection points of a muon with an outer detector cylinder.  '''
    
    entryPoint, exitPoint = False, False
    entryLabel, exitLabel = '',''

    detRadius = OD_RADIUS + 2 #YAML_PARAMS['roi_radius']
    detHeight = OD_HEIGHT + 4 #YAML_PARAMS['roi_height']
    #det_z_translation = detector.position[2]

    #We can parametrize the muon for simplicity:
    zenith, azimuth = vec['zenith'], vec['azimuth']
    mx = np.sin(zenith)*np.cos(azimuth)
    my = np.sin(zenith)*np.sin(azimuth)
    mz = -np.cos(zenith)        # Does this need to be negative?

    x0, y0, z0 = vec['init_x'], vec['init_y'], vec['init_z'] 

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

def make_phase_space_file():

    how_many = YAML_PARAMS['num_muons']
    filename = PATHS['workdir'] + 'muons' + str(SEED) + '.txt'
    FLUKA_JOB_FILES['muons'] = filename

    ###########################################################################################################
    # Choose the FLUKA numbers of the particles (positive to negative ratio)

    pos_neg = (np.random.random(how_many) > 0.72).astype(int)
    pos_neg += 10
    
    ############################################################################################################
    # Sample a number of zenith angles from the Mei and Hime parameterization using rejection sampling

    max_angle = minimize_scalar(lambda theta : -mei_and_hime_eq3_integrated(theta))['x'] # The maximum value of the Mei and Hime zenith angle function in our range
    max_y = mei_and_hime_eq3_integrated(max_angle)

    over_sample = int(3*how_many) # Roughly only 1/3 will be kept based on the rejection sampling

    random_ys = np.random.random(over_sample)*max_y    # Choose a bunch of random "y" values

    random_angles = np.random.random(over_sample)*1.5  # Choose a bunch of random "\theta" values

    evald_ys = mei_and_hime_eq3_integrated(random_angles) # Determine the value of the function at the angles

    keep_indices = np.where(evald_ys > random_ys)   # Where the point lies within the function, keep the point

    zenith_angles = random_angles[keep_indices][:how_many]    # Keep only the desired number of angles.

    ############################################################################################################
    # Sample an energy for each zenith angle from the Mei and Hime energy parameterization (also using rejection sampling)
    max_y = mei_hime_energy_intensity(0,0)
    energies = np.ones(how_many)*-1
    init_x, init_y  = np.empty(how_many), np.empty(how_many)
    init_z = np.ones(how_many)*(GEN_OFFSET + OD_HEIGHT/2) # Required because the OD is centred on (0,0,0)
    azimuths = np.random.random(how_many)*np.pi*2


    for i in range(how_many):

        energy = energies[i]
        zenith = zenith_angles[i]

        while energy == -1:

            random_ys = np.random.random(300)*max_y    # Choose a bunch of random "y" values
            random_energies = np.random.random(300)*(MAX_ENERGY - MIN_ENERGY) + MIN_ENERGY
            evald_ys = mei_hime_energy_intensity(zenith, random_energies)
            keep_indices = np.where(evald_ys > random_ys)   # Where the point lies within the function, keep the point
            if np.size(keep_indices) == 0:
                energy = -1
            else:
                energy = random_energies[keep_indices][0]

        energies[i] = energy


        ############################################################################################################
        # Sample points on a disk and the remaining direction cosines so that the muon intersects the target region.

        intersection = False
        while intersection == False:
            # Keeping in mind the coordinate transformations 
            rho = np.random.random()*(GEN_RADIUS**2)
            polar = np.random.random()*np.pi*2
            initial_x = np.sqrt(rho)*np.cos(polar)
            initial_y = np.sqrt(rho)*np.sin(polar)
            vec = {'zenith' : zenith_angles[i],
                   'azimuth': azimuths[i], 
                   'init_x' : initial_x, 
                   'init_y' : initial_y,
                   'init_z' : init_z[i]}
            
            if type(intersection_points(vec) is not bool):
                intersection = True
        
        init_x[i]   = initial_x
        init_y[i]   = initial_y


    # pos_neg, energy, init_x, init_y, init_z, cos_x, cos_y, -cos_z, weight

    cos_z = np.cos(zenith_angles)
    cos_x = np.sin(zenith_angles)*np.cos(azimuths)
    cos_y = np.sin(zenith_angles)*np.sin(azimuths)

    phase_space = np.column_stack((pos_neg, energies, init_x, init_y, init_z, cos_x, cos_y, cos_z, weights))

    with open(filename, 'w') as phase_space_file:
        for i in range(how_many):
            # Use string substitution with format specifications
            formatted_row = '{} {} {} {} {} {} {} {} {}'.format(pos_neg[i], energies[i], init_x[i], init_y[i], init_z[i], cos_x[i], cos_y[i], cos_z[i], 1)
            phase_space_file.write(formatted_row + '\n')

            




