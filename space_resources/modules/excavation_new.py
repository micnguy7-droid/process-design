"""

author: DD

Version 1.0

Testing mod

"""



import math as ma

from modules.transportation_onlyBeta import * # To import



import numpy as np



# Excavation parameters

depthM = 0.025  # m

trenchDepthM = 0.1  # m

radiusM = 0.15  # m

extAngle = 10  # deg

intAngle = math.degrees(phiVal) #45  # deg

cohCoeff = kcVal  # 2100  # Pa





def excavationMechanics(depth_m, trench_depth_m, radius_m, ext_angle, int_angle, coh_coeff, gVar, mRegolithVar):

    phi = int_angle*ma.pi/180   # Internal fiction angle (rad)

    delta = ext_angle*ma.pi/180   # External friction angle (rad)

    g = gVar  # 1.622   # Lunar gravity, m/s^2

    rho_N = 1800    # Max regolith density, kg/m^3

    rho_0 = 1100    # Min regolith density, kg/m^3

    H = 0.06   # Constant parameter for density function

    gamma = rho_N - (rho_N - rho_0) * \
        ma.exp(-(trench_depth_m+depth_m)/H)  # density calculation

    C = coh_coeff  # Cohesion coefficient (Pa)

    reg_vars = [phi, delta, gamma, C]

    reg_vars2 = g    # Lunar gravity, m/s^2



    r = (3.125*(1-ma.cos(62/180*ma.pi))+5.3125) / \
        5.3125*radius_m  # Scoop radius, m

    w = 0.100838  # Bucket drum width, m

    alpha_b = 90*ma.pi/180    # Leading edge angle, rad

    beta = 5*ma.pi/180        # Rake angle, rad

    s = 0.0032                # Side thickness, m

    e_b = 0.0015875           # Leading edge thickness, m

    exc_vars = [radius_m, r, depth_m]

    exc_vars2 = [w, alpha_b, beta, s, e_b]



    massReg = mRegolithVar  # 90  # Mass of regolith, kg, max of excavator



    # Excavation force calculation, N

    F_dig = Balovnev(reg_vars, reg_vars2, exc_vars, exc_vars2)

    global digDistance

    digDistance = massReg/(gamma*w*depth_m)

    E_dig = F_dig*digDistance/massReg  # Excavation energy per kg regolith

    dig_outputs = [F_dig, E_dig]

    return dig_outputs





def Balovnev(reg_vars, reg_vars2, exc_vars, exc_vars2):

    f = reg_vars[0]

    d = reg_vars[1]

    gam = reg_vars[2]

    c = reg_vars[3]

    g = reg_vars2



    R = exc_vars[0]

    r = exc_vars[1]

    depth = exc_vars[2]



    w = exc_vars2[0]

    ab = exc_vars2[1]

    b = exc_vars2[2]

    s = exc_vars2[3]

    eb = exc_vars2[4]



    r_scoop = 3.125/5.3125*R        # Length of side



    #bucket_excavation_arc = ma.acos((R - depth)/R)

    # bucket_drum_volume = (w*ma.pi*(R**2)*(bucket_excavation_arc/360) -

    #                1/2*w*(R**2)*ma.sin(bucket_excavation_arc)*

    #                cos(bucket_excavation_arc))

    # scoop_volume = (w*ma.pi*(r_scoop**2)*(62/360) -

    #                1/2*w*(r_scoop**2)*ma.sin(62*ma.pi/180)*ma.cos(62*ma.pi/180))

    #shovel_excavation_arc = ma.acos((r - depth)/r)

    # bucket_shovel_volume = (w*ma.pi*(r**2)*(shovel_excavation_arc) -

    #                1/2*w*(r**2)*ma.sin(shovel_excavation_arc)*

    #                ma.cos(shovel_excavation_arc))

    #prism_volume = bucket_shovel_volume - scoop_volume - bucket_drum_volume;

    #q = prism_volume*gam



    D = depth    # Length of side



    # Balovnev function calcs

    if b < 0.5 * (ma.asin(ma.sin(d)/ma.sin(f)) - d):

        A1 = A(b, f)

    else:

        A1 = A_prime(b, f, d)



    if ab < 0.5 * (ma.asin(ma.sin(d)/ma.sin(f)) - d):

        A2 = A(ab, f)

    else:

        A2 = A_prime(ab, f, d)



    if ma.pi/2 < 0.5 * (ma.asin(ma.sin(d)/ma.sin(f)) - d):

        A3 = A(ma.pi/2, f)

    else:

        A3 = A_prime(ma.pi/2, f, d)



    l = r_scoop  # Length of scoop sides

    ls = r_scoop    # Length of scoop sides



    # Balovnev equation

    dig_force = (w*D*A1*(1 + 1/ma.tan(b)*ma.tan(d))*(D*4*g*gam/2 + c*1/ma.tan(f) +

                                                     (D - l*ma.sin(b))*g*gam*(1-ma.sin(f))/(1+ma.sin(f))) +

                 w*eb*A2*(1 + ma.tan(d)*1/ma.tan(ab))*(eb*g*gam/2 + c*1/ma.tan(f) + D*g*gam *

                                                       (1-ma.sin(f))/(1+ma.sin(f))) +

                 D*A3*(2*s + 4*ls*ma.tan(d))*(D*g*gam/2 + c*1/ma.tan(f) +

                                              (D - ls*ma.sin(b))*g*gam*(1-ma.sin(f))/(1+ma.sin(f))))

    return dig_force



# Balovnev functions





def A(x, f):

    ans = ((1 - ma.sin(f)*ma.cos(2*x))/(1 - ma.sin(f)))

    return ans





def A_prime(x, f, d):

    ans = ((ma.cos(d)*(ma.cos(d) + ((ma.sin(f))**2 - (ma.sin(d))**2)**0.5)/(1 - ma.sin(f))) *

           ma.exp((2*x - ma.pi + d + ma.asin(ma.sin(d)/ma.sin(f)))*ma.tan(f)))

    return ans





digOutputs = excavationMechanics(

    depthM, trenchDepthM, radiusM, extAngle, intAngle, cohCoeff, gVal, mRegolith)

AlphaExc = 2.77778e-7*digOutputs[1]  # kWh/kg





# Part added by Baptiste for adding the "transportation energy requirement while excavating"



velocityExcavation = 0.01  # [m/s] TO CROSSCHECK (the velocity during excavation)



# Empty rover:

mRegolith = 0  # [kg]

MinSlipEmpty = slip_required(mRover + mRegolith, gVal, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal, nVal, kcVal, kphiVal, cVal, kVal, phiVal)  # The minimum slip when the rover is empty

EnergyReqTptPerDistEmpty = energy_requirements(MinSlipEmpty, mRover + mRegolith, velocityExcavation, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal, gVal, nVal, kcVal, kphiVal, cVal, kVal, phiVal, motor_efficiencyVal)  # [J/m] Energy requirements for transport during excavation (empty rover)

EnergyReqTptEmpty = EnergyReqTptPerDistEmpty * digDistance  # [J]

# Full rover:

mRegolith = 90  # [kg]

MinSlipFull = slip_required(mRover + mRegolith, gVal, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal, nVal, kcVal, kphiVal, cVal, kVal, phiVal)  # The minimum slip when the rover is full

EnergyReqTptPerDistFull = energy_requirements(MinSlipFull, mRover + mRegolith, velocityExcavation, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal, gVal, nVal, kcVal, kphiVal, cVal, kVal, phiVal, motor_efficiencyVal)  # [J/m] Energy requirements for transport during excavation (full rover)

EnergyReqTptFull = EnergyReqTptPerDistFull * digDistance  # [J]

# Mean of empty and full rover:

EnergyReqTpt = (EnergyReqTptEmpty+EnergyReqTptFull)/2  # [J]

EnergyReqTpt = EnergyReqTpt * 2.77778e-7  # [kWh]

AlphaTpt = EnergyReqTpt/mRegolith  # [kWh/kg]



# Alpha taking into account both the excavation and the transportation during excavation

Alpha = AlphaExc + AlphaTpt  # [kWh/kg]



# Just to have an idea of their relative importance:

print("AlphaExc:", AlphaExc, "[kWh/kg]")

print("AlphaTpt:", AlphaTpt, "[kWh/kg]")

print("Alpha:", Alpha, "[kWh/kg]")

