"""
Author: Baptiste Valentin
"""

"""
Transportation Energy Calculation Module
----------------------------------------

This module calculates the energy required to transport regolith between two given sites.

**Description**
This module relies on a semi-empirical method of Terramechanics to analyze the locomotion properties of a vehicle on the lunar surface.
The global energy requirement calculated takes a round trip into consideration: an outward trip with an empty rover, and a return
trip with a rover full of regolith.

**Main variables**
- Soil parameters:
    ° 'gVal': gravitational constant [m/s2]
    ° 'rhoVal': density [kg/m3]
    ° 'nVal': Exponent describing soil compressibility []
    ° 'kcVal': Cohesion modulus [N/m^n+1]
    ° 'kphiVal': Friction modulus [N/m^n+2]
    ° 'cVal': Cohesion [Pa]
    ° 'phiVal': Friction angle [rad]
    ° 'kVal': Shear modulus [m]
- Vehicle parameters:
    ° 'velocityVal': velocity of the vehicle [m/s]
    ° 'motor_efficiencyVal': motor efficiency []
    ° 'mRover': mass of the vehicle [kg]
    ° 'mRegolith": mass of the conveyed regolith (per trip) [kg]
    ° 'wheelbaseVal': distance between the front wheels' axle and the rear wheels' axle [m]
    ° 'heightCOGVal': height of the center of gravity [m]
- Wheel parameters:
    ° 'WheelRadiusVal': radius of the wheel(s) [m]
    ° 'WheelWidthVal': width of the wheel(s) [m] 
- Trip parameters:
    ° 'SlopeVal': Slope angle of the soil where the vehicle is traveling [rad]
    ° 'DistanceToTravel': distance between the excavation and the beneficiation site [m]
"""

########################################### imports #############################################################


from scipy.optimize import root_scalar
import scipy.integrate as integrate
import numpy as np
import math
import warnings


#################################### Definition of useful functions ##############################################


def Fz_front(mTotVar, gVar, slopeVar, lVar, hVar):
    """
    Calculates the vertical reaction force on one of the front wheels.

    Args:
        mTotVar (float): total mass of the vehicle [kg]
        gVar (float): gravitational constant [m/s2]
        slopeVar (float): slope angle of the soil where the vehicle is traveling [rad]
        lVar (float): distance between the front wheels' axle and the rear wheels' axle [m]
        hVar (float): height of the center of gravity [m]

    Returns:
        float: vertical reaction force on one of the front wheels [N].
    """

    if slopeVar >= 0:

        return (mTotVar / 2 * gVar * math.cos(slopeVar + math.atan(hVar / (lVar / 2))) * (math.sqrt((lVar / 2) ** 2 + hVar ** 2))) / (lVar)

        # return m_tot / 4 * g (for check purposes, should give the same when the slope is zero).

    else:

        return (mTotVar / 2 * gVar * math.cos(math.pi / 2 - slopeVar - math.atan((lVar / 2) / hVar)) * (math.sqrt((lVar / 2) ** 2 + hVar ** 2))) / (lVar)


def Fz_rear(mTotVar, gVar, slopeVar, lVar, hVar):
    """
    Calculates the vertical reaction force on one of the rear wheels.

    Args:
        mTotVar (float): total mass of the vehicle [kg]
        gVar (float): gravitational constant [m/s2]
        slopeVar (float): slope angle of the soil where the vehicle is traveling [rad]
        lVar (float): distance between the front wheels' axle and the rear wheels' axle [m]
        hVar (float): height of the center of gravity [m]

    Returns:
        float: vertical reaction force on one of the rear wheels [N].
    """
    if slopeVar >= 0:

        return (mTotVar / 2 * gVar * math.cos(math.pi / 2 - slopeVar - math.atan((lVar / 2) / hVar)) * (math.sqrt((lVar / 2) ** 2 + hVar ** 2))) / (lVar)

        # return m_tot / 4 * g (for check purposes, should give the same when the slope is zero).

    else:

        return (mTotVar / 2 * gVar * math.cos(slopeVar + math.atan(hVar / (lVar / 2))) * (math.sqrt((lVar / 2) ** 2 + hVar ** 2))) / (lVar)


def sigma(thetaVar, theta0Var, bVar, rVar, nVar, kcVar, kphiVar):
    """
    Calculates the normal/radial component of the contact stress between the wheel and the soil.

    Args:
        thetaVar: angular position along the wheel-terrain interface [rad]
        thetaOVar: entry angle (where the wheel first contacts the terrain [rad]
        bVar: width of the wheel(s) [m]
        rVar: radius of the wheel(s) [m]
        nVar: Exponent describing soil compressibility []
        kcVar: Cohesion modulus [N/m^n+1]
        kphiVar: Friction modulus [N/m^n+2]

    Returns:
        float: the normal/radial component of the contact stress between the wheel and the soil at a given angular position [Pa].
    """

    return (rVar ** nVar) * (kcVar / bVar + kphiVar) * (math.cos(thetaVar) - math.cos(theta0Var)) ** nVar


def delta_s(thetaVar, theta0Var, sVar, rVar):
    """
    Calculates the shear displacement (how much the soil deforms under the wheel due to shear forces).

    Args:
        thetaVar: angular position along the wheel-terrain interface [rad]
        thetaOVar: entry angle (where the wheel first contacts the terrain [rad]
        svar: slip ratio []
        rVar: radius of the wheel(s) [m]

    Returns:
        float: the shear displacement at a given angular position [m].
    """

    return rVar * ((theta0Var - thetaVar) - (1 - sVar) * (math.sin(theta0Var) - math.sin(thetaVar)))


def tau(thetaVar, theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):
    """
    Calculates the tangential component of the contact stress between the wheel and the soil.

    Args:
        thetaVar: angular position along the wheel-terrain interface [rad]
        thetaOVar: entry angle (where the wheel first contacts the terrain [rad]
        svar: slip ratio []
        bVar: width of the wheel(s) [m]
        rVar: radius of the wheel(s) [m]
        nVar: Exponent describing soil compressibility []
        kcVar: Cohesion modulus [N/m^n+1]
        kphiVar: Friction modulus [N/m^n+2]
        cVar: Cohesion [Pa]
        kVar: Shear modulus [m]
        phiVar: Friction angle [rad]

    Returns:
        float: the tangential component of the contact stress between the wheel and the soil at a given angular position [Pa].
    """

    return (cVar + sigma(thetaVar, theta0Var, bVar, rVar, nVar, kcVar, kphiVar) * math.tan(phiVar)) * (1 - math.exp(-delta_s(thetaVar, theta0Var, sVar, rVar) / kVar))


def fun_to_integrate_sigma(thetaVar, theta0Var, bVar, rVar, nVar, kcVar, kphiVar):

    return sigma(thetaVar, theta0Var, bVar, rVar, nVar, kcVar, kphiVar) * math.cos(thetaVar)


def fun_to_integrate_tau(thetaVar, theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):

    return tau(thetaVar, theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar) * math.sin(thetaVar)


def fun_to_integrate_traction(thetaVar, theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):

    return tau(thetaVar, theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar) * math.cos(thetaVar)


def fun_to_solve_front(theta0Var, sVar, mTotVar, gVar, slopeVar, lVar, hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):

    return Fz_front(mTotVar, gVar, slopeVar, lVar, hVar) - bVar * rVar * integrate.quad(fun_to_integrate_sigma, 0, theta0Var, args=(theta0Var, bVar, rVar, nVar, kcVar, kphiVar))[0] - bVar * rVar * integrate.quad(fun_to_integrate_tau, 0, theta0Var, args=(theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]


def fun_to_solve_rear(theta0Var, sVar, mTotVar, gVar, slopeVar, lVar, hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):

    return Fz_rear(mTotVar, gVar, slopeVar, lVar, hVar) - bVar * rVar * integrate.quad(fun_to_integrate_sigma, 0, theta0Var, args=(theta0Var, bVar, rVar, nVar, kcVar, kphiVar))[0] - bVar * rVar * integrate.quad(fun_to_integrate_tau, 0, theta0Var, args=(theta0Var, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]


def fNy(phiVar):
    """Function that extrapolates the value of Ny (soil weight bearing capacity factor) based on phiVar (friction angle)"""
    phiVec = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,

              29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

    NyVec = [0.0, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.0, 1.2, 1.4, 1.6, 1.9, 2.2, 2.5, 2.9, 3.3, 3.8, 4.4,

             5.1, 5.9, 6.8, 7.9, 9.2, 10.7, 12.5, 14.6, 17.1, 20.1, 23.7, 28.0, 33.3, 39.6, 47.3, 56.7, 68.1, 82.3,

             99.8]

    zNy = np.polyfit(phiVec, NyVec, 4)

    fNy = np.poly1d(zNy)

    NyVar = fNy(phiVar)

    return NyVar


def fNc(phiVar):
    """Function that extrapolates the value of Nc (cohesion bearing capacity factor) based on phiVar (friction angle)"""
    phiVec = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,

              29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

    NcVec = [5.7, 6.0, 6.3, 6.6, 7.0, 7.3, 7.7, 8.2, 8.6, 9.1, 9.6, 10.2, 10.8, 11.4, 12.1, 12.6, 13.7, 14.6, 15.5,

             16.6, 17.7, 18.9, 20.3, 21.7, 23.4, 25.1, 27.1, 29.2, 31.6, 34.2, 37.2, 40.4, 44.0, 48.1, 52.6, 57.8, 63.5,

             70.1, 77.5, 86]

    zNc = np.polyfit(phiVec, NcVec, 4)

    fNc = np.poly1d(zNc)

    NcVar = fNc(phiVar)

    return NcVar


def drawbar_pull(sVar, mTotVar, gVar, bVar, rVar, slopeVar, lVar, hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):
    """
    Calculates the drawbar pull of the rover (net tractive force available for moving the vehicle).

    Args:
        sVar: slip ratio []
        mTotVar: total mass of the vehicle [kg]
        gVar: gravitational constant [m/s2]
        bVar: width of the wheel(s) [m]
        rVar: radius of the wheel(s) [m]
        slopeVar: slope angle of the soil where the vehicle is traveling [rad]
        lVar: distance between the front wheels' axle and the rear wheels' axle [m]
        hVar: height of the center of gravity [m]
        nVar: Exponent describing soil compressibility []
        kcVar: Cohesion modulus [N/m^n+1]
        kphiVar: Friction modulus [N/m^n+2]
        cVar: Cohesion [Pa]
        kVar: Shear modulus [m]
        phiVar: Friction angle [rad]

    Returns:
        float: the drawbar pull [N].
    """

    sol_theta0_front = root_scalar(fun_to_solve_front, args=(sVar, mTotVar, gVar, slopeVar, lVar,
                                   hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

    z0_front = rVar * (1 - math.cos(sol_theta0_front.root)
                       )  # Sinkage of the front wheels [m]

    l0_front = sol_theta0_front.root * rVar

    # Verify it makes sense

    if l0_front < bVar:

        newb_front = min(l0_front, bVar)

        sol_theta0_front = root_scalar(fun_to_solve_front, args=(sVar, mTotVar, gVar, slopeVar, lVar, hVar,
                                       newb_front, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

        z0_front = rVar * (1 - math.cos(sol_theta0_front.root))

    # End verify

    sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(sVar, mTotVar, gVar, slopeVar, lVar,
                                  hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

    # Sinkage of the rear wheels [m]
    z0_rear = rVar * (1 - math.cos(sol_theta0_rear.root))

    l0_rear = sol_theta0_rear.root * rVar

    # Verify it makes sense

    if l0_rear < bVar:

        newb_rear = min(l0_rear, bVar)

        sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(sVar, mTotVar, gVar, slopeVar, lVar, hVar,
                                      newb_rear, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

        z0_rear = rVar * (1 - math.cos(sol_theta0_rear.root))

    # End verify

    Ft_front = bVar * rVar * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_front.root, args=(
        # Traction effort of a front wheel [N]
        sol_theta0_front.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    Ft_rear = bVar * rVar * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_rear.root, args=(
        # Traction effort of a rear wheel [N]
        sol_theta0_rear.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    Rc_front = bVar * (kcVar / bVar + kphiVar) * \
        (z0_front ** (nVar + 1)) / \
        (nVar + 1)  # Compaction resistance of a front wheel [N]

    Rc_rear = bVar * (kcVar / bVar + kphiVar) * \
        (z0_rear ** (nVar + 1)) / \
        (nVar + 1)  # Compaction resistance of a rear wheel [N]

    R_slope = mTotVar * gVar * math.sin(slopeVar)  # Slope resistance [N]

    db_pull = 2 * Ft_front + 2 * Ft_rear - 2 * Rc_front - 2 * \
        Rc_rear - R_slope  # Drawbar pull (of the whole vehicle) [N]

    return db_pull


def energy_requirements(sVar, mTotVar, vVar, bVar, rVar, slopeVar, lVar, hVar, gVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar, etaVar):
    """Function that computes the energy requirements for a rover with a given slip, mass, velocity, wheel width, wheel radius, and slope."""

    sol_theta0_front = root_scalar(fun_to_solve_front, args=(sVar, mTotVar, gVar, slopeVar, lVar,
                                   hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

    z0_front = rVar * (1 - math.cos(sol_theta0_front.root))

    l0_front = sol_theta0_front.root * rVar

    # Verify it makes sense

    if l0_front < bVar:

        newb_front = min(l0_front, bVar)

        sol_theta0_front = root_scalar(fun_to_solve_front, args=(sVar, mTotVar, gVar, slopeVar, lVar, hVar, newb_front, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar),

                                       method='toms748', bracket=[0, 1])

        z0_front = rVar * (1 - math.cos(sol_theta0_front.root))

    # End verify

    sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(sVar, mTotVar, gVar, slopeVar, lVar,
                                  hVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

    z0_rear = rVar * (1 - math.cos(sol_theta0_rear.root))

    l0_rear = sol_theta0_rear.root * rVar

    # Verify it makes sense

    if l0_rear < bVar:

        newb_rear = min(l0_rear, bVar)

        sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(sVar, mTotVar, gVar, slopeVar, lVar, hVar, newb_rear, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748',

                                      bracket=[0, 1])

        z0_rear = rVar * (1 - math.cos(sol_theta0_rear.root))

    # End verify

    Ft_front = bVar * rVar * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_front.root, args=(
        # Traction effort of a front wheel [N]
        sol_theta0_front.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    Ft_rear = bVar * rVar * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_rear.root, args=(
        # Traction effort of a rear wheel [N]
        sol_theta0_rear.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    Rc_front = bVar * (kcVar / bVar + kphiVar) * \
        (z0_front ** (nVar + 1)) / \
        (nVar + 1)  # Compaction resistance of a front wheel [N]

    Rc_rear = bVar * (kcVar / bVar + kphiVar) * \
        (z0_rear ** (nVar + 1)) / \
        (nVar + 1)  # Compaction resistance of a rear wheel [N]

    Tr_front = bVar * rVar ** 2 * integrate.quad(tau, 0, sol_theta0_front.root, args=(
        # Torque exerted on a front wheel [N.m]
        sol_theta0_front.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    Tr_rear = bVar * rVar ** 2 * integrate.quad(tau, 0, sol_theta0_rear.root, args=(
        # Torque exerted on a rear wheel [N.m]
        sol_theta0_rear.root, sVar, bVar, rVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar))[0]

    omega_r = vVar / (rVar * (1 - sVar))  # Rotation speed [rad/s]

    Power_front = Tr_front * omega_r  # Power required for a front wheel [W]

    Power_rear = Tr_rear * omega_r  # Power required for a rear wheel [W]

    Locomotion_power = (2 * Power_front + 2 * Power_rear) / (etaVar)  # [W]

    Locomotion_energy_km = Locomotion_power / vVar  # [J/m]

    return Locomotion_energy_km


def slip_required(mTotVar, gVar, bVar, rVar, slopeVar, lVar, hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):
    """
    Calculates the slip required so that the vehicle can move (by researching where the drawbar pull = 0).

    Args:
        mTotVar: total mass of the vehicle [kg]
        gVar: gravitational constant [m/s2]
        bVar: width of the wheel(s) [m]
        rVar: radius of the wheel(s) [m]
        slopeVar: slope angle of the soil where the vehicle is traveling [rad]
        lVar: distance between the front wheels' axle and the rear wheels' axle [m]
        hVar: height of the center of gravity [m]
        nVar: Exponent describing soil compressibility []
        kcVar: Cohesion modulus [N/m^n+1]
        kphiVar: Friction modulus [N/m^n+2]
        cVar: Cohesion [Pa]
        kVar: Shear modulus [m]
        phiVar: Friction angle [rad]

    Returns:
        float: the slip required so that the vehicle can move [].
    """

    # drawbar pull when slip = 0.
    Drawbar_pull_slip_0 = drawbar_pull(
        0, mTotVar, gVar, bVar, rVar, slopeVar, lVar, hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar)

    if Drawbar_pull_slip_0 > 0:

        Required_slip = 0

    else:

        Required_slip = root_scalar(drawbar_pull, args=(mTotVar, gVar, bVar, rVar, slopeVar, lVar,
                                    hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar), method='toms748', bracket=[0, 1])

        Required_slip = Required_slip.root

    if Required_slip > 0.5:

        Required_slip = np.nan

        warnings.warn("The required slip is too large: motion is impossible")

    return Required_slip


def compute_beta(mRoverVar, mRegolithVar, gVar, bVar, rVar, slopeVar, lVar, hVar, vVar, etaVar, distanceVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar):
    """
    Calculates the energy requirements for the transportation module.

    Returns:
        float: the energy requirements for transporting regolith between two given sites, per kilogram regolith and per kilometer [kJ/kg/km].
    """
    MassOutwardTrip = mRoverVar

    MassReturnTrip = mRoverVar + mRegolithVar

    # Outward trip (without regolith)

    RequiredSlipOutward = slip_required(
        MassOutwardTrip, gVar, bVar, rVar, slopeVar, lVar, hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar)

    EnergyPerDistanceOutward = energy_requirements(RequiredSlipOutward, MassOutwardTrip, vVar, bVar,
                                                   # [J/m]
                                                   rVar, slopeVar, lVar, hVar, gVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar, etaVar)

    EnergyOutward = EnergyPerDistanceOutward * distanceVar  # [J]

    # Return trip (with regolith)

    RequiredSlipReturn = slip_required(
        MassReturnTrip, gVar, bVar, rVar, slopeVar, lVar, hVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar)

    EnergyPerDistanceReturn = energy_requirements(RequiredSlipReturn, MassReturnTrip, vVar, bVar,
                                                  # [J/m]
                                                  rVar, slopeVar, lVar, hVar, gVar, nVar, kcVar, kphiVar, cVar, kVar, phiVar, etaVar)

    EnergyReturn = EnergyPerDistanceReturn * distanceVar  # [J]

    # Round trip

    EnergyRoundTrip = EnergyOutward + EnergyReturn

    Beta = EnergyRoundTrip/(mRegolithVar*distanceVar)  # [J/kg/m] or [kJ/kg/km]

    return Beta


################################### ASSUMPTIONS: regolith properties #############################################
# Regolith properties (from book "Introduction to the Mechanics of Space Robots").
gVal = 1.62  # 1.62  # Gravity (m/s2)

rhoVal = 1600  # Density (kg/m3)

nVal = 1  # Exponent

kcVal = 1100  # Cohesion modulus (N/m^n+1) '#1400

kphiVal = 820000  # Friction modulus (N/m^n+2)

cVal = 170  # Cohesion (Pa)

phiVal = math.radians(45)  # Friction angle (rad)

kVal = 18e-3  # Shear modulus (m)


############################################ Other ASSUMPTIONS ###################################################

# [m/s], maximum speed from "RASSOR, the reduced gravity excavator."
velocityVal = 0.49

motor_efficiencyVal = 0.6

mRover = 67  # [kg] Total mass of RASSOR2

# [kg] The rover is assumed to transport the maximum amount of regolith each time (from "RASSOR, the reduced gravity excavator.").
mRegolith = 90

WheelRadiusVal = 0.4318 / 2  # [m] Wheel radius (17 inches for the full wheel)

WheelWidthVal = 0.1  # [m] Wheel width

wheelbaseVal = 0.5  # [m]

heightCOGVal = 0.1  # [m] (The CoG is assumed to be centered)

SlopeVal = 0  # [rad] The soil is assumed to be flat

# [m] Distance between the excavation and the beneficiation site
DistanceToTravel = 1000

NcVal = fNy(phiVal*180/math.pi)  # Coefficient based on phi

NyVal = fNc(phiVal*180/math.pi)  # Coefficient based on phi


################################################ Compute Beta ##################################################


Beta = compute_beta(mRover, mRegolith, gVal, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal,
                    # [kJ/kg/km]
                    velocityVal, motor_efficiencyVal, DistanceToTravel, nVal, kcVal, kphiVal, cVal, kVal, phiVal)


def get_Beta(motor_efficiency, mRover_input):
    beta_return = compute_beta(mRover_input, mRegolith, gVal, WheelWidthVal, WheelRadiusVal, SlopeVal, wheelbaseVal, heightCOGVal,
                               velocityVal, motor_efficiency, DistanceToTravel, nVal, kcVal, kphiVal, cVal, kVal, phiVal)
    return beta_return/3600


Beta = Beta/3600  # [kWh/kg/km]


# print("Beta:", Beta, "[kWh/kg/km]")
