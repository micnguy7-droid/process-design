# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 18:24:21 2022

@author: fardi
"""
# imports
import math
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import warnings

# ASSUMPTIONS: regolith properties
# Regolith properties (from book "Introduction to the Mechanics of Space Robots").
g = 1.62  #1.62  # Gravity (m/s2)
rho = 1600  #1600  # Density (kg/m3), 1600 is taken as reference.
n = 1  # Exponent
kc = 2100  # Cohesion modulus (N/m^n+1) '#2100
kphi = 820000  # Friction modulus (N/m^n+2)
c = 170  # Cohesion (Pa)
phi = math.radians(45)  # Friction angle (rad) 'CHANGED IT TO 45. 37 before
k = 18e-3  # Shear modulus (m)
Nc = 70.1  # Coefficient based on phi
Ny = 68.1  # Coefficient based on phi

# Forces
def Fz_front(m_tot, slope, l, h): # l is the distance between the front and rear wheels, h is the height of the center of gravity
    if slope >= 0:
        return (m_tot / 2 * g * math.cos(slope + math.atan(h / (l / 2))) * (math.sqrt((l / 2) ** 2 + h **2))) / (l) # Vertical force on the front wheels
        #return m_tot / 4 * g (for check purposes, should give the same when the slope is zero).
    else:
        return (m_tot / 2 * g * math.cos(math.pi / 2 - slope - math.atan((l / 2 )/ h)) * (math.sqrt((l / 2) ** 2 + h **2))) / (l)


def Fz_rear(m_tot, slope, l, h):
    if slope >= 0:
        return (m_tot / 2 * g * math.cos(math.pi / 2 - slope - math.atan((l / 2 )/ h)) * (math.sqrt((l / 2) ** 2 + h **2))) / (l) # Vertical force on the rear wheels
        #return m_tot / 4 * g
    else:
        return (m_tot / 2 * g * math.cos(slope + math.atan(h / (l / 2))) * (math.sqrt((l / 2) ** 2 + h **2))) / (l)


# Tractive effort
def sigma(theta, theta0, b, r):
    return (r ** n) * (kc / b + kphi) * (math.cos(theta) - math.cos(theta0)) ** n


def delta_s(theta, theta0, s, r):
    return r * ((theta0 - theta) - (1 - s) * (math.sin(theta0) - math.sin(theta)))


def tau(theta, theta0, s, b, r):
    return (c + sigma(theta, theta0, b, r) * math.tan(phi)) * (1 - math.exp(-delta_s(theta, theta0, s, r) / k))


def fun_to_integrate_sigma(theta, theta0, b, r):
    return sigma(theta, theta0, b, r) * math.cos(theta)


def fun_to_integrate_tau(theta, theta0, s, b, r):
    return tau(theta, theta0, s, b, r) * math.sin(theta)


def fun_to_integrate_traction(theta, theta0, s, b, r):
    return tau(theta, theta0, s, b, r) * math.cos(theta)


def fun_to_solve_front(theta0, s, m_tot, slope, l, h, b, r):
    return Fz_front(m_tot, slope, l, h) - b * r * integrate.quad(fun_to_integrate_sigma, 0, theta0, args=(theta0, b, r))[0] - b * r * integrate.quad(fun_to_integrate_tau, 0, theta0, args=(theta0, s, b, r))[0]


def fun_to_solve_rear(theta0, s, m_tot, slope, l, h, b, r):
    return Fz_rear(m_tot, slope, l, h) - b * r * integrate.quad(fun_to_integrate_sigma, 0, theta0, args=(theta0, b, r))[0] - b * r * integrate.quad(fun_to_integrate_tau, 0, theta0, args=(theta0, s, b, r))[0]


# Function that computes the drawbar pull of the rover from its slip, mass, wheel width, wheel radius, and the slope value.
def drawbar_pull(slip, mass, b, r, slope, l, h):
    sol_theta0_front = root_scalar(fun_to_solve_front, args=(slip, mass, slope, l, h, b, r), method='toms748', bracket=[0, 1])
    z0_front = r * (1 - math.cos(sol_theta0_front.root))
    l0_front = sol_theta0_front.root * r
    # Verify it makes sense
    if l0_front < b:
        newb_front = min(l0_front, b)
        sol_theta0_front = root_scalar(fun_to_solve_front, args=(slip, mass, slope, l, h, newb_front, r), method='toms748', bracket=[0, 1])
        z0_front = r * (1 - math.cos(sol_theta0_front.root))
    # End verify
    #print("Theta0_front =", sol_theta0_front.root * 180 / math.pi, "[°] and z0_front =", z0_front * 10 ** (3), "[mm]")
    sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(slip, mass, slope, l, h, b, r), method='toms748', bracket=[0, 1])
    z0_rear = r * (1 - math.cos(sol_theta0_rear.root))
    l0_rear = sol_theta0_rear.root * r
    # Verify it makes sense
    if l0_rear < b:
        newb_rear = min(l0_rear, b)
        sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(slip, mass, slope, l, h, newb_rear, r), method='toms748', bracket=[0, 1])
        z0_rear = r * (1 - math.cos(sol_theta0_rear.root))
    # End verify
    #print("Theta0_rear =", sol_theta0_rear.root * 180 / math.pi, "[°] and z0_rear =", z0_rear * 10 ** (3), "[mm]")
    # Traction effort (MF010: page 306)
    Ft_front = b * r * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_front.root, args=(sol_theta0_front.root, slip, b, r))[0]
    Ft_rear = b * r * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_rear.root, args=(sol_theta0_rear.root, slip, b, r))[0]
    #print("Tractive effort on front/rear wheels:", Ft_front, "[N] /", Ft_rear, "[N]")

    # Compaction resistance with these sinkages
    Rc_front = b * (kc / b + kphi) * (z0_front ** (n + 1)) / (n + 1)  # [N]
    Rc_rear = b * (kc / b + kphi) * (z0_rear ** (n + 1)) / (n + 1)  # [N]
    #print("Compaction resistance on front/rear wheels:", Rc_front, "[N] /", Rc_rear, "[N]")

    # Slope resistance
    R_slope = m_tot * g * math.sin(slope)
    # print("Slope resistance (total):", R_slope, "[N]")

    # Drawbar pull (of the whole vehicle)
    db_pull = 2 * Ft_front + 2 * Ft_rear - 2 * Rc_front - 2 * Rc_rear - R_slope
    #print("Drawbar pull:", db_pull, "[N]")

    return db_pull


# Function that computes the energy requirements for a rover with a given slip, mass, velocity, wheel width, wheel radius, and slope.
def energy_requirements(slip, mass, velocity, b, r, slope, l, h):
    sol_theta0_front = root_scalar(fun_to_solve_front, args=(slip, mass, slope, l, h, b, r), method='toms748', bracket=[0, 1])
    z0_front = r * (1 - math.cos(sol_theta0_front.root))
    l0_front = sol_theta0_front.root * r
    # Verify it makes sense
    if l0_front < b:
        newb_front = min(l0_front, b)
        sol_theta0_front = root_scalar(fun_to_solve_front, args=(slip, mass, slope, l, h, newb_front, r),
                                       method='toms748', bracket=[0, 1])
        z0_front = r * (1 - math.cos(sol_theta0_front.root))
    # End verify
    # print("Theta0_front =", sol_theta0_front.root * 180 / math.pi, "[°] and z0_front =", z0_front * 10 ** (3), "[mm]")
    sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(slip, mass, slope, l, h, b, r), method='toms748', bracket=[0, 1])
    z0_rear = r * (1 - math.cos(sol_theta0_rear.root))
    l0_rear = sol_theta0_rear.root * r
    # Verify it makes sense
    if l0_rear < b:
        newb_rear = min(l0_rear, b)
        sol_theta0_rear = root_scalar(fun_to_solve_rear, args=(slip, mass, slope, l, h, newb_rear, r), method='toms748',
                                      bracket=[0, 1])
        z0_rear = r * (1 - math.cos(sol_theta0_rear.root))
    # End verify
    # print("Theta0_rear =", sol_theta0_rear.root * 180 / math.pi, "[°] and z0_rear =", z0_rear * 10 ** (3), "[mm]")
    # Traction effort (MF010: page 306)
    Ft_front = b * r * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_front.root, args=(sol_theta0_front.root, slip, b, r))[0]
    Ft_rear = b * r * integrate.quad(fun_to_integrate_traction, 0, sol_theta0_rear.root, args=(sol_theta0_rear.root, slip, b, r))[0]
    print("Tractive effort on front/rear wheels:", Ft_front, "[N] /", Ft_rear, "[N]")

    # Compaction resistance with these sinkages
    Rc_front = b * (kc / b + kphi) * (z0_front ** (n + 1)) / (n + 1)  # [N]
    Rc_rear = b * (kc / b + kphi) * (z0_rear ** (n + 1)) / (n + 1)  # [N]
    # print("Compaction resistance on front/rear wheels:", Rc_front, "[N] /", Rc_rear, "[N]")

    # Torque exerted on the wheel (MF010: page 306)
    # Solution: new computation of tau with the minimum of b and l0_rear/front
    Tr_front = b * r ** 2 * integrate.quad(tau, 0, sol_theta0_front.root, args=(sol_theta0_front.root, slip, b, r))[0]
    Tr_rear = b * r ** 2 * integrate.quad(tau, 0, sol_theta0_rear.root, args=(sol_theta0_rear.root, slip, b, r))[0]
    # print("Torque exerted on the front/rear wheels:", Tr_front, "[N] /", Tr_rear, "[N]")

    # Rotation speed
    omega_r = velocity / (r * (1 - slip))
    # print("Rotation speed wheels:", omega_r, "[rad/s]")

    # Power requirements (we will compute it only for the right slip = 0)
    Power_front = Tr_front * omega_r
    Power_rear = Tr_rear * omega_r
    Locomotion_power = (2 * Power_front + 2 * Power_rear) / (motor_efficiency)  # [W]
    Locomotion_energy_km = Locomotion_power / velocity  # [J/m]

    #    print("Total power requirement for locomotion:", Locomotion_power, "[W]")
    #    print("Energy/km:", Locomotion_energy_km, "[J/m]")

    return Locomotion_energy_km

# Function that computes the slip required so that the rover can move
def slip_required(mass, b, r, slope, l, h):
    # Research where the drawbar pull = 0
    Drawbar_pull_slip_0 = drawbar_pull(0, mass, b, r, slope, l, h) # drawbar pull when slip = 0.
    if Drawbar_pull_slip_0 > 0:
        Required_slip = 0
    else:
        Required_slip = root_scalar(drawbar_pull, args=(mass, b, r, slope, l, h), method='toms748', bracket=[0, 1])
        Required_slip = Required_slip.root
    if Required_slip > 0.5:
        Required_slip = np.nan
        warnings.warn("The required slip is too large: motion is impossible")
    return Required_slip


# ASSUMPTIONS: velocity, motor efficiency, total mass, wheel dimensions, wheelbase, height of CoG.
# Parameters (90 kg of regolith per trip)
velocity = 0.49  # [m/s], maximum speed from "RASSOR, the reduced gravity excavator."
motor_efficiency = 0.6  #0.6  #
m_tot = 67 + 90  # [kg] Total mass of RASSOR2 + mass of regolith carried (from "RASSOR, the reduced gravity excavator.")
WheelRadius = 0.4318 / 2  # [m] Wheel radius (17 inches for the full wheel)
WheelWidth = 0.1  # [m] Wheel width (from the picture?)
wheelbase = 0.5  # [m]
heightCOG = 0.1  # [m]
# The CoG is assumed to be centered

# Test and compare with TerraNOVA
#print("Wheel loading front: ", Fz_front(m_tot, 0, wheelbase, heightCOG))
#print("Wheel loading rear: ", Fz_rear(m_tot, 0, wheelbase, heightCOG))
#test = drawbar_pull(0.2, m_tot, WheelWidth, WheelRadius, 0, wheelbase, heightCOG)
#print("Total motor power:", velocity*energy_requirements(0.2, m_tot, velocity, WheelWidth, WheelRadius, 0, wheelbase, heightCOG), "[W]")

# What values do we want to compute?
continueBoolean = 0  # to get the energy requirement on a flat surface [J/m], with full regolith load onboard.
continueBoolean1 = 0  # to get a graph of energy requirements as a function of slope.
continueBoolean2 = 0  # to get a graph of energy requirement as a function of transported regolith mass.
continueBoolean3 = 1  # to get the beta values we want! [kWh/kg/km]
continueBoolean4 = 0  # to get an estimate of the influence of depth (to be continued)

if continueBoolean == 1:
    # Plot with drawbar pull as a function of slip on flat soil
    mPlatform = 67  # [kg]
    mPayload = 90  # [kg]
    mPayload = np.array(mPayload)
    mTot = mPlatform + mPayload
    # ASSUMPTION: slope value
    slopeval = math.radians(0) # [rad]
    slipval = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    db_vector = np.zeros(len(slipval))

    i = 0
    while i < len(slipval):
        db_vector[i] = drawbar_pull(slipval[i], mTot, WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
        i = i + 1

    # Research where the drawbar pull = 0
    Required_slip = slip_required(mTot, WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
    print("Required slip:", Required_slip, "[]")

    plot_mass = plt.figure(1)
    plt.plot(slipval, db_vector)
    plt.plot([0, 0.1], [0, 0], linestyle="--")
    plt.plot([Required_slip, Required_slip], [-20, 0], linestyle="--")
    plt.ylabel('Drawbar pull [N]')
    plt.xlabel('Slip []')
    plt.show()  # it has the expected shape!

    # Research where the drawbar pull = 0
    Required_slip = slip_required(mTot, WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
    print("Required slip:", Required_slip, "[]")

    # Reuse the function to compute power required at this slip
    EnergyPerKm = energy_requirements(Required_slip, mTot, velocity, WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
    print("Energy/km:", EnergyPerKm, "[J/m]")
    PowerConsumption = EnergyPerKm * velocity
    print("Power consumption:", PowerConsumption, "[W]")

if continueBoolean1 == 1:
    #  Full avec SLOPE qui varie
    mPlatform = 67  # [kg]
    mPayload = 90  # [kg]
    mPayload = np.array(mPayload)
    mTot = mPlatform + mPayload
    slopeval = np.arange(-10, 20, 1) * math.pi/180  # [rad]

    Required_slip_vec_slope = np.zeros(len(slopeval))
    Loc_energy_km_vec_slope = np.zeros(len(slopeval))

    i = 0

    while i < len(slopeval):
        Required_slip_vec_slope[i] = slip_required(mTot, WheelWidth, WheelRadius, slopeval[i], wheelbase, heightCOG)
        if not np.isnan(Required_slip_vec_slope[i]):
            Loc_energy_km_vec_slope[i] = energy_requirements(Required_slip_vec_slope[i], mTot, velocity, WheelWidth, WheelRadius, slopeval[i], wheelbase, heightCOG)
        i = i + 1

    print(Required_slip_vec_slope)
    plot_mass = plt.figure(1)
    plt.bar(slopeval*180/math.pi, Loc_energy_km_vec_slope)
    plt.ylabel('Locomotion energy requirements [J/m]')
    plt.xlabel('Climbing slope [°]')
    plt.show()  # it has the expected shape!


if continueBoolean2 == 1:
    #  Full avec mass regolith qui varie
    mPlatform = 67  # [kg]
    mass_regolith = np.arange(0, 90, 1)  # [kg]
    slopeval = 0

    Required_slip_vec_mass = np.zeros(len(mass_regolith))
    Loc_energy_km_vec_mass = np.zeros(len(mass_regolith))

    i = 0

    while i < len(mass_regolith):
        Required_slip_vec_mass[i] = slip_required(mPlatform+mass_regolith[i], WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
        if not np.isnan(Required_slip_vec_mass[i]):
            Loc_energy_km_vec_mass[i] = energy_requirements(Required_slip_vec_mass[i], mPlatform+mass_regolith[i], velocity, WheelWidth, WheelRadius, slopeval, wheelbase, heightCOG)
        i = i + 1
        # Why does the required energy per kilometer decreases when slope increases???

    print(Required_slip_vec_mass)
    plot_mass = plt.figure(1)
    plt.bar(mass_regolith, Loc_energy_km_vec_mass)
    plt.ylabel('Locomotion energy requirements [J/m]')
    plt.xlabel('Regolith mass [kg]')
    plt.show()  # it has the expected shape!


if continueBoolean3 == 1:
    velocity = 0.49  # [m/s], from "RASSOR, the reduced gravity excavator."
    motor_efficiency = 0.6  #
    mRover = 67  # [kg] Total mass of RASSOR2
    mRegolith = 90  # [kg] The rover is assumed to transport the maximum amount of regolith each time (from "RASSOR, the reduced gravity excavator.").
    WheelRadius = 0.4318 / 2  # [m] Wheel radius (17 inches for the full wheel)
    WheelWidth = 0.1  # [m] Wheel width (from the picture?)
    wheelbase = 0.5  # [m]
    heightCOG = 0.1  # [m]
    Slope = 0  # [rad] The soil is assumed to be flat
    # The CoG is assumed to be centered

    DistanceToTravel = 1000  # [m]
    MassOutwardTrip = mRover
    MassReturnTrip = mRover + mRegolith

    # Outward trip (without regolith)
    RequiredSlipOutward = slip_required(MassOutwardTrip, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)
    EnergyPerDistanceOutward = energy_requirements(RequiredSlipOutward, MassOutwardTrip, velocity, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)  # [J/m]
    EnergyOutward = EnergyPerDistanceOutward * DistanceToTravel  # [J]

    # Return trip (with regolith)
    RequiredSlipReturn = slip_required(MassReturnTrip, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)
    EnergyPerDistanceReturn = energy_requirements(RequiredSlipReturn, MassReturnTrip, velocity, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)  # [J/m]
    EnergyReturn = EnergyPerDistanceReturn * DistanceToTravel  # [J]

    # Round trip
    EnergyRoundTrip = EnergyOutward + EnergyReturn  # [J] for a mass of regolith of 20 [kg] (max load), and a distance of 1 [km]
    Beta = EnergyRoundTrip/(mRegolith*DistanceToTravel*3600)  # [kWh/kg/km]
    #print("The total energy required for the round trip is:", EnergyRoundTrip/1000, "[kJ] (", EnergyOutward/1000, "[kJ] for the the outward trip,", EnergyReturn/1000, "[kJ] for the return trip).")
    #print("Hence, beta =", EnergyRoundTrip/(mRegolith*DistanceToTravel*3600), "[kWh/kg/km]")  # Looks good since Dorian got 5.4 [kJ/kg/km] which might be a slight over-estimate.
    
    # N.B.: transporting 40 kg of regolith requires:
            # either two rovers to transport it at the same time (energy required times 2),
            # either the same rover to perform two round trips (energy required times 2 because the distance is multiplied by 2).

    # Graph as a function of regolith mass (assume a distance of 1 [km])
    """mRegolithVector = np.arange(0, 50, 1) * 90  # [kg]
    plt.plot(mRegolithVector, Beta * mRegolithVector, label='Assumptions: \n ')
    plt.ylabel('Locomotion energy requirements [kJ/km]')
    plt.xlabel('Regolith mass [kg]')
    plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
    plt.show()  # it has the expected shape!"""


if continueBoolean4 == 1:
    # Here, I need to ask the variation of soil density we are looking at (to determine whether we can neglect it in the computation of transport).
    # ATTENTION: apparently, the soil density has no influence on beta... Logical since it doesn't appear in the equations BUT to control!!!

    Depth = 1  # [m]
    HoleRadius = 2  # [m]
    SlopeHole = math.atan(Depth/HoleRadius)  # [rad]
    print("The slope in the hole is:", SlopeHole*180/math.pi, "[°]")

    velocity = 0.2  # [m/s], see "Regolith Advanced Surface Systems Operations Robot (RASSOR)"
    motor_efficiency = 0.6  #
    mRover = 66  # [kg] Total mass of RASSOR2 + mass of regolith carried
    mRegolith = 20  # [kg] The rover is assumed to transport the maximum amount of regolith each time.
    WheelRadius = 0.4318 / 2  # [m] Wheel radius (17 inches for the full wheel)
    WheelWidth = 0.1  # [m] Wheel width (from the picture?)
    wheelbase = 0.5  # [m]
    heightCOG = 0.1  # [m]
    Slope = 0  # [rad]

    DistanceToTravel = 1000  # [m]
    MassOutwardTrip = mRover
    MassReturnTrip = mRover + mRegolith

    # Outward trip: flat segment (without regolith)
    RequiredSlipOutwardFlat = slip_required(MassOutwardTrip, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)
    EnergyPerDistanceOutwardFlat = energy_requirements(RequiredSlipOutwardFlat, MassOutwardTrip, velocity, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)  # [J/m]
    EnergyOutwardFlat = EnergyPerDistanceOutwardFlat * (DistanceToTravel - HoleRadius)  # [J]

    # Outward trip: hole (without regolith)
    RequiredSlipOutwardHole = slip_required(MassOutwardTrip, WheelWidth, WheelRadius, -SlopeHole, wheelbase, heightCOG)
    EnergyPerDistanceOutwardHole = energy_requirements(RequiredSlipOutwardHole, MassOutwardTrip, velocity, WheelWidth, WheelRadius, SlopeHole, wheelbase, heightCOG)  # [J/m]
    EnergyOutwardHole = EnergyPerDistanceOutwardHole * (HoleRadius/math.cos(SlopeHole))  # [J]

    # Return trip: hole (with regolith)
    RequiredSlipReturnHole = slip_required(MassReturnTrip, WheelWidth, WheelRadius, SlopeHole, wheelbase, heightCOG)
    EnergyPerDistanceReturnHole = energy_requirements(RequiredSlipReturnHole, MassReturnTrip, velocity, WheelWidth, WheelRadius, SlopeHole, wheelbase, heightCOG)  # [J/m]
    EnergyReturnHole = EnergyPerDistanceReturnHole * (HoleRadius/math.cos(SlopeHole))  # [J]

    # Return trip: flat (with regolith)
    RequiredSlipReturnFlat = slip_required(MassReturnTrip, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)
    EnergyPerDistanceReturnFlat = energy_requirements(RequiredSlipReturnFlat, MassReturnTrip, velocity, WheelWidth, WheelRadius, Slope, wheelbase, heightCOG)  # [J/m]
    EnergyReturnFlat = EnergyPerDistanceReturnFlat * (DistanceToTravel - HoleRadius)  # [J]

    # Round trip
    EnergyOutward = EnergyOutwardFlat + EnergyOutwardHole
    EnergyReturn = EnergyReturnFlat + EnergyReturnHole
    EnergyRoundTrip = EnergyOutwardFlat + EnergyOutwardHole + EnergyReturnFlat + EnergyReturnHole  # [J] for a mass of regolith of 20 [kg] (max load), and a distance of 1 [km]
    Beta = EnergyRoundTrip/(mRegolith*DistanceToTravel)  # [J/kg/m] or [kJ/kg/km]
    print("The total energy required for the round trip is:", EnergyRoundTrip/1000, "[kJ] (", EnergyOutward/1000, "[kJ] for the the outward trip,", EnergyReturn/1000, "[kJ] for the return trip).")
    print("Hence, beta =", EnergyRoundTrip/(mRegolith*DistanceToTravel), "[kJ/kg/km]")  # Looks good since Dorian got 5.4 [kJ/kg/km] which might be a slight over-estimate.