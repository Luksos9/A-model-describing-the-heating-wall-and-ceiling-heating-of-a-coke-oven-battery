from sympy import symbols, solve, Eq
import numpy

""" Initial conditions"""

a_solar = 0.7  # dimensionless
E_sky = 0.8  # dimensionless
a_sky = 0.8  # dimensionless
A = 1  # [m2] - area
T_comb, T1 = 1300, 1300  # [K] - Temperature of coke gas (Also it's T1 - temperature at the surface (gas-ceiling)
b = 5.67 * 10 ** (-8)  # [W/(m^2*K^4 )]
hs = 20  # [W/(m^2*K)]
h1 = 1000  # [W/(m^2*K)]
La = 1  # m
Ka = 0.1
Lb = 1
Kb = 0.1
Ta = 308  # [K]
Gs = 1373  # [W/m^2]
cos_o = (2 ** 1/2) / 2  # degrees
T_sky = 280  # [K]
T_sun = 5700  # [K]
T_surface = 300  # [K] - Ceiling surface temperature (one we want to calculate)

""" Energy balance components"""
Ra, Rb = La / Ka, Lb / Kb
R_total = 3.28
print(R_total)

E_cond = (T_comb - T_surface) / R_total  # Ts is unknown

E_conv = hs * (T_surface - Ta)

G_solar = Gs * cos_o
E_solar_abs = a_solar * G_solar

G_sky = b * T_sky ** 4
E_sky_abs = a_sky * G_sky

E_emit = E_sky * b * T_surface ** 4

for i in range(20):
    T_surface_upper = (A * R_total * a_solar * Gs * cos_o +
                       A * R_total * E_sky * b * T_sky ** 4 + T_comb + A * R_total * hs * Ta)
    T_surface_lower = (A * R_total * E_sky * b * T_surface ** 3 + A * R_total * hs + 1)
    T_surface_i = T_surface_upper / T_surface_lower
    print('T_surface_i: {}, T_surface: {}'.format(T_surface_i, T_surface))
    T_surface = T_surface_i * 0.2 + T_surface * 0.8
    print('T_surface - T_surface_i: {}'.format((T_surface - T_surface_i) ** 2))
    #print("Iteration {}: {:.3f} K".format(i, T_surface))


Ts = (-T_comb - (hs * Ta * R_total) + (a_solar * G_solar * R_total) + (a_sky * R_total * b * T_sky ** 4)) / (
            -1 - (hs * R_total) - (E_sky * b * R_total))


