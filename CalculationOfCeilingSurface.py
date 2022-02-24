from sympy import symbols, solve, Eq
import numpy

""" Initial conditions"""

a_solar = 0.7  # dimensionless
E_sky = 0.8  # dimensionless
a_sky = 0.8  # dimensionless
A = 1  # [m2] - area
T_comb, T1 = 1300, 1300  # [K] - Temperature of coke gas (Also it's T1 - temperature at the surface (gas-ceiling)
b = 5.67 * 10 ** (-8)  # [W/(m^2*K^4 )]
hs = 25  # [W/(m^2*K)]
Ta = 308  # [K]
Gs = 1373  # [W/m^2]
cos_o = (2 ** 1 / 2) / 2  # degrees
T_sky = 280  # [K]
T_sun = 5700  # [K]
T_surface = 300  # [K] - Ceiling surface temperature (one we want to calculate)

""" Total Ceiling Resistance made of 6 layers """
L1, L2, L3, L4, L5, L6 = 0.27, 0.25, 0.6, 0.45, 0.2, 0.06  # [m] - Length
k1, k2, k3, k4, k5, k6 = 0.5, 0.9, 0.9, 0.32, 0.6, 1
R1, R2, R3, R4, R5, R6 = L1 / k1, L2 / k2, L3 / k3, L4 / k4, L5 / k5, L6 / k6
R_total = 20
print('R total: {}'.format(R_total))

""" Energy balance """

# E_cond - E_conv + E_solar_abs + E_sky_abs - E_emit = 0
# E_cond + E_solar_abs + E_sky_abs = E_conv + E_emit


""" Energy balance components"""

# E_cond = (T_comb - T_surface) / R_total
#
# E_conv = hs * (T_surface - Ta)
#
# G_solar = Gs * cos_o
# E_solar_abs = a_solar * G_solar
#
# G_sky = b * T_sky ** 4
# E_sky_abs = a_sky * G_sky
#
# E_emit = E_sky * b * T_surface ** 4

for i in range(10):
    print('new T surface: {}'.format(T_surface))
    T_surface_upper = (A * R_total * a_solar * Gs * cos_o +
                       A * R_total * E_sky * b * T_sky ** 4 + T_comb + A * R_total * hs * Ta)
    T_surface_lower = (A * R_total * E_sky * b * T_surface ** 3 + A * R_total * hs + 1)

    T_surface_i = T_surface_upper / T_surface_lower
    print('T_surface - T_surface_i: {}'.format((T_surface - T_surface_i) ** 2))
    T_surface = T_surface_i
    print("Iteration {}: T_surface: {:.3f} K".format(i + 1, T_surface))
