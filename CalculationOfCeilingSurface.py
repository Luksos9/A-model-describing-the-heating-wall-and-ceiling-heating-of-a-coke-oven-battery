import matplotlib.pyplot as plt
import numpy as np

""" Initial conditions"""

a_solar = 0.7  # dimensionless
E_sky = 0.8  # dimensionless
a_sky = 0.8  # dimensionless
A = 1  # [m2] - area
T_comb, T1 = 1300, 1300  # [K] - Temperature of coke gas (Also it's T1 - temperature at the surface (gas-ceiling)
b = 5.67 * 10 ** (-8)  # [W/(m^2*K^4 )] - Stefan Boltzmann constant
hs = 25  # [W/(m^2*K)] - heat transfer coefficient
Ta = 308  # [K] - Temperature around surface (air)
Gs = 1373  # [W/m^2] - radiation
cos_o = (2 ** 1 / 2) / 2  # degrees
T_sky = 283.15  # [K] - sky temperature
T_sun = 5700  # [K] - sun temperature
T_surface = 320  # [K] - assumed ceiling surface temperature (one we want to calculate)

""" Total Ceiling Resistance made of 6 layers """
L1, L2, L3, L4, L5, L6 = 0.27, 0.25, 0.6, 0.45, 0.2, 0.06  # [m] - Length
k1, k2, k3, k4, k5, k6 = 0.5, 0.9, 0.9, 0.32, 0.6, 1  # [W/(m · K)] - Thermal conductivity
R1, R2, R3, R4, R5, R6 = L1 / k1, L2 / k2, L3 / k3, L4 / k4, L5 / k5, L6 / k6  # [K/W] - Resistance
R_total = R1 + R2 + R3 + R4 + R5 + R6

""" Energy balance """
# E_cond - E_conv + E_solar_abs + E_sky_abs - E_emit = 0

""" Energy balance components"""
# E_cond = (T_comb - T_surface) / R_total
# E_conv = hs * (T_surface - Ta)

# G_solar = Gs * cos_o
# E_solar_abs = a_solar * G_solar

# G_sky = b * T_sky ** 4
# E_sky_abs = a_sky * G_sky

# E_emit = E_sky * b * T_surface ** 4

""" Solution """
for i in range(10):
    print('new T surface: {}'.format(T_surface))
    T_surface_upper = (A * R_total * a_solar * Gs * cos_o +
                       A * R_total * E_sky * b * T_sky ** 4 + T_comb + A * R_total * hs * Ta)
    T_surface_lower = (A * R_total * E_sky * b * T_surface ** 3 + A * R_total * hs + 1)

    T_surface_i = T_surface_upper / T_surface_lower
    print('T_surface - T_surface_i: {}'.format((T_surface - T_surface_i) ** 2))
    T_surface = T_surface_i
    print("Iteration {}: T_surface: {:.3f} K".format(i + 1, T_surface))

""" Temperature distribution """
Q = (T1 - T_surface) / R_total
T2 = T1 - (Q * R1)
T3 = T2 - (Q * R2)
T4 = T3 - (Q * R3)
T5 = T4 - (Q * R4)
T6 = T5 - (Q * R5)

print(T1, T2, T3, T4, T5, T6, T_surface)

""" Temperature / thickness plot"""
xpoints = np.array(
    [0, L1, L1 + L2, L1 + L2 + L3, L1 + L2 + L3 + L4, L1 + L2 + L3 + L4 + L5, L1 + L2 + L3 + L4 + L5 + L6])
ypoints = np.array([T1, T2, T3, T4, T5, T6, T_surface])

y2points = [1300,
            1284.28901759065,
            1268.57803518131,
            1252.86705277196,
            1237.15607036261,
            1221.44508795327,
            1205.73410554392,
            1190.02312313457,
            1174.31214072523,
            1158.60115831588,
            1142.89017590653,
            1134.80839483177,
            1126.72661375700,
            1118.64483268223,
            1110.56305160747,
            1102.48127053270,
            1094.39948945793,
            1086.31770838317,
            1078.23592730840,
            1070.15414623363,
            1062.07236515887,
            1042.67609057943,
            1023.27981599998,
            1003.88354142054,
            984.487266841104,
            965.090992261664,
            945.694717682224,
            926.298443102784,
            906.902168523343,
            887.505893943903,
            868.109619364463,
            827.195602673456,
            786.281585982449,
            745.367569291442,
            704.453552600435,
            663.539535909428,
            622.625519218422,
            581.711502527415,
            540.797485836408,
            499.883469145401,
            458.969452454394,
            449.271315164674,
            439.573177874954,
            429.875040585234,
            420.176903295513,
            410.478766005793,
            400.780628716073,
            391.082491426353,
            381.384354136633,
            371.686216846913,
            361.988079557193,
            360.242414845043,
            358.496750132894,
            356.751085420744,
            355.005420708594,
            353.259755996445,
            351.514091284295,
            349.768426572146,
            348.022761859996,
            346.277097147846,
            344.531432435697]

x2points = [0,
            0.0270000000000000,
            0.0540000000000000,
            0.0810000000000000,
            0.108000000000000,
            0.135000000000000,
            0.162000000000000,
            0.189000000000000,
            0.216000000000000,
            0.243000000000000,
            0.270000000000000,
            0.295000000000000,
            0.320000000000000,
            0.345000000000000,
            0.370000000000000,
            0.395000000000000,
            0.420000000000000,
            0.445000000000000,
            0.470000000000000,
            0.495000000000000,
            0.520000000000000,
            0.580000000000000,
            0.640000000000000,
            0.700000000000000,
            0.760000000000001,
            0.820000000000001,
            0.880000000000001,
            0.940000000000001,
            1.00000000000000,
            1.06000000000000,
            1.12000000000000,
            1.16500000000000,
            1.21000000000000,
            1.25500000000000,
            1.30000000000000,
            1.34500000000000,
            1.39000000000000,
            1.43500000000000,
            1.48000000000000,
            1.52500000000000,
            1.57000000000000,
            1.59000000000000,
            1.61000000000000,
            1.63000000000000,
            1.65000000000000,
            1.67000000000000,
            1.69000000000000,
            1.71000000000000,
            1.73000000000000,
            1.75000000000000,
            1.77000000000000,
            1.77600000000000,
            1.78200000000000,
            1.78800000000000,
            1.79400000000000,
            1.80000000000000,
            1.80600000000000,
            1.81200000000000,
            1.81800000000000,
            1.82400000000000,
            1.83000000000000]

plt.title("Ceiling temperature distribution analytical and numeric")
plt.xlabel("Thickness [m]")
plt.ylabel("Temperature [K]")

plt.plot(x2points, y2points)
plt.plot(xpoints, ypoints, marker='o', ms=3)
plt.grid(linestyle='--')
plt.legend(['numeric', 'analytical'])


