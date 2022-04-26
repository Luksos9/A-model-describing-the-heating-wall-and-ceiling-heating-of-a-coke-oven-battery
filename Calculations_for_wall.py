""" Temperatures """

T_initial = 300  # [K]
T_end = 1300  # [K]
T_gas = 1800  # [K]

""" Data """

A = 1  # area [m^2]
L1, L2 = 2 * 0.15, 2 * 0.3  # Thickness [m]
L = L1 + L2  # Total Thickness of both layers
h = 1000  # [W/(m^2*K)] - Heat transfer coefficient

""" Combined coefficients"""

density_silica, density_coal = 2340, 1900  # kg / (m^3)
density_combined = density_silica * 1 / 3 + density_coal * 2 / 3
specific_heat_capacity_silica, specific_heat_capacity_coal = 920, 1800  # J/kg * K
specific_heat_capacity_combined = specific_heat_capacity_silica * 1 / 3 + specific_heat_capacity_coal * 2 / 3
K1, K2 = 0.5, 2  # Thermal conductivities of Silica brick (1) and Coal/coke bed (2)
K_combined = K1 * 1 / 3 + K2 * 2 / 3  # as 0.15 is 1/3 of 0.45 and 0.3 is 2/3 of 0.45 (total length)

a = K_combined / (density_combined * specific_heat_capacity_combined)  # thermal diffusivity

Bi = (h * L) / K_combined  # Way more than 0.1 so no lumped heat transfer, we need to use Fourier

dimensionless_temperature = (T_end - T_gas) / (T_initial - T_gas)  # 0. 4

""" After knowing dimensionless_temperature and 1/Bi we can determine Fourier Number (Fo) """
Fo = 0.3  # value rode from charts

# Fo = (a * t) / (L ** 2) => t = (Fo * (L ** 2)) / a

def calc_t(Fo, L, a):
    return (Fo * (L ** 2)) / a

final = calc_t(Fo, L, a)
print("1300 K will be reached after: {} hours".format(final / 3600))

