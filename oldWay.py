""" Initial Conditions """
A = 1  # [m2] - area
T_comb, T1 = 1300, 1300  # [K] - Temperature of coke gas (Also it's T1 - temperature at the surface (gas-ceiling)
h1 = 25  # [W/(m^2*K)]
hrad = 6.05  # [W/(m^2*K)]
T_surrounding = 308.15  # [K]
# T_surface = 300  # [K] - Ceiling surface temperature (one we want to calculate)

""" Total Ceiling Resistance made of 6 layers """
L1, L2, L3, L4, L5, L6 = 0.27, 0.25, 0.6, 0.45, 0.2, 0.06  # [m] - Length
k1, k2, k3, k4, k5, k6 = 0.85, 0.9, 0.9, 0.32, 1, 1
R1, R2, R3, R4, R5, R6, R_conv, R_rad = L1 / k1, L2 / k2, L3 / k3, L4 / k4, L5 / k5, \
                                        L6 / k6, 1 / h1, 1 / hrad  # A is equal to 1 m, so not included

""" Two parallel resistances R_rad and R_conv can be replaced by an equivalent resistance R_equiv """

R_equiv = R_conv * R_rad / (R_rad + R_conv)

R_total = R1 + R2 + R3 + R4 + R5 + R6 + R_equiv
print('R total: {}'.format(R_total))

""" After calculating Rtotal we know values of all the variables that are needed to
 determine rate of heat transfer (Equation 3.1) for our case: """

Q = (T1 - T_surrounding) / R_total

""" Now knowing Q we can determine temperatures from T2 – T7 using Equation 3.5: """

T2 = T1 - (Q * R1)
T3 = T2 - (Q * R2)
T4 = T3 - (Q * R3)
T5 = T4 - (Q * R4)
T6 = T5 - (Q * R5)
T7 = T6 - (Q * R6)

print('T1: {}, T2: {}, T3: {}, T4: {}, T5: {}, T6: {}, T7: {}'.format(T1, T2, T3, T4, T5, T6, T7))
