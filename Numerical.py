from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, \
    Viewer
from fipy.tools import numerix
import matplotlib.pyplot as plt
import numpy as np

Ta = 308  # [K] - Temperature around surface (air)
T_comb, T1 = 1300, 1300  # [K] - Temperature of coke gas (Also it's T1 - temperature at the surface (gas-ceiling)

""" Lengths of ceiling's layers """
L1, L2, L3, L4, L5, L6 = 0.27, 0.25, 0.6, 0.45, 0.2, 0.06  # [m] - Layer's thickness
# to simplify calculations average thermal conductivity, density and
# specific heat capacity will be calculated. In order to do so their
# % contribution in total ceiling's thickness will be calculated

# 1. Total ceiling's thickness
L_total = L1 + L2 + L3 + L4 + L5 + L6
print(L_total)

""" Thermal conductivities of ceiling's layers """
k1, k2, k3, k4, k5, k6 = 0.5, 0.9, 0.9, 0.32, 0.6, 1

""" Resistances """
R1, R2, R3, R4, R5, R6 = L1 / k1, L2 / k2, L3 / k3, L4 / k4, L5 / k5, L6 / k6
R_total = R1 + R2 + R3 + R4 + R5 + R6

""" Densities of ceiling's layers """
Ro_1, Ro_2, Ro_3, Ro_4, Ro_5, Ro_6 = 680, 1000, 1000, 1200, 800, 1000

""" Heat Capacity of ceiling's layers """
C1, C2, C3, C4, C5, C6 = 2340, 1800, 1800, 1900, 1700, 2300

""" diffusion coefficients """
D1 = k1 / (Ro_1 * C1)  # k/Ro*Cp
D2 = k2 / (Ro_2 * C2)  # k/Ro*Cp
D3 = k3 / (Ro_3 * C3)  # k/Ro*Cp
D4 = k4 / (Ro_4 * C4)  # k/Ro*Cp
D5 = k5 / (Ro_5 * C5)  # k/Ro*Cp
D6 = k6 / (Ro_6 * C6)  # k/Ro*Cp

""" Temperature across layers chart """
Q = (T_comb - Ta) / R_total
T2 = T1 - (Q * R1)
T3 = T2 - (Q * R2)
T4 = T3 - (Q * R3)
T5 = T4 - (Q * R4)
T6 = T5 - (Q * R5)
T7 = T6 - (Q * R6)

""" Temperature / thickness plot"""
xpoints = np.array(
    [0, L1, L1 + L2, L1 + L2 + L3, L1 + L2 + L3 + L4, L1 + L2 + L3 + L4 + L5, L1 + L2 + L3 + L4 + L5 + L6])
ypoints = np.array([T1, T2, T3, T4, T5, T6, T7])

plt.title("Ceiling's temperature diversity")
plt.xlabel("Thickness [m]")
plt.ylabel("Temperature [K]")

plt.plot(xpoints, ypoints)
plt.grid(linestyle='--')
plt.show()

nx = 500
dx = L_total / nx  # dla roznych dxów

""" jak przejsc z aproxymacji 2 pochodnej do macierzy """
""" aproxymacja drugiej pochodnej """
""" macierz A * wektor X """  # Macierz A bedzie macierza kwadratowa o rozmiarze nx * nx
""" macierz tridiagonalna """

mesh = Grid1D(nx=nx, dx=dx)

phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=0.)

phiAnalytical = CellVariable(name="analytical value",
                             mesh=mesh)

L = nx * dx

D = FaceVariable(mesh=mesh, value=D1)
X = mesh.faceCenters[0]
D.setValue(D1, where=(0 <= X) & (X <= L1))
D.setValue(D2, where=(L1 < X) & (X <= L1 + L2))
D.setValue(D3, where=(L1 + L2 < X) & (X <= L1 + L2 + L3))
D.setValue(D4, where=(L1 + L2 + L3 < X) & (X <= L1 + L2 + L3 + L4))
D.setValue(D5, where=(L1 + L2 + L3 + L4 < X) & (X <= L1 + L2 + L3 + L4 + L5))
D.setValue(D6, where=(L1 + L2 + L3 + L4 + L5 < X) & (X <= L1 + L2 + L3 + L4 + L5 + L6))

valueLeft = 0
fluxRight = 1

phi = CellVariable(mesh=mesh)
phi.faceGrad.constrain([fluxRight], mesh.facesRight)
phi.constrain(valueLeft, mesh.facesLeft)

phi.setValue(0)
DiffusionTerm(coeff=D).solve(var=phi)

""" Analityczne dla D podstaw"""
x = mesh.cellCenters[0]
phiAnalytical.setValue(x)
phiAnalytical.setValue(10 * x - 9. * L / 4.,
                       where=(L / 4. <= x) & (x < 3. * L / 4.))
phiAnalytical.setValue(x + 18. * L / 4.,
                       where=3. * L / 4. <= x)
print(phi.allclose(phiAnalytical, atol=1e-8, rtol=1e-8))

"""phiAnalytical.setValue(D1,
                       where=(0 <= x) & (x <= L1))
phiAnalytical.setValue(D2,
                       where=(L1 < x) & (x <= L1 + L2))
phiAnalytical.setValue(D3,
                       where=(L1 + L2 < x) & (x <= L1 + L2 + L3))
phiAnalytical.setValue(D4,
                       where=(L1 + L2 + L3 < x) & (x <= L1 + L2 + L3 + L4))
phiAnalytical.setValue(D5,
                       where=(L1 + L2 + L3 + L4 < x) & (x <= L1 + L2 + L3 + L4 + L5))
phiAnalytical.setValue(D6,
                       where=(L1 + L2 + L3 + L4 + L5 < x) & (x <= L1 + L2 + L3 + L4 + L5 + L6))"""

""" drugi wykres z temperatura wprost (values wziete z rozwiazania analityznyego wczesniej)"""

# 1 . policzone  temp
# 2. narysowany wykres temp od x

from fipy import input

if __name__ == '__main__':
    Viewer(vars=(phi)).plot()
    input("Non-uniform steady-state diffusion. Press <return> to proceed...")
