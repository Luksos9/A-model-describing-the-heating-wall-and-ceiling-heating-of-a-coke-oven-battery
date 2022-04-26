from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, \
    Viewer, input
from fipy.tools import numerix
from builtins import range

""" Lengths of ceiling's layers """
L1, L2, L3, L4, L5, L6 = 0.27, 0.25, 0.6, 0.45, 0.2, 0.06  # [m] - Layer's thickness
# to simplify calculations average thermal conductivity, density and
# specific heat capacity will be calculated. In order to do so their
# % contribution in total ceiling's thickness will be calculated

# 1. Total ceiling's thickness
L_total = L1 + L2 + L3 + L4 + L5 + L6

# 2. Each layer contribution
L1_contribution = L1/L_total
L2_contribution = L2/L_total
L3_contribution = L3/L_total
L4_contribution = L4/L_total
L5_contribution = L5/L_total
L6_contribution = L6/L_total

""" Thermal conductivities of ceiling's layers """
k1, k2, k3, k4, k5, k6 = 0.5, 0.9, 0.9, 0.32, 0.6, 1
# 3. Total thermal conductivity
k_total = k1 * L1_contribution + k2 * L2_contribution + k3 * L3_contribution + \
    k4 * L4_contribution + k5 * L5_contribution + k6 * L6_contribution

""" Densities of ceiling's layers """
Ro_1, Ro_2, Ro_3, Ro_4, Ro_5, Ro_6 = 680, 1000, 1000, 1200, 800, 1000
# 4. Total density
Ro_total = Ro_1 * L1_contribution + Ro_2 * L2_contribution + Ro_3 * L3_contribution + \
    Ro_4 * L4_contribution + Ro_5 * L5_contribution + Ro_6 * L6_contribution

""" Heat Capacity of ceiling's layers """
C1, C2, C3, C4, C5, C6 = 2340, 1800, 1800, 1900, 1700, 2300
# 5. Total density
C_total = C1 * L1_contribution + C2 * L2_contribution + C3 * L3_contribution + \
    C4 * L4_contribution + C5 * L5_contribution + C6 * L6_contribution

nx = 50
dx = 1  # dx
mesh = Grid1D(nx=nx, dx=dx)
print(mesh)

phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=0.)


D = k_total / (Ro_total * C_total)  # k/Ro*Cp

valueLeft = 1300  # [K] - Left Temperature (where coke gas is burning)
valueMid = 300  # Right Temperature
valueRight = 330  # Right Temperature

phi.constrain(valueRight, mesh.facesRight)
phi.constrain(valueLeft, mesh.facesLeft)
phi.constrain(valueLeft, mesh.facesLeft)

eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D)  # zostaw tylko D bedzie trzeba

""" Stability condition """
timeStepDuration = 0.9 * dx ** 2 / (2 * D)
steps = 100  # maksymalna liczba kroków

phiAnalytical = CellVariable(name="analytical value",
                             mesh=mesh)  # zostaw

if __name__ == '__main__':
    viewer = Viewer(vars=(phi, phiAnalytical),
                    datamin=273., datamax=1400.)  # min i max to wartosci y
    viewer.plot()


x = mesh.cellCenters[0]  # wartosci w centrum siatki
t = timeStepDuration * steps  # maksymalny czas

""" Analityczne """
try:
    from scipy.special import erf

    phiAnalytical.setValue(1 - erf(x / (2 * numerix.sqrt(D * t))))
except ImportError:
    print("The SciPy library is not available to test the solution to \
the transient diffusion equation")

for step in range(steps):
    eqX.solve(var=phi,
              dt=timeStepDuration)
    if __name__ == '__main__':
        viewer.plot()

from builtins import range

for step in range(steps):
    eqX.solve(var=phi,
              dt=timeStepDuration)
    if __name__ == '__main__':
        viewer.plot()

print(phi.allclose(phiAnalytical, atol=7e-4))

from fipy import input

if __name__ == '__main__':
    input("Explicit transient diffusion. Press <return> to proceed...")
