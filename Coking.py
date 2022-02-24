from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, \
    Viewer, input
from fipy.tools import numerix
from builtins import range

nx = 50  # L
dx = 1.
mesh = Grid1D(nx=nx, dx=dx)

phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=0.)

D = 1.  # k/Ro*Cp

valueLeft = 1300  # temp po lewej
valueRight = 290  # temp po prawej

phi.constrain(valueRight, mesh.facesRight)  # zostaw
phi.constrain(valueLeft, mesh.facesLeft)  # zostaw

eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D)  # zostaw tylko D bedzie trzeba

""" Warunek stabilnosci """
timeStepDuration = 0.9 * dx ** 2 / (2 * D)
steps = 100  # maksymalna liczba kroków

phiAnalytical = CellVariable(name="analytical value",
                             mesh=mesh)  # zostaw

if __name__ == '__main__':
    viewer = Viewer(vars=(phi, phiAnalytical),
                    datamin=0., datamax=1400.)  # min i max to wartosci y
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

for step in range(steps):  # do tej petli wrzucic zeby zmienial sie warunek brzegowy
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

print(phi.allclose(phiAnalytical, atol = 7e-4))

from fipy import input
if __name__ == '__main__':
    input("Explicit transient diffusion. Press <return> to proceed...")