from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, Viewer
from fipy.tools import numerix
from fipy import input

nx = 50
dx = 1.
mesh = Grid1D(nx=nx, dx=[0, 0.27, 0.25, 0.6, 0.45, 0.2])

phi = CellVariable(name="solution variable",
                   mesh=mesh,
                   value=0.)

D = 1.

valueLeft = 1
valueRight = 0

phi.constrain(valueRight, mesh.facesRight)
phi.constrain(valueLeft, mesh.facesLeft)

eqX = TransientTerm() == ExplicitDiffusionTerm(coeff=D)

timeStepDuration = 0.9 * dx**2 / (2 * D)
steps = 100 # krok powinnien byc raczej duzy (0.25 etc)


print(mesh)

phiAnalytical = CellVariable(name="analytical value",
                             mesh=mesh)
if __name__ == '__main__':
    viewer = Viewer(vars=(phi, phiAnalytical),
                    datamin=0., datamax=1.)
    viewer.plot()

x = mesh.cellCenters[0]
t = timeStepDuration * steps

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

print(phi.allclose(phiAnalytical, atol = 7e-4))

if __name__ == '__main__':
    input("Explicit transient diffusion. Press <return> to proceed...")



