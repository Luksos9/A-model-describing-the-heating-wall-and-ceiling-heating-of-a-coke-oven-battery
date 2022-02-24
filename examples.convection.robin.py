from fipy import CellVariable, FaceVariable, Grid1D, DiffusionTerm, PowerLawConvectionTerm, ImplicitSourceTerm, Viewer
from fipy.tools import numerix
nx = 100
dx = 1.0 / nx

mesh = Grid1D(nx=nx, dx=dx)
C = CellVariable(mesh=mesh)
D = 2.0
P = 3.0
C.faceGrad.constrain([0], mesh.facesRight)

convectionCoeff = FaceVariable(mesh=mesh, value=[P])
convectionCoeff[..., mesh.facesLeft.value] = 0.
diffusionCoeff = FaceVariable(mesh=mesh, value=1.)
diffusionCoeff[..., mesh.facesLeft.value] = 0.

eq = (PowerLawConvectionTerm(coeff=convectionCoeff)
      == DiffusionTerm(coeff=diffusionCoeff) - ImplicitSourceTerm(coeff=D)
      - (P * mesh.facesLeft).divergence)
A = numerix.sqrt(P**2 + 4 * D)
x = mesh.cellCenters[0]
CAnalytical = CellVariable(mesh=mesh)
CAnalytical.setValue(2 * P * numerix.exp(P * x / 2)
                     * ((P + A) * numerix.exp(A / 2 * (1 - x))
                        - (P - A) * numerix.exp(-A / 2 *(1 - x)))
                     / ((P + A)**2*numerix.exp(A / 2)
                        - (P - A)**2 * numerix.exp(-A / 2)))

if __name__ == '__main__':
    C.name = '$C$'
    CAnalytical.name = '$C_{analytical}$'
    viewer = Viewer(vars=(C, CAnalytical))
if __name__ == '__main__':
    restol = 1e-5
    anstol = 1e-3
else:
    restol = 0.5
    anstol = 0.15
res = 1e+10
while res > restol:
    res = eq.sweep(var=C)
    if __name__ == '__main__':
        viewer.plot()
print(C.allclose(CAnalytical, rtol=anstol, atol=anstol))