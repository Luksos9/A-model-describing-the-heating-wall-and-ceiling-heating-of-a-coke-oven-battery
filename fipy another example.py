import fipy as fp

nx = 100.
Lx = 1.
mesh = fp.Grid1D(dx=0.01, nx=100.)
C = fp.CellVariable(mesh=mesh, name='C')
A = 1.0
D = 2.0
P = 3.0

C.faceGrad.constrain(0, mesh.facesRight)
DiffFace = fp.FaceVariable(mesh=mesh, value=A)
DiffFace.setValue(0., where=mesh.facesLeft)

x = mesh.faceCenters
convValue = (P, 0)
Convface = fp.FaceVariable(mesh=mesh, value=convValue)
Convface.setValue(0., where=mesh.facesLeft)
mask = (mesh.facesLeftmesh.faceNormals).divergence
Af = mesh._cellAreas[mesh.facesLeft][0]
dPR = mesh._cellDistances[mesh.facesLeft][0]
mask = maskAAf / (1 + dPR)
eq = 0 == fp.DiffusionTerm(coeff=DiffFace, var=C) - \
fp.ExponentialConvectionTerm(Convface, var=C) - \
fp.ImplicitSourceTerm(D, var=C) + fp.ImplicitSourceTerm(coeff=Pmask, var=C) - P * mask
res = 1
restol = 1e-4
viewer = fp.Viewer(vars=C)
while res > restol:
    res = eq.sweep(var=C)
print('residual: %s' % res)
viewer.plot()
