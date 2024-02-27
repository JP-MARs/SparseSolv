import os, sys
os.system('taskkill.exe /f /im i_view64.exe')

import SparseSolvPy
from netgen.occ import *
from ngsolve import *
from numpy import *

cyl = Cylinder(Pnt(0,0,0), Z, r=1.0, h=1.0).bc("S_side")
cyl.mat("cyl")
cyl.faces.Min(Z).name="S_bottom"
cyl.faces.Max(Z).name="S_top"

geo = OCCGeometry(cyl)
mesh = Mesh(geo.GenerateMesh(maxh=0.2))
mesh.Curve(4)
mesh.GetMaterials(), mesh.GetBoundaries()

order = 1
fes = HCurl(mesh, order=order, dirichlet="S_side|S_bottom|S_top", type1=True, nograds = False)
print ("Hcurl_ndof =", fes.ndof)

u = fes.TrialFunction()
v = fes.TestFunction()

r = 1.0
h = 0.2
mu = 4*pi*1e-7
J = CoefficientFunction((0,0,1))

a = BilinearForm(fes)
a += 1/mu*curl(u)*curl(v)*dx
a += 1e-6/mu*u*v*dx

f = LinearForm(fes)
f += v*J * dx

gfA = GridFunction(fes)

with TaskManager():
	c = Preconditioner(a, type="local")
	a.Assemble()
	f.Assemble()

xx = linspace(-r,r,100)
solvers.CG(sol=gfA.vec, rhs=f.vec, mat=a.mat, pre=c.mat, tol=1e-12, printrates=False, maxsteps=1000)
H = curl(gfA)/mu
Hy1 = array([H[1](mesh(xi, 0.0, 0.1)) for xi in xx])

import scipy.sparse as sp
A = sp.csr_matrix (a.mat.CSR())
Acut = A[:,fes.FreeDofs()][fes.FreeDofs(),:]
fcut = array(f.vec.FV())[fes.FreeDofs()]
ucut = array(f.vec.FV(), copy=True)[fes.FreeDofs()]

rows, cols = Acut.nonzero()
vals = Acut[rows, cols]
vals = ravel(vals)
mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals)

tol = 1e-8
max_iter = 2000
acc = 1.15

is_log = True
is_save_best = True
is_diag_scale = True
DivergeType = 1
BadDivVal = 100.0
BadDivCount = 100
log = SparseSolvPy.MatSolversICCG.solveICCGPy(len(fcut), tol, max_iter, acc, mat, fcut, ucut, True)

print(type(log))
for arg in log:
	print(arg)

array(gfA.vec.FV(), copy=False)[fes.FreeDofs()] = ucut

H = curl(gfA)/mu
Hy2 = array([H[1](mesh(xi, 0.0, 0.1)) for xi in xx])

from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})
plt.figure(figsize=(3,3),dpi=400)
ax = plt.axes([0.20,0.15,0.75,0.75])
ax.set_xlim(-r,r)
ax.set_ylim(-0.6,0.6)
ax.set_xlabel('${\\it x}$ (mm)',fontname='times new roman')
ax.set_ylabel('${\\it y}$ (mm)',fontname='times new roman')

ax.plot(xx,xx/2,'k-',linewidth=0.5);
ax.plot(xx,Hy1,'r--',linewidth=0.5);
ax.plot(xx,Hy2,'b:',linewidth=0.5);

ax.xaxis.grid(True, which='major', linestyle=':', linewidth=0.7)
ax.xaxis.grid(True, which='minor', linestyle=':', linewidth=0.5)
ax.yaxis.grid(True, which='major', linestyle=':', linewidth=0.5)
ax.yaxis.grid(True, which='minor', linestyle=':', linewidth=0.5)
ax.tick_params(left=True, right=True, top=True, bottom=True, which="major", direction='in', length=2 ,width=0.4, pad=3)

plt.setp(ax.get_xticklabels(),fontname='times new roman')
plt.setp(ax.get_yticklabels(),fontname='times new roman')

ax.legend(['Analytic', 'CG Solver', 'Shifted ICCG'], frameon=False, bbox_to_anchor=(0.35,0.85),loc='center',framealpha=0.0,fancybox='none',edgecolor='none',prop={"family":'Times New Roman',"size":8})


FileName = os.path.basename(__file__).replace(".py",".png")
plt.savefig(FileName)
os.startfile(FileName)
