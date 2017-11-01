from math import *
import random
from scipy.stats import norm
from scipy import signal
#from scipy.fftpack import fft2,ifft2,fftshift
from scipy.fftpack import fftn,ifftn,fftshift
import scipy.stats as stats
import numpy as np
import Sequences
import troot as tr
from ROOT import *
from copy import deepcopy
import array

# aliases
SH  = np.fft.fftshift
FT  = lambda x:  ( fftn( x ) ) #         Fourier transform
IFT = lambda x:  fftshift( ifftn( x )  ) # Inverse Fourier transform

gunits = 11
gunit=1.
gdsigma = 1.
gesigma = 1.5

#gelsigma = 5.
#gdifsigma = 5.

#--- generate grid
#
# Assume that the active volume is located at (0,0,0) 
#                     and have a length (20.,20.,20.) cm
#					  grid is 2x2x2 mm2
#


#---- arrays defauls

def arnull(ndim=2,size=10):
	""" create a ndim with nbins null array 
	"""
	kk = [0.]*size
	if (nn==3): kk = [kk,]*size
	ar = np.array([kk,]*size)
	return ar

def argaus2d(ndim=2,size=10,sigma=1.):
	""" create a nn dimensional array with a gausian of sigma which mean is in the cente of the array
	"""
	xsigma = sigma/float(size)
	ar = arnull(ndim,size)
	i0 = int(size/2)
	for i in range(size):
		x = norm.pdf(float(i-i0)/xsigma)
		for j in range(size):
			y = norm.pdf(float(j-i0)/xsigma)
			for k in range(size):
				z = norm.pdf(float(k-i0)/xsigma)
				ar[i][j][k] = x*y*z
	return ar

def argaus2d(sigma=1.):
	xsigma = sigma/gunit
	ar = arnull(nn=2)
	i0 = int(gunits/2)
	for i in range(gunits):
		x = norm.pdf(float(i-i0)/xsigma)
		for j in range(gunits):
			y = norm.pdf(float(j-i0)/xsigma)
			ar[i][j] = x*y
	return ar

def arnoise(ar,noise=0,poisson=False):
	br = deepcopy(ar)
	nx,ny = ar.shape
	for i in range(nx):
		for j in range(ny):
			mm = ar[i][j]
			if (noise>0):
				br[i][j] = mm + stats.poisson.rvs(noise)
			if (poisson): 
				if (mm>0):
					br[i][j] = stats.poisson.rvs(mm)
	return br

def arclean(ar,nlevel=0,nvalue=0.):
	""" 
	"""
	br = deepcopy(ar)
	nx,ny = ar.shape
	for i in range(nx):
		for j in range(ny):
			if (br[i][j]<=nlevel): br[i][j]=nvalue
	return br

#------ fourier utilities

def arlast(ar):
	nxs = ar.shape
	br = deepcopy(ar)
	cr = deepcopy(ar)
	#print ' br ',cr
	n = len(ar)
	br[0]=cr[-1]
	for i in range(n-1):
		br[i+1]=cr[i]
	cr = deepcopy(br)
	#print 'br (i) ',br
	if (len(nxs)<=1): return br
	for i in range(n):
		br[i][0]=cr[i][-1]
		for j in range(n-1):
			br[i][j+1]=cr[i][j]
	if (len(nxs)<=2): return br
	cr = deepcopy(br)
	#print ' br (j) ',br
	for i in range(n):
		for j in range(n):
			br[i][j][0] = cr[i][j][-1]
			for k in range(n-1):
				br[i][j][k+1] = cr[i][j][k]
	#print ' br (k) ',br
	return br 

def armesfourier(sig,res):
	fsig = FT(sig)
	fres = FT(res)
	fmes = Sequences.Multiply(fsig,fres)
	mes = IFT(fmes)
	mes = Sequences.Apply(mes,lambda x: x.real)
	mes = arlast(mes)
	return mes

def armesconvolve(sig,res):
	mes = signal.convolve2d(sig,res,mode='same')
	return mes

def arrecfourier(mes,res):
	fmes = FT(mes)
	fres = FT(res)
	fresinv = Sequences.Apply(fres,lambda x: 1/x)
	fsig = Sequences.Multiply(fmes,fresinv)
	sig = IFT(fsig)
	sig = Sequences.Apply(sig,lambda x: x.real)
	return sig

def arrecfouriercut(mes,res,gamma=0.):
	fmes = FT(mes)
	fres = FT(res)
	def ffilter(x):
		if (abs(x)>gamma):
			y = 1/x
		else:
			y = abs(x)/x*(1/gamma)
		return y
	fresinv = Sequences.applyto(fres,ffilter)
	fsig = Sequences.Multiply(fmes,fresinv)
	sig = IFT(fsig)
	sig = Sequences.Apply(sig,lambda x: x.real)
	return sig

#-------- plotting tools ---

def thisto3dar(ar):
	nx,ny,nz = len(ar),len(ar[0]),len(ar[0][0])
	h = TH3F('thisto3d','thisto3d',nx,1.*0,1.*nx,ny,0,1.*ny,nz,0.,1.*nz)
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				h.Fill(i+0.5,j+0.5,k+0.5,ar[i][j][k])
	return h

def thisto4dar(ar):
	n = len(ar)
	xs,ys,zs,ts = [],[],[],[]
	for i in range(n):
		for j in range(n):
			for k in range(n):
				xs.append(i)
				ys.append(j)
				zs.append(k)
				ts.append(ar[i][j][k])
	t = thisto4dlist(xs,ys,zs,ts)
	return t

def thisto4dlist( x, y, z, t, markerstyle = 20, markersize = 1 ):
    '''
        Plot a 3D dataset (x,y,z) with an extra color coordinate (t).
    '''
    data = array.array( 'd', [0.] * 4 )
    tree = TTree('DummyTree','DummyTree')
    tree.Branch('xyzt', data, 'x/D:y:z:t')

    for datai in zip(x,y,z,t):
        data[0], data[1], data[2], data[3] = datai
        tree.Fill()
    tree.SetMarkerStyle( markerstyle )
    tree.SetMarkerSize( markersize )
    tree.Draw('x:y:z:t','','zcol')
    return tree

def plot4d(ars,nx=2,ny=2):
	c = TCanvas()
	c.Divide(nx,ny)
	ts = []
	for i in range(len(ars)):
		c.cd(i+1)
		t = thisto4dar(ars[i])
		ts.append(t)
	return c,ts
