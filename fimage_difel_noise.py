from math import *
import random
from scipy.stats import norm
from scipy import signal
#from scipy.fftpack import fft2,ifft2,fftshift
from scipy.fftpack import fftn,ifftn,fftshift
from scipy.stats import poisson
import numpy as np
import Sequences
import troot.thisto as thisto
from ROOT import *
from copy import deepcopy
import array

import arimage as ari

gsigma=1.
gelsigma=1.
gpoisson=True
gnoise=2
ggamma=0.1

def diffusion(sigma=1.,npoints=10,poisson=False):
	res3d = ari.arres3d(sigma)
	#print ' res ',res3d
	sig3d = ari.arsig3d(npoints)
	#print ' sig ',sig3d
	mes3d = ari.armesfourier(sig3d,res3d)
	#print ' mes ',mes3d
	#rec3d = arrecfourier(mes3d,res3d)
	mes3dn = []
	nn = len(mes3d)
	for i in range(nn):
		ar = mes3d[i]
		arn = ari.arnoise(ar,poisson=poisson)
		mes3dn.append(arn)
	mes3dn = np.array(mes3dn)
	return res3d,sig3d,mes3dn

def el(mes3d,sigma=1.,noise=0.):
	print ' sigma ',sigma
	res2d = ari.arres2d(sigma)
	print ' res2d !!!'
	elmes = []
	elmesn = []
	nn = len(mes3d)
	for i in range(nn):
		ielmesn = elplane(mes3d[i],res2d,noise)
		elmesn.append(ielmesn)
	#elmes = np.array(ielmes)
	elmesn = np.array(elmesn)
	return res2d,elmesn

def elplane(mes2d,res2d,noise=0.):
	mes2d = ari.armesfourier(mes2d,res2d)
	mes2d = ari.arnoise(mes2d,noise)
	return mes2d

def generate():
	res3d,sig3d,dif3d = diffusion(sigma=gsigma,poisson=gpoisson)
	res2d,el3d = el(dif3d,sigma=gelsigma,noise=gnoise)
	return res3d,res2d,sig3d,dif3d,el3d

def elrecoplane(mes2d,res2d,gamma=0.):
	rec2d = ari.arrecfouriercut(mes2d,res2d,gamma=gamma)
	return rec2d	

def elreco(mes3d,res2d,gamma=0.):
	elrec3d = []
	nn = len(mes3d)
	for i in range(nn):
		ielrec = elrecoplane(mes3d[i],res2d,gamma)
		elrec3d.append(ielrec)
	elrec3d = np.array(elrec3d)
	x0 = Sequences.mapA3(min,elrec3d)
	print ' x0 ',x0
	#if (x0<0):
	#	x0add = lambda x: x-x0
	#	elrec3d = Sequences.applyto(elrec3d,x0add)
	#	x0 = Sequences.mapA3(min,elrec3d)
	#	print ' x0 add ',x0
	return elrec3d

def difreco(mes3d,rec3d,gamma=0.):
	rec3d = ari.arrecfouriercut(mes3d,res3d,gamma=gamma)
	x0 = Sequences.mapA3(min,rec3d)
	print ' x0 ',x0
	#if (x0<0):
	#	x0add = lambda x: x-x0
	#	rec3d = Sequences.applyto(rec3d,x0add)
	#	x0 = Sequences.mapA3(min,rec3d)
	#	print ' x0 add ',x0
	return rec3d

def reco(mes3d,res3d,res2d,gamma=0.):
	recel3d = elreco(mes3d,res2d,gamma=ggamma)
	recdif3d = difreco(recel3d,res3d,gamma=ggamma)
	return recdif3d,recel3d


res3d,res2d,sig3d,dif3d,el3d = generate()
recdif3d,recel3d = reco(el3d,res3d,res2d)

print " diffusion + EL "
c,ts = ari.plot4d([sig3d,dif3d,el3d,recel3d,recdif3d],nx=3)
raw_input('enter key')


def plotslices(xs,ms,rs,zs):
	ncon=20
	def ranges(xx):
		x0 = Sequences.mapA3(min,xx)
		xf = Sequences.mapA3(max,xx)
		print 'ranges ',x0,xf
		con = np.linspace(x0,xf,ncon)
		return con
	nn=len(xs)
	hxs = map(thisto.thisto2darray,xs)
	hms = map(thisto.thisto2darray,ms)
	hrs = map(thisto.thisto2darray,rs)
	hzs = map(thisto.thisto2darray,zs)
	t = TCanvas()
	t.Divide(nn,4)
	ii=0
	def iplot(h,ii):
		ii=ii+1
		t.cd(ii)
		h.SetContour(ncon,con)
		h.Draw('col')
		return ii
	con = ranges(xs)
	for h in hxs: ii=iplot(h,ii)
	con = ranges(ms)
	for h in hms: ii=iplot(h,ii)
	con = ranges(rs)
	for h in hrs: ii=iplot(h,ii)
	con = ranges(zs)
	for h in hzs: ii=iplot(h,ii)
	return t,hxs,hms,hrs,hzs

print 'slices '
cs = plotslices(sig3d,el3d,recel3d,recdif3d) 
raw_input('enter_key')
