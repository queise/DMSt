import sys
import numpy as np
import math as m
from scipy.integrate import odeint
from lib.constants import *
from lib.unitconv import *
from lib.functions import *

def F_calcNxnum(i, t, dm):
    """ Nx(t) is found solving dN/dt through odeint. The result is a list, we are interested in the [0]
        scipy.integrate.odeint(func, y0, t, args=(), Dfun=None, col_deriv=0, full_output=0, ml=None, mu=None, rtol=None, atol=None, tcrit=None,... """
    # Case 1: capture and annihilation, with rx(t)
    print '-- Nx(t):'
    tuple1 = dm.Capt_v[i], dm.anncs, t.rx, t.time                 # extra arguments for function dN1dt
    Nx_C_A = np.array([item[0] for item in odeint(dNdt_C_A, dm.Nxinit, t.time, args=(tuple1))])
    print "   ... Nx_C_A(t) calculated"
    # Case 2: capture, selfcapture and annihilation, with rx(t)
    Nx_C_sC_A = 0.
    if dm.SELFCAPT:
        tuple2 = dm.Capt_v[i], dm.anncs, t.rx, t.time, dm.selfCapt_v[i]         # extra arguments for function dNdt_C_sC_A_rxconst
        Nx_C_sC_A = np.array([item[0] for item in odeint(dNdt_C_sC_A, dm.Nxinit, t.time, args=(tuple2))])
        print "   ... Nx_C_sC_A(t) calculated"
        Nx = Nx_C_sC_A
    else:
        Nx = Nx_C_A
    return Nx, Nx_C_A, Nx_C_sC_A

def dNdt_C_A_rxconst(Nx, t, capt, ann): # return derivatives of the array Nx
    """ dN/dt = capt - ann*(N(t)**2) """
    return capt - ann*(Nx[0]**2) # np.array([ capt-ann*(Nx[0]**2) ])

def dNdt_C_A(Nx, t, capt, anncs, rxvect, time): # return derivatives of the array Nx
    """ dN/dt = capt - ann(t)*(N(t)**2) """
    # first obtains rx at t by interpolating:
    rx = np.interp(t,time,rxvect)
    # recalcs ann calling funtion for rx
    ann = F_ann(anncs,rx)
    return capt - ann*(Nx[0]**2)

def dNdt_C_sC_A_rxconst(Nx, t, capt, selfcapt, ann):
    """ dN/dt = capt + selfCapt*N(t) - ann*(N(t)**2) """
    return capt + selfcapt*Nx[0] - ann*(Nx[0]**2)

def dNdt_C_sC_A(Nx, t, capt, anncs, rxvect, time, selfcapt):
    """ dN/dt = capt + selfCapt*N(t) - ann(t)*(N(t)**2) """
    # first obtains rx at t by interpolating:
    rx = np.interp(t,time,rxvect)
    # recalcs ann calling funtion for rx
    ann = F_ann(anncs,rx)
    return capt + selfcapt*Nx[0] - ann*(Nx[0]**2)

def dNdt_C_sC_A_geoxx(Nx, t, capt, anncs, rxvect, time, selfcapt, Nx21geoxx, kgeoxx):
    """ dN/dt = capt + selfCapt*N(tgeo)*(rx(t)2/rxgeo2) - ann(t)*(N(t)**2) """
    # first obtains rx at t by interpolating:
    rx = np.interp(t,time,rxvect)
    # recalcs ann calling funtion for rx
    ann = F_ann(anncs,rx)
    # obtains rx at tgeoxx from vector:
    rxgeoxx = rxvect[kgeoxx]
    return capt + selfcapt*Nx21geoxx*m.pow(rx/rxgeoxx,2.) - ann*(Nx[0]**2)

def dN31dt(Nx, t, capt, anncs, rxvect, time, kbec, annbec):
    """ dN/dt = capt - ann*(N(t)**2) """
    if t<time[kbec]:
        # obtains rx at t by interpolating:
        rx = np.interp(t,time,rxvect)
        # recalcs ann calling funtion for rx
        ann = F_ann(anncs,rx)
    else:
        ann = annbec
    return capt - ann*(Nx[0]**2)

def dN32dt(Nx, t, capt, anncs, rxvect, time, selfcapt, Nx21geoxx, kgeoxx, kbec, annbec):
     """ dN/dt = capt + selfCapt*N(tgeo)*(rx(t)2/rxgeo2) - ann*(N(t)**2) """
     # obtains rx at t by interpolating:
     rx = np.interp(t,time,rxvect)
     if t<time[kbec]:
         # recalcs ann calling funtion for rx
         ann = F_ann(anncs,rx)
     else:
         ann = annbec
     # obtains rx at tgeoxx from vector:
     rxgeoxx = rxvect[kgeoxx]
     return capt + selfcapt*Nx21geoxx*m.pow(rx/rxgeoxx,2.) - ann*(Nx[0]**2)

def F_checkNxanalytical(i, t, dm, rxth, Nx_C_A):
    """ Calc Nx_C_A(t) from analytical solution """
    Nx_C_AA = [ m.sqrt(dm.Capt_v[i]/(0.5*F_ann(dm.anncs,rxth)))*m.tanh(m.sqrt(dm.Capt_v[i]*(0.5*F_ann(dm.anncs,rxth)))*t) for t in t.time]
    for it1,it2 in zip(Nx_C_A,Nx_C_AA):
        print "%.2e %.2e" % (it1,it2)
#    --> the test shows that both solution are identical

def d2rdt2(Y, t, Massx, nada):
    """ the 2nd order d2r/dt2=-GM/r2 is divided in Y = [Y[0], Y[1]] = [u, w], where r(t)=u and r'(t)=w, so u'=w and w'=-GM/u2 """
    uprime = Y[1]
    wprime = -Grav*Massx/m.pow(Y[0],2.)
    Yprime = np.array([uprime,wprime])
    return Yprime

def F_ann(anncs, rx):
    """ calculation of DM annihilation rate """
    try:
        return anncs / ((4./3.)*m.pi*(rx**3.))
    except ZeroDivisionError:
        print "ZERO DIVISION, rx=",rx

