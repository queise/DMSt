import numpy as np
from scipy.integrate import quad,romberg
import math as m
from lib.unitconv import *
from lib.constants import *

#####################
#    Main function
#####################

def F_capture(i, st, dm):
    print "-- Capture:"
    capt = F_captBert(dm.rhox, dm.mx_v[i], dm.vx, dm.sigx, st.Ms, st.Rs, st.An)      # 1/s
#   capt = F_captGould(dm.rhox,dm.mx_v[i],dm.vx,dm.sigx_m,vs,st.Ms,st.Rs,st.mn,st.An) # (GeV/cm3,GeV,m/s,m2,m/s,kg,m,kg, )
    print "   capt=",capt #,"( captG=",captG," captB=",captB," )"
    # self-Capture calculation (see assumptions in function):
    selfcapt = 0.
    if dm.SELFCAPT:
#       selfcapt = F_selfCapt(dm.rhox,dm.mx_v[i],dm.sigxx_m,dm.vx,vs,st.Ms,st.Rs)  # errors probably in units
        selfcapt_yr = 1.06e-3 * (dm.sigxx/1.e-24) / dm.mx_v[i] * dm.rhox    # yr-1
        selfcapt    = c_s2yr(selfcapt_yr)     # s-1
        print "   selfcapt=",selfcapt,"s-1 or ",selfcapt_yr," yr"
    dm.Capt_v[i]     = capt
    dm.selfCapt_v[i] = selfcapt

#####################
#    Gould (87)  (adapted from Moskalenko)
#####################
def F_captGould(rhox, mx, vx, sigx, vs, Ms, Rs, mn, An):
# should be   GeV/cm3,GeV,m/s,m2,m/s,kg,m,kg, )
    mxkg   = c_GeV2kg(mx)
    rhoxGeVm3 = rhox*1.e+6
    nu     = 3.*m.pow(vs,2)/2./m.pow(vx,2)  # nu is probably wrong in here and in moskalenko. correct would be sqrt(nu)
    print '*** Remember that nu in F_captGould may be wrong. ***'
    sys.exit() 
    mu     = mxkg / mn
    print "nu:",nu,"mu:",mu
    muminus= (mu - 1.)/2.
    sigeff = min( sigx*m.pow(An,4.)*Ms/mn , m.pi*m.pow(Rs,2.) )
    print "sigeff=",sigeff
    tuple = rhoxGeVm3,mx,vx,vs,Ms,Rs,nu,mu,muminus
#    integr = quad(F_r2dCdV,0.,Rs,args=(tuple)) # scipy.integrate.quad(func, a, b, args=(), full_output=0, epsabs=1.49e-08, epsrel=1.49e-08, limit=50, points=None, weight=None, wvar=None, wopts=None, maxp1=50, limlst=50)	
    integr1= quad(F_r2dCdV,0.,Rs/3.,args=(tuple))
    integr2= quad(F_r2dCdV,Rs/3.,2*Rs/3.,args=(tuple))
    integr3= quad(F_r2dCdV,2*Rs/3.,Rs,args=(tuple))
#    print "0 - ",Rs/100.," - ",Rs/10.,"- ",Rs,": ",integr1[0],"+",integr2[0],"+",integr3[0]
#    integr = romberg(F_r2dCdV,0.,Rs,args=(tuple))
#    print "integr:",integr[0]
#    return 4.*m.pi*integr[0]
    return 4.*m.pi*sigeff*(integr1[0]+integr2[0]+integr3[0])

def F_r2dCdV(rad, rhoxGeVm3, mx, vx, vs, Ms, Rs, nu, mu, muminus):
    vesc   = F_vesc(rad,Ms,Rs)
    A      = m.sqrt( 3*m.pow(vesc,2.)*mu / 2./m.pow(vx,2)/m.pow(muminus,2.) )
    Aplus  = A + nu
    Aminus = A -nu
#    print "A,Aplus,Aminus:",A,Aplus,Aminus
    term1 = m.sqrt(6./m.pi)/(4./3.*m.pi*m.pow(Rs,3.))* rhoxGeVm3/mx \
	        * m.pow(vesc/vx,2.)*vx/2./nu/m.pow(A,2.)
    term2 = ((Aplus*Aminus)-0.5)*(F_chi(-nu,nu)-F_chi(Aminus,Aplus))
    term3 = (0.5*Aplus*m.exp(-m.pow(Aminus,2.))) - (0.5*Aminus*m.exp(-m.pow(Aplus,2.))) - \
	        (nu*m.exp(-m.pow(nu,2.)))
#    print "rad,vesc,term1,term2,term3:",rad,vesc,term1,term2,term3
    return pow(rad,2.) * term1 * (term2+term3)

def F_vesc(rad, Ms, Rs):
    return m.sqrt( Grav*Ms/Rs*(3.-m.pow(rad/Rs,2.)) )

def F_chi(a, b):
    return m.sqrt(m.pi/2) * (m.erf(b)-m.erf(a))


##################################
#    Bertone & Fairbairn (2008)
##################################
def F_captBert(rhox, mx, vx, sigx, Ms, Rs, An):
    vesc = m.sqrt(2*Grav*Ms/Rs)     # m/s
#    print "vesc=",vesc
    siggeo = m.pi * m.pow(Rs*1.e+2,2) # cm2
    sigDMp = sigx*Ms/mp*m.pow(An,3.) # cm2
    sigeff = min( sigDMp, siggeo)  # cm2
    if sigeff == siggeo:
        print "   geolimit active in capture. Max sigx is:",siggeo/(Ms/mp*m.pow(An,3.))
#        print "n=",sigx*Ms/mp*m.pow(An,3.)," geo=",m.pi * m.pow(Rs*1.e+2,2)," sigeff=",sigeff, "sigxmax=",m.pi*m.pow(Rs*1.e+2,2)/(Ms/mp*m.pow(An,3.))
    return ( (8./(3.*m.pi))**0.5 ) *rhox*(vx*100.)/mx * ( 3.*(vesc**2.)/(2.*(vx**2.)) ) * sigeff

def F_captBert2(rhox, mx, vx, sigx, Ms, Rs, An):	# INCLUDES GRAV ENHANCEMENT
    vesc = m.sqrt(2*Grav*Ms/Rs)     # m/s
    print "F_captBert2           vesc=",vesc
    siggeo = m.pi * m.pow(Rs*1.e+2,2) # cm2
    sigDMp = sigx*Ms/mp*m.pow(An,3.) # cm2
    sigeff = min( sigDMp, siggeo)  # cm2
    if sigeff == siggeo:
        print "   geolimit active in capture. Max sigx is:",siggeo/(Ms/mp*m.pow(An,3.))
    sigeff2 = sigeff * ( 1. + m.pow(vesc/vx,2.) )
    return ( (8./(3.*m.pi))**0.5 ) *rhox*(vx*100.)/mx * ( 3.*(vesc**2.)/(2.*(vx**2.)) ) * sigeff2


##################################
#    adapted from Guver et al. (2012, arXiv:1201.2400) originally from Zentner (2009) 0907.3448,
#    but with account of diff vs and vx
##################################
def F_selfCapt(rhox, mx, sigxx, vx, vs, Ms, Rs):
#	       GeV/cm3,GeV,m2,m/s,m/s,kg,m
    avpot 	= 1	# see eq.25, =1 is conservative 
    nu	= m.sqrt( 3.*m.pow(vs,2)/2./m.pow(vx,2) )
    erfnn 	= m.erf(nu)/nu
    print "erfnn=",erfnn
    vesc = m.sqrt(2.*Grav*Ms/Rs)	# m/s
    print "vesc=",vesc
    print "(1./(1.-(m.pow(vesc,2))))",(1./(1.-(m.pow(vesc,2))))
    return m.sqrt(3./2.)*(rhox/mx)*sigxx*(vesc/vx)*avpot*erfnn*(1./(1.-(m.pow(vesc,2))))

#####################
#    Kouvaris
#####################
def F_captKouv(rhox, mx):
    return 1.25e+24 * rhox * (100./mx)

#####
# Sun
#####
def F_captSun(rhox, mx, vx, sigsi):
    CwptSI = 1.24e+20 * (rhox/0.3) * ((270./vx)**3) * 5.4*(sigsi/1.e-42) * ((100./mx)**2)
    return CwptSI
