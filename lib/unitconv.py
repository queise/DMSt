#!/usr/bin/python
import math as m

# MASS #########################

def c_GeV2g(mass):
        return mass * 1.782662e-24

def c_GeV2kg(mass):
	return mass * 1.782662e-27

def c_GeV2Msun(mass):
	return c_kg2Msun(c_GeV2kg(mass))

def c_Msun2kg(mass):
	return mass * 1.98855e+30

def c_Msun2GeV(mass):
	return mass * 1.98855e+30 / 1.782662e-27	

def c_kg2Msun(mass):
	return mass / 1.98855e+30

# DISTANCE #####################

def c_Rsun2m(dist):
	return dist * 6.955e+8
def c_m2Rsun(dist):
        return dist / 6.955e+8
def c_pc2m(dist):
	return dist * 3.08567758066631e+16
def c_cm2pc(dist):
	return dist / 3.08567758066631e+18

def c_ly2pc(dist):
	return dist * 0.306601

# TIME #########################

def c_s2yr(time):
	return time/3600./24./365.

def c_yr2s(time):
	return time*3600.*24.*365.

# ENERGY #######################

def c_GeV2erg(ener):
	return ener * 1.602176565e-10 * 1.e+7	# GeV to J to erg

# OTHERS #######################

def c_ergs2Lsun(lum):
	return lum / 3.846e+33

def c_Msunpc32GeVcm3(rho):
	return c_Msun2GeV(rho)/m.pow(c_pc2m(1.),3)/1.e+6

def c_GeVcm32Msunpc3(rho):
	return c_GeV2Msun(rho)/m.pow(c_cm2pc(1.),3)
