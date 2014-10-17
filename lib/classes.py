import sys
import numpy as np
from lib.unitconv import *
from lib.constants import mp, kb

class cl_star():
    def __init__(self,typeSt):
        self.typeSt = typeSt
        if typeSt == 'NS':
            self.Ms = c_Msun2kg(1.44)
            self.Rs = 10600.     # m
            self.Temp = 1.e+5    # K
            self.densStars = 5.e-4       # NumNS/pc3
            self.An = 1
            self.dist = 200	# pc 
        elif typeSt == 'WD':
            self.Ms = c_Msun2kg(1.4)  #1.4-0.6 # from Msun to kg
            self.Rs = c_Rsun2m(0.0045) #0.0045-0.02 # from Rsun to m
            self.Temp = 5.e+6  # K
            self.densStars = 5.e-3       # NumWD/pc3
            self.An = 12
            self.dist = 200    # pc
        else:
            print "Type of star not recognized. Only 'NS' and 'WD' accepted. You tried: '%s'" % typeSt
            sys.exit()
        self.rhoc = self.Ms*1.e+3/ (4./3.*m.pi*m.pow(self.Rs*1.e+2,2.)) # g/cm3
        self.vs = 22000.    # m/s
        self.mn = self.An * mp    # kg
        self.nb = self.rhoc*1e-3 / self.mn # numb.bary / cm3
#        self.pF = 426.e+6 * 5.344286e-28 # from eV to kg m/s, Fermi momentum 
        self.Eth= (3./2.)*kb*self.Temp

    def prints(self):
        print "--- Stellar characteristics"
        print "    type of star: ",self.typeSt
        print "    Ms = %.2f Msun,  Rs = %.5f Rsun,  rhoc = %.2e g/cm3,  Temp = %.2e K,   An = %i,  vs = %.2e m/s "  % (c_kg2Msun(self.Ms),c_m2Rsun(self.Rs),self.rhoc,self.Temp,self.An,self.vs)

class cl_dm():
    def __init__(self, sigx, anncs):
        # inputs:
        self.rhox  	= 0.4     					# GeV/cm3
        self.mx_v 	= [ 10.**i for i in range(8,8+1)]  		# GeV
#        self.mx_v	= [ 10.**(2*i-2) for i in range(2,6+1)]# (1,3+1)]  	# GeV
        self.mxkg_v     = [ c_GeV2kg(mx) for mx in self.mx_v ]    				# kg
        self.mxg_v  	= [ c_GeV2g(mx) for mx in self.mx_v ]      				# g
        self.anncs 	= anncs # 1.e-50  						# cm3/s
        self.sigx  	= sigx    ; self.sigx_m  = self.sigx/1.e+4  	# m2
        self.sigxx 	= 1.e-30  ; self.sigxx_m = self.sigxx/1.e+4    	# m2
        self.vx    	= 27000.  					# m/s
        self.SELFCAPT 	= False
        self.BOSON    	= False
	self.DEG 	= False		# ?????????????????????
        if self.BOSON==False:
            self.BEC = False
	# initializing other variables, with different values for each DM mass:
        self.Nsg_v 	= np.zeros(len(self.mx_v))
	self.tNsg_v	= np.zeros(len(self.mx_v))
        self.notes_v	= ['' for k in range(0,len(self.mx_v))]
        self.notes_t    = [[] for k in range(0,len(self.mx_v))]
#        self.Flux_v	= np.zeros(len(self.mx_v))
	self.FluxEmit_v = np.zeros(len(self.mx_v))
        self.FluxEarth_v= np.zeros(len(self.mx_v))
        self.Capt_v	= np.zeros(len(self.mx_v))
        self.selfCapt_v = np.zeros(len(self.mx_v))
	self.SELFCOLL   = [False]*len(self.mx_v)
        self.kscoll     = np.zeros(len(self.mx_v),dtype=int)
        # more:
	self.Nxinit = 0.	# np.array([0.])                   # initial values N(t=0)
	
    def prints(self):
        print "--- Dark matter:"
        print "    rhox = %.2e GeV/cm3,  anncs = %.2e cm3/s,  sigx = %.2e cm2,  vx = %.2e m/s"   % (self.rhox,self.anncs,self.sigx,self.vx)
        print "    selfinteractions:",self.SELFCAPT," (sigxx = %.2e cm2)" % self.sigxx
        print "    bosonic DM:",self.BOSON

    def printOutput(self):
        # prepares last column for output file:
        codecol = []
        for i in xrange(len(self.mx_v)):
            if 0 in self.notes_t[i]:
                codecol.append(2)		# case 2: all cases which I am not confindent and should be rechecked
            elif '2b' in self.notes_t[i]:
                codecol.append(2)           # case 2: all cases which I am not confindent and should be rechecked
            elif 3 in self.notes_t[i]:
                codecol.append(2)           # case 2: all cases which I am not confindent and should be rechecked
            elif 4 in self.notes_t[i]:
                codecol.append(2)           # case 2: all cases which I am not confindent and should be rechecked
            elif 5 in self.notes_t[i]:
                codecol.append(2)           # case 2: all cases which I am not confindent and should be rechecked
            elif 6 in self.notes_t[i]:
                codecol.append(2)           # case 2: all cases which I am not confindent and should be rechecked
            elif 1 in self.notes_t[i]:
                codecol.append(0)		# case 0: no DM collapse, continuous neutrino flux
            elif self.SELFCOLL[i]:
                codecol.append(1)           # case 1: DM self-collapse, neutrino flashes

            f = open('DMStNeutrinoFlux.dat', 'a')
            print ""
            print "DM collapse and neutrino flash from self-annihilation burst or continuous neutrino emission:"
            print "sig_el_x   sig_ann_x  mx(GeV)  Nx(collapse) period(yr) Fluxemitt(GeV/yr) Flux@Earth(GeV/m2/yr) notes"
#            for it1,it2,it3,it4,it5 in zip(self.mx_v,self.tNsg_v,self.Nsg_v,self.FluxEarth_v,self.notes_t):
#                print " %.2e  %.2e    %.2e      %.2e              %s" % (it1,c_s2yr(it2),it3,it4,it5)
            for it1,it2,it3,it4,it5,it6,it7 in zip(self.mx_v, self.Nsg_v, self.tNsg_v, self.FluxEmit_v, self.FluxEarth_v, self.notes_t, codecol):
                print "%.2e   %.2e   %.2e   %.2e   %.2e   %.2e          %.2e              %s" % (self.sigx, self.anncs, it1, it2, c_s2yr(it3), it4, it5, it6)
                f.write("%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %i\n" % (self.sigx, self.anncs, it1, it2, c_s2yr(it3), it4, it5, it7) )
            f.close()

            print """\nNotes:
    0:  Error (check Nx or annrate in last time step)
    1:  No DM self-collapse (Nsg not reached within the age of the Universe)  
    2:  Equilibrium capture = annihilation rates reached, and (a) lasted or (b) disappeared (prob error in rx transition)"
    3:  Degenerate Dark Star formed (Nsg>Nde and >Nfe)
    4:  Geometrical limit in self-capture was reached: attention that Nx_C_sC_A_geoxx is not correct
    5:  Bose-Einstein condensate formed, check if BEC is correctly implemented
    6:  BEC occurs before geoxx limit, check if implementation is correct in this case
             """

class cl_timevars():
    def __init__(self):
        self.time     	= np.array([ (10.**(k/3.))-1. for k in range(0,54)])  # s
        
#        time  = np.linspace(0.,5.e+12,20)     # numpy.linspace(start, stop, num=50, endpoint=True, retstep=False)
        self.rx		= np.zeros(len(self.time))
        self.rxtag    	= ['yes' for k in range(0,len(self.time))]
        self.ann	= np.zeros(len(self.time))
        self.annirate 	= np.zeros(len(self.time))
	self.Nx		= np.zeros(len(self.time))
	self.Nsg	= np.zeros(len(self.time))
	self.Nx22     	= np.zeros(len(self.time))
        self.Nx3      	= np.zeros(len(self.time))
        self.Nx31     	= np.zeros(len(self.time))
        self.Nx32     	= np.zeros(len(self.time))
