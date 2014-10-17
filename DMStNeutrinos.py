#!/usr/bin/python
import numpy as np
from scipy.integrate import odeint
import math as m
from lib.DMcapture import F_capture
from lib.functions import *
from lib.unitconv import *
from lib.classes import *
from lib.constants import *

def calcALL(sigx, anncs):

    print " ********************************** "
    print "               INPUT " 
    print " ********************************** "
    # Input and initialize dark matter:
    dm = cl_dm(sigx, anncs)
    dm.prints()

    # Input star:
    st = cl_star('WD')  # 'NS' or 'WD' for Neutron Star or White Dwarf
    st.prints()

    print " ********************************** "
    print "            CALCULATIONS "
    print " ********************************** "
    for i in xrange(len(dm.mx_v)):
	  print " ********************************** "
	  print "          mx = %.2e " % dm.mx_v[i]
	  print " ********************************** "

          # Initializes lists that evolve with time:
          t = cl_timevars()

	  # check bullet cluster constraints:
          F_checksigxxconstraints(i, dm)

          # Characteristic time scales
          t1, tth = F_calctimescales(i, st, dm)

          # Characteristic radius of dark matter distribution, sets t.rx:
          rxco, rxth = F_calcDMradius(i, t, st, dm, t1, tth)

	  # Characteristic DM Numbers:
          Nsg,Nsgvect,Nfe,Nde,Ngeoxx,Ngeoxxvect,Nbec,Nchadeg = F_calcDMnumbers(i, t, st, dm, rxth)

          # Capture rates calculation:
          F_capture(i, st, dm)

	  # Calc Nx(t) numerically (solves dN/dt = capt + selfCapt*N(t) - ann*(N(t)**2)):
          t.Nx, Nx_C_A, Nx_C_sC_A = F_calcNxnum(i, t, dm)

          # Check if numerical solution agrees with analytical:
#          F_checkNxanalytical(i, t, dm, rxth, Nx_C_A)

	  # Check if self-collapse is produced, if geometrical limit in selfcapture is reached, if degenerate core and Bose-Einstein condensate are formed:
          kgeoxx, kbec  = F_checkSelfColl_Deg_geoxx_Bec_eq(i, t, dm, t.Nx, Nsgvect, Nde, Nfe, Ngeoxxvect, Nbec)

	  # Time when self-collapse starts and time of BEC formation:
          tNsg, tNbec = F_calctimeSelfColl_Bec(t, dm, t.Nx, Nsg, Nbec, kbec)

	  # In the special case when BEC is formed, calculates rxbec and annbec, and updates t.rx :
	  if dm.BEC: rxbec, annbec = F_calcBECcharact(i, t, dm, kbec)

	  # Recalcs Nx_C_sC_A(t) if selfcapture geoxx is reached (NOT CORRECT)
#          if dm.SELFCAPT:
#              Nx_C_sC_A = F_recalcNx_C_sC_A_with_geoxx(i, t, dm, Nx_C_sC_A, kgeoxx )
#              t.Nx = Nx_C_sC_A

          # Updates kbec and t.rx with rxbec (may have changed with new Nx_C_sC_A)
          if dm.BEC and dm.SELFCAPT: kbec = F_recalckbec(i, t, dm, Nx_C_sC_A, Nbec, rxbec)

	  # Calcs Nx_C_A_BEC(t) (capture and annihilation, with BEC after t>tNbec):
          if dm.BEC and not dm.SELFCAPT:
              Nx_C_A_BEC = F_calcNx_C_A_BEC(i, t, dm, kbec, annbec)
              t.Nx = Nx_C_A_BEC

          # Calcs Nx_C_sC_A_BEC(t) (capture, selfcapture and annihilation, with BEC after for t>tNbec):
          if dm.BEC and dm.SELFCAPT:
              Nx_C_sC_A_BEC = F_calcNx_C_sC_A_BEC(i, t, dm, Nx_C_sC_A, kgeoxx, kbec, annbec)
              t.Nx = Nx_C_sC_A_BEC

          # Recheck if self-collapse is produced with new Nx(t):
          if not (not dm.SELFCAPT and not dm.BEC): F_checkSelfColl_Deg(i, t, dm, t.Nx, Nsgvect, Nde, Nfe)

          # Calcs DM annihilation rate before collapse and checks if equilibrium capt=ann is reached:
          F_calcAnnrate_eq(i, t, dm)

	  # Computes the collapse of the DM cloud:
          timecoll, rxcoll ,vxcoll, Nxcoll, annrate, tff, kend = F_computeCollapse(i, t, dm, t.Nx)

          # Prints time evolution of parameters on screen:
          F_printallevol(i, t, dm, Nsgvect, kend, tff, timecoll, rxcoll, vxcoll, Nxcoll, annrate)

 	  # Computes emitted neutrino fluxes:   
          F_computeNeutFluxes(i, st, dm)
          
    # OUTPUT:
    dm.printOutput()

def main():
    for sigx in [1.e-30, 1.e-32, 1.e-34, 1.e-36, 1.e-38, 1.e-40, 1.e-42, 1.e-44, 1.e-48, 1.e-50 ]:
        for anncs in [1.e-26, 1.e-30, 1.e-35, 1.e-40, 1.e-45, 1.e-50, 1.e-55, 1.e-60]:
            calcALL(sigx, anncs)
  
if __name__ == "__main__":
    main()

