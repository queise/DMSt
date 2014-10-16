import sys
import numpy as np
import math as m
from scipy.integrate import odeint
from lib.constants import *
from lib.unitconv import *
from lib.NxFunctions import *

def F_checksigxxconstraints(i, dm):
    if((dm.sigxx/1.e-24/dm.mx_v[i])>2.):
        print """   *******************************

                 sigxx not allowed by Bullet cluster constraints

                       *********************************"""
        sys.exit() 

#def F_rxco(t,mxkg,sigx):
#	if t==0:
#	  t=1.e-10
#	return (mxkg/mp) * (2.8e+10/(sigx/1.e-55)/c_s2yr(t))	# cm

def F_rxco2(t, t1, mr, nb, sigx_m, Rs, mxkg, pF):
    """ radius of DM orbits """
    if (t-t1)<=0.:
        tt1=1.e-10
    else:
        tt1=t-t1
    term = 1.+ ( 8.*m.pi*m.sqrt(2.)*m.pow(mr,3.)*m.pow(nb,2.)*sigx_m*Grav*m.pow(Rs,2.)*tt1/3./mxkg/pF )
#    print "%.2e   %.2e   %.2e" % (t,tt1,term)
    return Rs / m.sqrt(term) 

def F_rxth(mx, rhoc, Temp):
    """ radius of DM distribution after thermalization with stellar matter """
    return 64. * m.pow(Temp/1.e+5,0.5) * m.pow(1.e+14/rhoc,0.5) * m.pow(100./mx,0.5) # cm

def F_calctimescales(i, st, dm):
    """ Calculates characteristic time scales """
    mr = st.mn*dm.mxkg_v[i]/(st.mn+dm.mxkg_v[i])        # reduced mass, kg
    # containment time (orbits within Rs)
    t1   = c_yr2s( 2.7e-2 * m.pow(dm.mxkg_v[i]/st.mn,1.5) / (dm.sigx/1.e-55) )       # s          
    # therm time
    tth1 = c_yr2s( 2.5e+5 * m.pow(dm.mxkg_v[i]/st.mn,2.) * m.pow(st.mn/mr,3.) / (dm.sigx/1.e-55) )    # s
#   tth2 = m.pow(dm.mxkg_v[i]/st.mn,2.) * pF / 6. / m.sqrt(2.) / st.Temp / st.nb / dm.sigx * m.pow(st.mn/mr,3.) 
    tth2 = m.pow(dm.mxkg_v[i],2.)*st.mn*pF /4./m.sqrt(2.)/(st.nb*1.e+6)/dm.sigx_m/m.pow(st.mn,3.)/st.Eth
    tth  = tth2
    print "-- Time scales: t1=%.2e , tth=%.2e, tth1=%.2e, tth2=%.2e" % (t1,tth,tth1,tth2)
    return t1, tth 

def F_calcDMradius(i, t, st, dm, t1, tth):
    """ Calculates the DM sphere radius before and after thermalization """
    mr = st.mn*dm.mxkg_v[i]/(st.mn+dm.mxkg_v[i])        # reduced mass, kg
    # before thermalization (cooling), rx changes with time:
    rxco  = np.array([ F_rxco2(tim,t1,mr,(st.nb*1.e+6),dm.sigx_m,st.Rs,dm.mxkg_v[i],pF) for tim in t.time ]) # cm
    print "-- Radius: rxco at t1 = ",F_rxco2(t1+0.1,t1,mr,(st.nb*1.e+6),dm.sigx_m,st.Rs,dm.mxkg_v[i],pF)
    # after thermalization:
    rxth1 = F_rxth(dm.mx_v[i],st.rhoc,st.Temp)  # cm (formula)
    rxth2 = np.interp(tth,t.time,rxco)    	# cm (rxco(tth))
    rxth  = rxth1
    print "   rxth=%.2e , rxth1=%.2e , rxth2=%.2e" % (rxth,rxth1,rxth2)
    for k in xrange(len(t.time)):
        if t.time[k]<t1:
            t.rxtag[k]  = 'Rs  '
            t.rx[k] = st.Rs*1.e+2
        elif t.time[k]<tth:
            t.rxtag[k]  = 'rxco'
            t.rx[k] = rxco[k]
        elif t.time[k]>=tth:
            t.rxtag[k]  = 'rxth'
            t.rx[k] = rxth
    return rxco, rxth

def F_calcDMnumbers(i, t, st, dm, rxth):
    Nsg    = (4./3.)*m.pi*m.pow(rxth,3.)*st.rhoc/dm.mxg_v[i]      # DM selg-gravitates for N>Nsg
    Nsgvect= np.array([(4./3.)*m.pi*m.pow(it,3.)*st.rhoc/dm.mxg_v[i] for it in t.rx])
#   print "-- Nsgvect=",Nsgvect
    Nfe   = 3.e+37 * m.pow(st.Temp/1.e+5,3.) * m.pow(1.e+14/st.rhoc,1.5)
    Nde   = 2.5e+39 * m.pow(st.Temp/1.e+5,21/10) * m.pow(1.e+14/st.rhoc,9/10) * m.pow(100/dm.mx_v[i],1.5)
    Ngeoxx= m.pi*m.pow(rxth,2.)/dm.sigxx          # geolimin of self-capture reached for N>Ngeoxx
    Ngeoxxvect = np.array([m.pi*m.pow(it,2.)/dm.sigxx for it in t.rx])
    Nbec  = 2.7e+36*(m.pow(st.Temp/1.e+5,3.))     # bosons form a BEC for N>Nbec
    Nchadeg  = 1.e+51 *(m.pow(100./dm.mx_v[i],3.))
    print "-- N: Nsg(fix rx)=%.2e Nfe=%.2e Nde=%.2e Ngeoxx=%.2e Nbec=%.2e Nchadeg=%.2e" % (Nsg,Nfe,Nde,Ngeoxx,Nbec,Nchadeg)
#   print "      Ngeoxxvect=",Ngeoxxvect
    return Nsg,Nsgvect,Nfe,Nde,Ngeoxx,Ngeoxxvect,Nbec,Nchadeg

def F_checkSelfColl_Deg_geoxx_Bec_eq(i, t, dm, Nx, Nsgvect, Nde, Nfe, Ngeoxxvect, Nbec):
    """ Checks if DM self-collapse is produced and if degenerate star forms (through F_checkSelfColl_Deg),
        and also checks if Bose-Einstein condensate forms, and if equilibrium between DM capture and annihilation is reached. """
    kgeoxx = 0 ; kbec = 0
    kscoll = F_checkSelfColl_Deg(i, t, dm, Nx, Nsgvect, Nde, Nfe)
    if dm.SELFCAPT:
        for k in xrange(len(t.time)):
            if Nx[k]>Ngeoxxvect[k]:
                print "Geometrical limit in self-capture was reached at t=%.2e." % t.time[k]
                dm.notes_t[i].append(4)
                kgeoxx = k
                break
    if dm.BOSON:
        for k in xrange(len(t.time)):
            if Nx[k]>Nbec:
                print "Nbec(=%e) was reached at t=%.2e, check if BEC was correctly taken into account" % (Nbec,t.time[k])
                dm.notes_t[i].append(5) 
                kbec = k
                dm.BEC = True
                break
        if Nx[len(t.time)-1]<Nbec:
            dm.BEC = False
            print "DM is boson but BEC is not formed (Nxmax=%.2e and Nbec=%.2e)" % (Nx[len(t.time)-1],Nbec)
        if dm.BEC and dm.SELFCAPT and kbec <= kgeoxx:
            print " BEC occurs before geoxx limit, check if implementation is correct in this case"
            dm.notes_t[i].append(6)
    return kscoll, kgeoxx, kbec

def F_calctimeSelfColl_Bec(t, dm, Nx, Nsg, Nbec, kbec):
    """ More characteristic time scales.
        Selfcollapse time, interpolates to find the time at which Nx = Nsg
        if Nsg is not within Nx(t) the np.interp wrongly gives the last value, but error note 1 was already set to detect that. """
    tNsg = np.interp(Nsg,Nx,t.time) # fixed rx case, not used
    if dm.BEC:
        tNbec= np.interp(Nbec,Nx,t.time)
        print " **************** "
        print " **************** "
        print " **************** "
        print " **************** CHECK THIS:"
        print " tNbec=",tNbec,"  t.time[kbec]=",t.time[kbec]
    else: tNbec = None
    return tNsg, tNbec

def F_calcBECcharact(i, t, dm, kbec):
    assert dm.BEC == True, "Error, function F_calcBECcharact should only be called if BEC is formed"
    rxbec = 1.5e-5 * m.pow(100./dm.mx_v[i],0.5)   # cm
    annbec= 0.5 * F_ann(dm.anncs,rxbec)
    print "   rxbec=%.2e , annbec=%2.e " % (rxbec,annbec)
    # after BEC is reached, rx is reset to rxbec
    for k in xrange(len(t.time)):
        if k>=kbec:
            t.rxtag[k]  = 'Rbec'
            t.rx[k]     = rxbec
    return rxbec, annbec

def F_recalcNx_C_sC_A_with_geoxx(i, t, dm, Nx_C_sC_A, kgeoxx ):
    """ Update case 2: capture, selfcapture and annihilation, with rx(t) and with geoxx limit.
        The Nx_C_sC_A_geoxx SEEMS NOT CORRECT """
    assert dm.SELFCAPT== True, "Error, function F_recalcNx_C_sC_A_with_geoxx should only be called if SELFCAPT"
    tuple22 = dm.Capt_v[i], dm.anncs, t.rx, t.time, dm.selfCapt_v[i], Nx_C_sC_A[kgeoxx], kgeoxx         # extra arguments for function dNdt_C_sC_A_geo$
    Nx_C_sC_A_geoxx = np.array([item[0] for item in odeint(dNdt_C_sC_A_geoxx,dm.Nxinit,t.time,args=(tuple22))])
    for k in xrange(len(t.time)):
#        print 'Nx_C_sC_A_geoxx[k]=',Nx_C_sC_A_geoxx[k]
        if k>=kgeoxx:
            Nx_C_sC_A[k] = Nx_C_sC_A_geoxx[k]       # changes Nx_C_sC_A after geoxx is reached
    print "   ... Nx_C_sC_A(t) re-calculated to account for geo limit in selfcapture"
    return Nx_C_sC_A

def F_calcNx_C_A_BEC(i, t, dm, kbec, annbec):
    """ Update case 1: capture and annihilation, with BEC after t>tNbec """
    assert (dm.BOSON and dm.BEC and not dm.SELFCAPT), "Error, function F_calcNx_C_A_BEC should not have been called"
    tuple31 = dm.Capt_v[i], dm.anncs, t.rx, t.time, kbec, annbec         # extra arguments for function dN3dt
    Nx_C_A_BEC = np.array([item[0] for item in odeint(dN31dt,dm.Nxinit,t.time,args=(tuple31))])
    print "   ... Nx31(t) calculated"
    return Nx_C_A_BEC

def F_recalckbec(i, t, dm, Nx_C_sC_A, Nbec, rxbec):
    assert (dm.BOSON and dm.BEC and dm.SELFCAPT), "Error, function F_recalckbec should not have been called"
    # kbec is recalculated with N22
    print "... kbec recalculated with new Nx_C_sC_A(t)"
    if 4 in dm.notes_t[i]: dm.notes_t[i].remove(4)
    for k in xrange(len(t.time)):
#        print 'Nx_C_sC_A[k]=',Nx_C_sC_A[k],' Nbec=',Nbec
        if Nx_C_sC_A[k]>Nbec:
            print "Nbec(=%e) was reached at t=%.2e and should be taken into account" % (Nbec,t.time[k])
            if 4 not in dm.notes_t[i]: dm.notes_t[i].append(4) #dm.notes_v[i]=dm.notes_v[i]+' 6'
            kbec = k
            break
    # rvect is reset acordingly  
    for k in xrange(len(t.time)):
        if k>=kbec:
            t.rxtag[k]  = 'Rbec'
            t.rx[k] = rxbec
    return kbec

def F_calcNx_C_sC_A_BEC(i, t, dm, Nx_C_sC_A, kgeoxx, kbec, annbec):
    """ capture, selfcapture and annihilation, with BEC after for t>tNbec """
    assert (dm.BOSON and dm.BEC and dm.SELFCAPT), "Error, function F_calcNx_C_sC_A_BEC should not have been called"
    tuple32 = dm.Capt_v[i], dm.anncs, t.rx, t.time, dm.selfCapt_v[i], Nx_C_sC_A[kgeoxx], kgeoxx, kbec, annbec         # extra arguments for function d$
#    print "tuple32= capt, dm.anncs, t.rx, t.time, selfcapt, Nx_C_sC_A[kgeoxx], kgeoxx, kbec, annbec"
#    print "tuple32=",tuple32
    Nx_C_sC_A_BECtemp = np.array([item[0] for item in odeint(dN32dt,dm.Nxinit,t.time,args=(tuple32))])
    Nx_C_sC_A_BEC     = np.zeros(len(t.time))
    for k in xrange(len(t.time)):
        if k<kbec:
            Nx_C_sC_A_BEC[k] = Nx_C_sC_A[k]       # before bec is reached
        else:
            Nx_C_sC_A_BEC[k] = Nx_C_sC_A_BECtemp[k]      # after bec
    print "   ... Nx_C_sC_A_BEC(t) calculated"
    return Nx_C_sC_A_BEC


def F_checkSelfColl_Deg(i, t, dm, Nx, Nsgvect, Nde, Nfe):
    print "... checking conditions for self-collapse and deg hold:"
    kscoll = 0
    if 1 in dm.notes_t[i]: dm.notes_t[i].remove(1)
    for k in xrange(len(t.time)):
        if k>=1 and (Nx[k] >= 1.e+100 or Nx[k] <= 1.e-100):
            print "-- Error at k=%i, t=%.2e " % (k,t.time[k])
            dm.notes_t[i].append(0)
            dm.SELFCOLL[i]  = False
            kscoll = k # well, is not k of selfcollapse, but k of error, but for printing output it work the same
            break
        if Nx[k] >= Nsgvect[k]:
            kscoll = k
            print "-- DM reaches self gravitation (self-collapse) at t=%.2e" % t.time[k]
            print "t.time[kscoll]=%.2e, Nx[kscoll]=%.2e "% (t.time[kscoll],Nx[kscoll])
            dm.SELFCOLL[i]  = True
            break
        elif k == (len(t.time)-1):
            print "-- No DM self-collapse (Nsg not reached within the age of the Universe)."
            dm.notes_t[i].append(1) 
            kscoll = (len(t.time)-1)
            dm.SELFCOLL[i]  = False
    if 3 in dm.notes_t[i]: dm.notes_t[i].remove(3)
    for k in xrange(len(t.time)):
        if (Nx[k] >= Nde and Nx[k] >= Nfe and Nx[k] < Nx[kscoll]):
            print "-- Attention: degenerate Dark Star formed at t=%.2e (before self-collapse)" % t.time[k]
            if 3 not in dm.notes_t[i]: dm.notes_t[i].append(3)
            dm.DEG = True
            break

    return kscoll

def F_computeCollapse(i, t, dm, Nx, kscoll):
    if dm.SELFCOLL[i]:
        print "-- calculating DM self-collapse"
        dm.tNsg_v[i] = t.time[kscoll]
        dm.Nsg_v[i]  = Nx[kscoll]
        rxsg         = t.rx[kscoll]         # cm
        vxsg         = (t.rx[kscoll]-t.rx[kscoll-1])/(t.time[kscoll]-t.time[kscoll-1])      # cm/s
        print "   the self-collapse was produced at t=%.2e when Nx=%2.e" % (dm.tNsg_v[i],dm.Nsg_v[i])
        # free-fall time:
        Massx  = c_GeV2kg(dm.Nsg_v[i]*dm.mx_v[i])
        densff = Massx/4.*3./m.pi/m.pow(rxsg*1.e-2,3.)
        tff    = m.sqrt(3.*m.pi/32./Grav/densff)
        # timecoll list should be of smaller than tff
        timecoll = [ (tff*10.**(k/10.)) for k in range(-20,10)]  # s
#       timecoll = [ (10.**(k/10.)) for k in range(-40,0)]
#        timecoll = [ (0.99*tff)+(k*tff/500000.) for k in xrange(1,5000000)]
        # Solving 2nd order ordinary diff eq: d2r/dt2 = -GM/r2, converting it to 1rst order: rxvx = [r(t) , v(t)] see function d2rdt2
        rxinit  = np.array([rxsg*1.e-2, vxsg*1.e-2])  # initial r and v
        nada    = 1.
        tuplerx = Massx, nada                 # extra arguments for function d2rdt2
        rxvx    = odeint(d2rdt2,rxinit,timecoll,args=(tuplerx))
        rxcoll  = np.array([item[0]*1.e+2 for item in rxvx])        # cm
        vxcoll  = np.array([item[1] for item in rxvx])              # m/s
        # Checking when rx gets to zero:
        kend = len(timecoll)-1
        for k in xrange(len(timecoll)):
            if rxcoll[k]==0. or rxcoll[k]>1.e+15 or rxcoll[k]<1.e-50:
                kend=k-1
                break
        # re-size vectors to avoid the rx zeros:
        rxcoll   = np.resize(rxcoll,kend+1)
        timecoll = np.resize(timecoll,kend+1)
        vxcoll   = np.resize(vxcoll,kend+1)
        # Evolution of Nx(t) during collapse, only capture and annihilation, with rx(t), neglecting selfcapture
        tuple1  = dm.Capt_v[i], dm.anncs, rxcoll, timecoll
        Nxcoll  = np.array([item[0] for item in odeint(dNdt_C_A,dm.Nsg_v[i],timecoll,args=(tuple1))])
        # Energy from annihilated particles:
        Volx    = np.zeros(len(timecoll))
        annrate = np.zeros(len(timecoll))
        for k in xrange(0,kend+1):
            Volx[k]    = 4./3.*m.pi*m.pow(rxcoll[k],3.)
            try:
                annrate[k] = dm.anncs/Volx[k] * m.pow(Nxcoll[k],2.)
            except OverflowError:
                annrate[k] = 9.e+99
    else:
        timecoll = rxcoll = vxcoll = Nxcoll = annrate = tff = kend = 0

    return timecoll, rxcoll ,vxcoll, Nxcoll, annrate, tff, kend

def F_calcAnnrate_eq(i, t, dm):
    eps = 0.01
    ktemp = 0
    for k in xrange(len(t.time)):
        try:
            t.annirate[k] = F_ann(dm.anncs,t.rx[k]) * m.pow(t.Nx[k],2.)
        except OverflowError:
            t.annirate[k] = 9.e+99
        if ktemp==0 and (t.annirate[k]*(1-eps) < dm.Capt_v[i] < t.annirate[k]*(1+eps)):
            print "-- Equilibrium between capture and annihilation rates reached at t=%.2e" % t.time[k]
            dm.FluxEmit_v[i] = t.annirate[k] * dm.mx_v[i] * c_yr2s(1.)
            ktemp=1
        if k == (len(t.time)-1):
            if (t.annirate[k]*(1-eps) < dm.Capt_v[i] < t.annirate[k]*(1+eps)):
                print "   Equilibrium lasted until the end of time."
                dm.notes_t[i].append('2a') 
            elif ktemp==1:
                print "   Equilibrium disappeared (probably because of transition from rxco to rth, that should be rechecked)"
                dm.notes_t[i].append('2b')

def F_printallevol(i, t, dm, Nsgvect, kscoll, kend, tff, timecoll, rxcoll, vxcoll, Nxcoll, annrate):
    print" time(s)   time(yr)     Nx        Nsg      rx(type)          ann   "
    for k in xrange(0,kscoll+1):
        print "%.2e   %.2e   %.2e   %.2e   %.2e (%s)   %.2e   " % (t.time[k], c_s2yr(t.time[k]), t.Nx[k],
                                                                Nsgvect[k], t.rx[k], t.rxtag[k], t.annirate[k])
    if dm.SELFCOLL[i]:
        print" now DM self collapses in a timescale < %.2e s" % tff
        print " t         rx(cm)       vx/c       Nx      annrate(part/s)"
        for k in xrange(0,kend+1):
            print "%.2e   %.2e   %.2e   %.2e   %.2e" % (timecoll[k],rxcoll[k],vxcoll[k]/clight,Nxcoll[k],annrate[k])

def F_computeNeutFluxes(i, st, dm):
    # Flux of energy emitted:
    # if not SELFCOLL, dm.FluxEmit_v[i] was already calculated
    if dm.SELFCOLL[i]:
        dm.FluxEmit_v[i]  = dm.Nsg_v[i] * dm.mx_v[i] / c_s2yr(dm.tNsg_v[i])               # GeV/yr
    # Flux of neutrinos at Earth:
    dm.FluxEarth_v[i] = dm.FluxEmit_v[i] /4./m.pi/m.pow(c_pc2m(st.dist),2.)               # GeV/m2/yr

