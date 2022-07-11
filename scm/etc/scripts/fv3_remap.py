#!/usr/bin/env python

import numpy as np
import math

r3 = 1./3.
r23 = 2./3.
r12 = 1./12.

def ppm_limiters(dm, a4, itot, lmt):
    # INPUT PARAMETERS: 
    #real , intent(in):: dm(*)     !< Linear slope
    #integer, intent(in) :: itot      !< Total Longitudes
    #integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
    #                                      !< (Lin) 2: Positive definite constraint 
    #                                      !< 3: do nothing (return immediately)
    #! INPUT/OUTPUT PARAMETERS:
    #     real , intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
    # ! LOCAL VARIABLES:
    #      real  qmp
    #      real  da1, da2, a6da
    #      real  fmin
    #      integer i

    #! Developer: S.-J. Lin
    if (lmt == 3):
        return a4
        
    if (lmt == 0):
    #! Standard PPM constraint
        for i in range(0,itot):
            if(dm[i] == 0.):
                a4[1,i] = a4[0,i]
                a4[2,i] = a4[0,i]
                a4[3,i] = 0.
            else:
                da1  = a4[2,i] - a4[1,i]
                da2  = da1*da1
                a6da = a4[3,i]*da1
                if(a6da < -da2):
                    a4[3,i] = 3.*(a4[1,i]-a4[0,i])
                    a4[2,i] = a4[1,i] - a4[3,i]
                elif(a6da > da2):
                    a4[3,i] = 3.*(a4[2,i]-a4[0,i])
                    a4[1,i] = a4[2,i] - a4[3,i]
    elif (lmt == 1):
    #! Improved full monotonicity constraint (Lin 2004)
    #! Note: no need to provide first guess of A6 <-- a4(4,i)
        for i in range(0,itot):
            qmp = 2.*dm[i]
            a4[1,i] = a4[0,i]-np.sign(qmp)*np.abs(np.min([np.abs(qmp),np.abs(a4[1,i]-a4[0,i])]))
            a4[2,i] = a4[0,i]+np.sign(qmp)*np.abs(np.min([np.abs(qmp),np.abs(a4[2,i]-a4[0,i])]))
            a4[3,i] = 3.*( 2.*a4[0,i] - (a4[1,i]+a4[2,i]) )
    elif (lmt == 2):
    #! Positive definite constraint
         for i in range(0,itot):
             if( np.abs(a4[2,i]-a4[1,i]) < -a4[3,i] ):
                 fmin = a4[0,i]+0.25*(a4[2,i]-a4[1,i])**2/a4[3,i]+a4[3,i]*r12
                 if( fmin < 0.):
                    if(a4[0,i] < a4[2,i] and a4[0,i] < a4[1,i]):
                        a4[2,i] = a4[0,i]
                        a4[1,i] = a4[0,i]
                        a4[3,i] = 0.
                    elif(a4[2,i] > a4[1,i]):
                        a4[3,i] = 3.*(a4[1,i]-a4[0,i])
                        a4[2,i] = a4[1,i] - a4[3,i]
                    else:
                        a4[3,i] = 3.*(a4[2,i]-a4[0,i])
                        a4[1,i] = a4[2,i] - a4[3,i]
    return a4

def cs_limiters(im, extm, a4, iv):
    #integer, intent(in) :: im
    #integer, intent(in) :: iv
    #logical, intent(in) :: extm(im)
    #real , intent(inout) :: a4(4,im)   !< PPM array
    #! LOCAL VARIABLES:
    #real  da1, da2, a6da
    #integer i
    
    if (iv == 0):
        #! Positive definite constraint
        for i in range(0,im):
            if (a4[0,i] <= 0.):
                a4[1,i] = a4[0,i]
                a4[2,i] = a4[0,i]
                a4[3,i] = 0.
            else:
                if (np.abs(a4[2,i]-a4[1,i]) < -a4[3,i]):
                    if ((a4[0,i]+0.25*(a4[2,i]-a4[1,i])**2/a4[3,i]+a4[3,i]*r12) < 0.):
                        #! local minimum is negative
                        if (a4[0,i] < a4[2,i] and a4[0,i] < a4[1,i]):
                            a4[2,i] = a4[0,i]
                            a4[1,i] = a4[0,i]
                            a4[3,i] = 0.
                        elif (a4[2,i] > a4[1,i]):
                            a4[3,i] = 3.*(a4[1,i]-a4[0,i])
                            a4[2,i] = a4[1,i] - a4[3,i]
                        else:
                            a4[3,i] = 3.*(a4[2,i]-a4[0,i])
                            a4[1,i] = a4[2,i] - a4[3,i]
    elif (iv == 1):
        for i in range(0,im):
            if ((a4[0,i]-a4[1,i])*(a4[0,i]-a4[2,i]) >= 0.):
                a4[1,i] = a4[0,i]
                a4[2,i] = a4[0,i]
                a4[3,i] = 0.
            else:
                da1  = a4[2,i] - a4[1,i]
                da2  = da1**2
                a6da = a4[3,i]*da1
                if (a6da < -da2):
                    a4[3,i] = 3.*(a4[1,i]-a4[0,i])
                    a4[2,i] = a4[1,i] - a4[3,i]
                elif (a6da > da2):
                    a4[3,i] = 3.*(a4[2,i]-a4[0,i])
                    a4[1,i] = a4[2,i] - a4[3,i]
    else:
        #! Standard PPM constraint
        for i in range(0,im):
            if (extm[i]):
                a4[1,i] = a4[0,i]
                a4[2,i] = a4[0,i]
                a4[3,i] = 0.
            else:
                da1  = a4[2,i] - a4[1,i]
                da2  = da1**2
                a6da = a4[3,i]*da1
                if (a6da < -da2):
                    a4[3,i] = 3.*(a4[1,i]-a4[0,i])
                    a4[2,i] = a4[1,i] - a4[3,i]
                elif (a6da > da2):
                    a4[3,i] = 3.*(a4[2,i]-a4[0,i])
                    a4[1,i] = a4[2,i] - a4[3,i]
    return a4

def ppm_profile(a4, delp, km, i1, i2, iv, kord):

     #! INPUT PARAMETERS:
     #integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others iv = 2: temp (if remap_t) and w (iv=-2)
     #integer, intent(in):: i1      !< Starting longitude
     #integer, intent(in):: i2      !< Finishing longitude
     #integer, intent(in):: km      !< Vertical dimension
     #integer, intent(in):: kord    !< Order (or more accurately method no.):
     #real , intent(in):: delp(i1:i2,km)     !< Layer pressure thickness
     #!INPUT/OUTPUT PARAMETERS:
     #real , intent(inout):: a4(4,i1:i2,km)  !< Interpolated values
     #! DESCRIPTION:
     #!
     #!   Perform the piecewise parabolic reconstruction
     #! 
     #! !REVISION HISTORY: 
     #! S.-J. Lin   revised at GFDL 2007
     #!-----------------------------------------------------------------------
     #! local arrays:
     it = i2 - i1 + 1
     
     dc   = np.zeros((it,km))
     h2   = np.zeros((it,km))
     delq = np.zeros((it,km))
     df2  = np.zeros((it,km))
     d4   = np.zeros((it,km))
     #real    dc(i1:i2,km)
     #real    h2(i1:i2,km)
     #real  delq(i1:i2,km)
     #real   df2(i1:i2,km)
     #real    d4(i1:i2,km)

     #! local scalars:
     #integer i, k, km1, lmt, it
     #real  fac
     #real  a1, a2, c1, c2, c3, d1, d2
     #real  qm, dq, lac, qmp, pmp
     
     
     
     km1 = km - 1
     
     for k in range(1,km):
         for i in range(i1-1,i2):
             delq[i,k-1] =   a4[0,i,k] - a4[0,i,k-1]
             d4[i,k  ]   = delp[i,k-1] + delp[i,k]
     
     for k in range(1,km1):
         for i in range(i1-1,i2):
             c1  = (delp[i,k-1]+0.5*delp[i,k])/d4[i,k+1]
             c2  = (delp[i,k+1]+0.5*delp[i,k])/d4[i,k]
             df2[i,k] = delp[i,k]*(c1*delq[i,k] + c2*delq[i,k-1]) / (d4[i,k]+delp[i,k+1])
             dc[i,k] = np.sign(df2[i,k])*np.abs(np.min([np.abs(df2[i,k]), np.max([a4[0,i,k-1],a4[0,i,k],a4[0,i,k+1]])-a4[0,i,k], a4[0,i,k]-np.min([a4[0,i,k-1],a4[0,i,k],a4[0,i,k+1]])]))
     
    #!-----------------------------------------------------------
    #! 4th order interpolation of the provisional cell edge value
    #!-----------------------------------------------------------

     for k in range(2,km1):
         for i in range(i1-1,i2):
             c1 = delq[i,k-1]*delp[i,k-1] / d4[i,k]
             a1 = d4[i,k-1] / (d4[i,k] + delp[i,k-1])
             a2 = d4[i,k+1] / (d4[i,k] + delp[i,k])
             a4[1,i,k] = a4[0,i,k-1] + c1 + 2./(d4[i,k-1]+d4[i,k+1]) * (delp[i,k]*(c1*(a1 - a2)+a2*dc[i,k-1]) - delp[i,k-1]*a1*dc[i,k])
     
    #! Area preserving cubic with 2nd deriv. = 0 at the boundaries
    #! Top
     for i in range(i1-1,i2):
         d1 = delp[i,1]
         d2 = delp[i,2]
         qm = (d2*a4[0,i,0]+d1*a4[0,i,1]) / (d1+d2)
         dq = 2.*(a4[0,i,1]-a4[0,i,0]) / (d1+d2)
         c1 = 4.*(a4[1,i,2]-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4[1,i,1] = qm - 0.25*c1*d1*d2*(d2+3.*d1)
         #! Top edge:
         #!-------------------------------------------------------
         a4[1,i,0] = d1*(2.*c1*d1**2-c3) + a4[1,i,1]
         #!-------------------------------------------------------
         #!        a4[2,i,1] = (12./7.)*a4[1,i,1]-(13./14.)*a4[1,i,2]+(3./14.)*a4[1,i,3]
         #!-------------------------------------------------------
         #! No over- and undershoot condition
         a4[1,i,1] = np.max([a4[1,i,1], np.min([a4[0,i,0], a4[0,i,1]])])
         a4[1,i,1] = np.min([a4[1,i,1], np.max([a4[0,i,0], a4[0,i,1]])])
         dc[i,0] =  0.5*(a4[1,i,1] - a4[0,i,0])

     #! Enforce monotonicity  within the top layer
     
     if (iv == 0):
        for i in range(i1-1,i2):
           a4[1,i,0] = np.max([0., a4[1,i,0]])
           a4[1,i,1] = np.max([0., a4[1,i,1]])
     elif (iv == -1):
         for i in range(i1-1,i2):
             if (a4[1,i,0]*a4[0,i,0] <= 0. ):
                  a4[1,i,0] = 0.
     elif (np.abs(iv) == 2):
         for i in range(i1-1,i2):
             a4[1,i,0] = a4[0,i,0]
             a4[2,i,0] = a4[0,i,0]
      
     #! Bottom
     #! Area preserving cubic with 2nd deriv. = 0 at the surface
     for i in range(i1-1,i2):
         d1 = delp[i,km-1]
         d2 = delp[i,km1-1]
         qm = (d2*a4[0,i,km-1]+d1*a4[0,i,km1-1]) / (d1+d2)
         dq = 2.*(a4[0,i,km1-1]-a4[0,i,km-1]) / (d1+d2)
         c1 = (a4[1,i,km1-1]-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4[1,i,km-1] = qm - c1*d1*d2*(d2+3.*d1)
         #! Bottom edge:
         #!-----------------------------------------------------
         a4[2,i,km-1] = d1*(8.*c1*d1**2-c3) + a4[1,i,km-1]
         #!        dc[i,km] = 0.5*(a4[3,i,km] - a4[1,i,km])
         #!-----------------------------------------------------
         #!        a4[3,i,km] = (12./7.)*a4[1,i,km]-(13./14.)*a4[1,i,km-1]+(3./14.)*a4[1,i,km-2]
         #! No over- and under-shoot condition
         a4[1,i,km-1] = np.max([a4[1,i,km-1], np.min([a4[0,i,km-1], a4[0,i,km1-1]])])
         a4[1,i,km-1] = np.min([a4[1,i,km-1], np.max([a4[0,i,km-1], a4[0,i,km1-1]])])
         dc[i,km-1] = 0.5*(a4[0,i,km-1] - a4[1,i,km-1])

     #! Enforce constraint on the "slope" at the surface

     ##ifdef BOT_MONO
     #     do i=i1,i2
     #        a4(4,i,km) = 0
     #        if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
     #        d1 = a4(1,i,km) - a4(2,i,km)
     #        d2 = a4(3,i,km) - a4(1,i,km)
     #        if ( d1*d2 < 0. ) then
     #             a4(2,i,km) = a4(1,i,km)
     #             a4(3,i,km) = a4(1,i,km)
     #        else
     #             dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
     #             a4(2,i,km) = a4(1,i,km) - dq
     #             a4(3,i,km) = a4(1,i,km) + dq
     #        endif
     #     enddo
     ##else
     if (iv == 0):
         for i in range(i1-1,i2):
            a4[1,i,km-1] = np.max([0.,a4[1,i,km-1]])
            a4[2,i,km-1] = np.max([0.,a4[2,i,km-1]])
     elif (iv < 0):
         for i in range(i1-1,i2):
             if (a4[0,i,km-1]*a4[2,i,km-1] <= 0.):
                 a4[2,i,km-1] = 0.
     ##endif
     
     for k in range(0,km1):
         for i in range(i1-1,i2):
            a4[2,i,k] = a4[1,i,k+1]

     #!-----------------------------------------------------------
     #! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
     #!-----------------------------------------------------------
     #! Top 2 and bottom 2 layers always use monotonic mapping
     for k in range(0,2):
         for i in range(i1-1,i2):
             a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
         a4[:,:,k] = ppm_limiters(dc[:,k], a4[:,:,k], it, 0)

     if (kord >= 7):
         #!-----------------------
         #! Huynh's 2nd constraint
         #!-----------------------
         for k in range(1,km1):
             for i in range(i1-1,i2):
                 #! Method#1
                 #!           h2[i,k] = delq[i,k] - delq[i,k-1]
                 #! Method#2 - better
                 h2[i,k] = 2.*(dc[i,k+1]/delp[i,k+1] - dc[i,k-1]/delp[i,k-1]) / (delp[i,k]+0.5*(delp[i,k-1]+delp[i,k+1])) * delp[i,k]**2
                 #! Method#3
                 #!!!         h2[i,k] = dc[i,k+1] - dc[i,k-1]
         fac = 1.5           #! original quasi-monotone
         
         for k in range(2,km-2):
             for i in range(i1-1,i2):
                 #! Right edges
                 #!        qmp   = a4[1,i,k] + 2.0*delq[i,k-1]
                 #!        lac   = a4[1,i,k] + fac*h2[i,k-1] + 0.5*delq[i,k-1]
                 pmp   = 2.*dc[i,k]
                 qmp   = a4[0,i,k] + pmp
                 lac   = a4[0,i,k] + fac*h2[i,k-1] + dc[i,k]
                 a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], qmp, lac])]), np.max([a4[0,i,k], qmp, lac])])
                 #! Left  edges
                 #!        qmp   = a4[1,i,k] - 2.0*delq[i,k]
                 #!        lac   = a4[1,i,k] + fac*h2[i,k+1] - 0.5*delq[i,k]
                 #!
                 qmp   = a4[0,i,k] - pmp
                 lac   = a4[0,i,k] + fac*h2[i,k+1] - dc[i,k]
                 a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], qmp, lac])]), np.max([a4[0,i,k], qmp, lac])])
                 #!-------------
                 #! Recompute A6
                 #!-------------
                 a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
             #! Additional constraint to ensure positivity when kord=7
             if (iv == 0 and kord >= 6):
                 a4[:,:,k] = ppm_limiters(dc[:,k], a4[:,:,k], it, 2)
     else:
        lmt = kord - 3
        lmt = np.max([0, lmt])
        if (iv == 0):
            lmt = np.min([2, lmt])

        for k in range(2,km-2):
            if( kord != 4):
                for i in range(i1-1,i2):
                    a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
             
            if(kord != 6):
                 a4[:,:,k] = ppm_limiters(dc[:,k], a4[:,:,k], it, lmt)
        
     for k in range(km1-1,km):
         for i in range(i1-1,i2):
             a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
         a4[:,:,k] = ppm_limiters(dc[:,k], a4[:,:,k], it, 0)

     return a4

def scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin):
    #! Optimized vertical profile reconstruction:
    #! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
    #integer, intent(in):: i1, i2
    #integer, intent(in):: km      !< vertical dimension
    #integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others
    #integer, intent(in):: kord
    #real, intent(in)   ::   qs(i1:i2)
    #real, intent(in)   :: delp(i1:i2,km)     !< Layer pressure thickness
    #real, intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
    #real, intent(in):: qmin
    #!-----------------------------------------------------------------------
    im = i2 - i1 + 1
    extm = np.zeros([im,km],dtype=bool)
    ext5 = np.zeros([im,km],dtype=bool)
    ext6 = np.zeros([im,km],dtype=bool)
    
    gam = np.zeros([im,km])
    q   = np.zeros([im,km+1])
    d4  = np.zeros([im])
    
    #logical, dimension(i1:i2,km):: extm, ext5, ext6
    #real  gam(i1:i2,km)
    #real    q(i1:i2,km+1)
    #real   d4(i1:i2)
    #real   bet, a_bot, grat 
    #real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
    #integer i, k, im
    
    if (iv == -2):
        for i in range(0,im):
            gam[i,1] = 0.5
            q[i,0] = 1.5*a4[0,i,0]
        for k in range(1,km-1):
            for i in range(0,im):
                grat = delp[i,k-1] / delp[i,k]
                bet =  2. + grat + grat - gam[i,k]
                q[i,k] = (3.*(a4[0,i,k-1]+a4[0,i,k]) - q[i,k-1])/bet
                gam[i,k+1] = grat / bet
        for i in range(0,im):
           grat = delp[i,km-2] / delp[i,km-1]
           q[i,km-1] = (3.*(a4[0,i,km-2]+a4[0,i,km-1]) - grat*qs[i] - q[i,km-2]) / (2. + grat + grat - gam[i,km-1])
           q[i,km] = qs[i]
        for k in range(km-2,-1,-1):
            for i in range(0,im):
                q[i,k] = q[i,k] - gam[i,k+1]*q[i,k+1]
    else:
        for i in range(0,im):
            grat = delp[i,1] / delp[i,0]   #! grid ratio
            bet = grat*(grat+0.5)
            q[i,0] = ((grat+grat)*(grat+1.)*a4[0,i,0] + a4[0,i,1]) / bet
            gam[i,0] = ( 1. + grat*(grat+1.5) ) / bet
        for k in range(1,km):
            for i in range(0,im):
                d4[i] = delp[i,k-1] / delp[i,k]
                bet =  2. + d4[i] + d4[i] - gam[i,k-1]
                q[i,k] = ( 3.*(a4[0,i,k-1]+d4[i]*a4[0,i,k]) - q[i,k-1] )/bet
                gam[i,k] = d4[i] / bet
        for i in range(0,im):
            a_bot = 1. + d4[i]*(d4[i]+1.5)
            q[i,km] = (2.*d4[i]*(d4[i]+1.)*a4[0,i,km-1]+a4[0,i,km-2]-a_bot*q[i,km-1]) / ( d4[i]*(d4[i]+0.5) - a_bot*gam[i,km-1])
        for k in range(km-1,-1,-1):
            for i in range(0,im):
                q[i,k] = q[i,k] - gam[i,k]*q[i,k+1]
    

    #!----- Perfectly linear scheme --------------------------------
    if (np.abs(kord) > 16):
        for k in range(0,km):
            for i in range(0,im):
                a4[1,i,k] = q[i,k]
                a4[2,i,k] = q[i,k+1]
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        return a4

    #!----- Perfectly linear scheme --------------------------------
    #!------------------
    #! Apply constraints
    #!------------------
    
    #! Apply *large-scale* constraints 
    for i in range(0,im):
        q[i,1] = np.min([q[i,1], np.max([a4[0,i,0], a4[0,i,1]])])
        q[i,1] = np.max([q[i,1], np.min([a4[0,i,0], a4[0,i,1]])])
    
    for k in range(1,km):
        for i in range(0,im):
            gam[i,k] = a4[0,i,k] - a4[0,i,k-1]
    
    #! Interior:
    for k in range(2,km-1):
        for i in range(0,im):
            if (gam[i,k-1]*gam[i,k+1] > 0.):
                #! Apply large-scale constraint to ALL fields if not local max/min
                q[i,k] = np.min([q[i,k], np.max([a4[0,i,k-1],a4[0,i,k]])])
                q[i,k] = np.max([q[i,k], np.min([a4[0,i,k-1],a4[0,i,k]])])
            else:
                if (gam[i,k-1] > 0):
                    #! There exists a local max
                    q[i,k] = np.max([q[i,k], np.min([a4[0,i,k-1],a4[0,i,k]])])
                else:
                    #! There exists a local min
                    q[i,k] = np.min([q[i,k], np.max([a4[0,i,k-1],a4[0,i,k]])])
                    if (iv == 0):
                        q[i,k] = np.max([0., q[i,k]])

    #! Bottom:
    for i in range(0,im):
        q[i,km-1] = np.min([q[i,km-1], np.max([a4[0,i,km-2], a4[0,i,km-1]])])
        q[i,km-1] = np.max([q[i,km-1], np.min([a4[0,i,km-2], a4[0,i,km-1]])])

    for k in range(0,km):
        for i in range(0,im):
            a4[1,i,k] = q[i,k  ]
            a4[2,i,k] = q[i,k+1]
    
    for k in range(0,km):
        if (k == 0 or k == km-1):
            for i in range(0,im):
                extm[i,k] = (a4[1,i,k]-a4[0,i,k]) * (a4[2,i,k]-a4[0,i,k]) > 0.
        else:
            for i in range(0,im):
                extm[i,k] = gam[i,k]*gam[i,k+1] < 0.
        if ( np.abs(kord) > 9 ):
            for i in range(0,im):
                x0 = 2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k])
                x1 = np.abs(a4[1,i,k]-a4[2,i,k])
                a4[3,i,k] = 3.*x0
                ext5[i,k] = np.abs(x0) > x1
                ext6[i,k] = np.abs(a4[3,i,k]) > x1

    #!---------------------------
    #! Apply subgrid constraints:
    #!---------------------------
    #! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
    #! Top 2 and bottom 2 layers always use monotonic mapping

    if (iv == 0):
        for i in range(0,im):
            a4[1,i,0] = np.max([0., a4[1,i,0]])
    elif (iv == -1):
        for i in range(0,im):
            if ( a4[1,i,0]*a4[0,i,0] <= 0. ):
                a4[1,i,0] = 0.
    elif (iv == 2):
        for i in range(0,im):
            a4[1,i,0] = a4[0,i,0]
            a4[2,i,0] = a4[0,i,0]
            a4[3,i,0] = 0.
            
    if (iv != 2):
        for i in range(0,im):
            a4[3,i,0] = 3.*(2.*a4[0,i,0] - (a4[1,i,0]+a4[2,i,0]))
        a4[:,:,0] = cs_limiters(im, extm[:,0], a4[:,:,0], 1)
    
    #! k=1
    for i in range(0,im):
        a4[3,i,1] = 3.*(2.*a4[0,i,1] - (a4[1,i,1]+a4[2,i,1]))
    a4[:,:,1] = cs_limiters(im, extm[:,1], a4[:,:,1], 2)
    
    #!-------------------------------------
    #! Huynh's 2nd constraint for interior:
    #!-------------------------------------
    for k in range(2,km-2):
        if (np.abs(kord) < 9):
            for i in range(0,im):
                #! Left  edges
                pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                lac_1 = pmp_1 + 1.5*gam[i,k+2]
                a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                #! Right edges
                pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                lac_2 = pmp_2 - 1.5*gam[i,k-1]
                a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])

                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        elif (np.abs(kord) == 9):
            for i in range(0,im):
                if (extm[i,k] and extm[i,k-1]):
                    #! grid-scale 2-delta-z wave detected
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                    a4[3,i,k] = 0.
                elif (extm[i,k] and extm[i,k+1]):
                    #! grid-scale 2-delta-z wave detected
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                    a4[3,i,k] = 0.
                elif (extm[i,k] and a4[0,i,k] < qmin):
                    #! grid-scale 2-delta-z wave detected
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                    a4[3,i,k] = 0.
                else:
                    a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
                    #! Check within the smooth region if subgrid profile is non-monotonic
                    if(np.abs(a4[3,i,k]) > np.abs(a4[1,i,k]-a4[2,i,k])):
                        pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                        lac_1 = pmp_1 + 1.5*gam[i,k+2]
                        a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                        pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                        lac_2 = pmp_2 - 1.5*gam[i,k-1]
                        a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                        a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))

        elif (np.abs(kord) == 10):
            for i in range(0,im):
                if (ext5[i,k]):
                    if (ext5[i,k-1] or ext5[i,k+1]):
                        a4[1,i,k] = a4[0,i,k]
                        a4[2,i,k] = a4[0,i,k]
                    elif (ext6[i,k-1] or ext6[i,k+1]):
                        pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                        lac_1 = pmp_1 + 1.5*gam[i,k+2]
                        a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                        pmp_2 = a4[1,i,k] + 2.*gam[i,k]
                        lac_2 = pmp_2 - 1.5*gam[i,k-1]
                        a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                elif (ext6[i,k]):
                    if (ext5[i,k-1] or ext5[i,k+1]):
                        pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                        lac_1 = pmp_1 + 1.5*gam[i,k+2]
                        a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                        pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                        lac_2 = pmp_2 - 1.5*gam[i,k-1]
                        a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
            for i in range(0,im):
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        elif (np.abs(kord) == 12):
            for i in range(0,im):
                if (extm[i,k]):
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                    a4[3,i,k] = 0.
                else:        #! not a local extremum
                    a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])
                    #! Check within the smooth region if subgrid profile is non-monotonic
                    if (np.abs(a4[3,i,k]) > np.abs(a4[1,i,k]-a4[2,i,k])):
                        pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                        lac_1 = pmp_1 + 1.5*gam[i,k+2]
                        a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                        pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                        lac_2 = pmp_2 - 1.5*gam[i,k-1]
                        a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                        a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])
        elif (np.abs(kord) == 13):
            for i in range(0,im):
                if (ext6[i,k]):
                    if (ext6[i,k-1] and ext6[i,k+1]):
                        #! grid-scale 2-delta-z wave detected
                        a4[1,i,k] = a4[0,i,k]
                        a4[2,i,k] = a4[0,i,k]
            for i in range(0,im):
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        elif (np.abs(kord) == 14):
            for i in range(0,im):
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        elif (np.abs(kord) == 15):   #! Revised abs(kord)=9 scheme
            for i in range(0,im):
                if (ext5[i,k] and ext5[i,k-1]):
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                elif (ext5[i,k] and ext5[i,k+1]):
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                elif (ext5[i,k] and a4[0,i,k] < qmin):
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                elif (ext6[i,k]):
                    pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                    lac_1 = pmp_1 + 1.5*gam[i,k+2]
                    a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                    pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                    lac_2 = pmp_2 - 1.5*gam[i,k-1]
                    a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
            for i in range(0,im):
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        elif (np.abs(kord) == 16):
            for i in range(0,im):
                if (ext5[i,k]):
                    if (ext5[i,k-1] or ext5[i,k+1]):
                        a4[1,i,k] = a4[0,i,k]
                        a4[2,i,k] = a4[0,i,k]
                    elif (ext6[i,k-1] or ext6[i,k+1]):
                        #! Left  edges
                        pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                        lac_1 = pmp_1 + 1.5*gam[i,k+2]
                        a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                        #! Right edges
                        pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                        lac_2 = pmp_2 - 1.5*gam[i,k-1]
                        a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
            for i in range(0,im):
                a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        else:      #! kord = 11, 13
            for i in range(0,im):
                if (ext5[i,k] and (ext5[i,k-1] or ext5[i,k+1] or a4[0,i,k] < qmin)):
                    #! Noisy region:
                    a4[1,i,k] = a4[0,i,k]
                    a4[2,i,k] = a4[0,i,k]
                    a4[3,i,k] = 0.
                else:
                    a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
                    
        #! Additional constraint to ensure positivity
        if (iv == 0):
            a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 0)

    ####end for k in range(3,km-2)

    #!----------------------------------
    #! Bottom layer subgrid constraints:
    #!----------------------------------
    if (iv == 0):
        for i in range(0,im):
            a4[2,i,km-1] = np.max([0., a4[2,i,km-1]])
    elif (iv == -1):
        for i in range(0,im): 
            if (a4[2,i,km-1]*a4[0,i,km-1] <= 0.):
                a4[2,i,km-1] = 0.

    for k in range(km-2,km):
        for i in range(0,im):
            a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
        if (k == (km-2)):
            a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 2)
        if (k == km-1):
            a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 1)

    return a4

def map_scalar(km, pe1, q1, qs, kn, pe2, i1, i2, iv, kord, q_min):
    #! iv=1
    #integer, intent(in) :: i1                !< Starting longitude
    #integer, intent(in) :: i2                !< Finishing longitude
    #integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == temp 2 == remap temp with cs scheme
    #integer, intent(in) :: kord              !< Method order
    #integer, intent(in) :: km                !< Original vertical dimension
    #integer, intent(in) :: kn                !< Target vertical dimension
    #real, intent(in) ::   qs(i1:i2)       !< bottom BC
    #real, intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
    #real, intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
    #real, intent(in) ::    q1(ibeg:iend,km) !< Field input
    #! INPUT/OUTPUT PARAMETERS:
    #real, intent(inout)::  q2(ibeg:iend,kn) !< Field output
    
    im = i2 - i1 + 1
    q2 = np.zeros([im,kn])
    #real, intent(in):: q_min

    #! DESCRIPTION:
    #! IV = 0: constituents
    #! pe1: pressure at layer edges (from model top to bottom surface)
    #!      in the original vertical coordinate
    #! pe2: pressure at layer edges (from model top to bottom surface)
    #!      in the new vertical coordinate
    #! LOCAL VARIABLES:
    dp1 = np.zeros([im,km])
    q4  = np.zeros([4,im,km])
    #real    dp1(i1:i2,km)
    #real   q4(4,i1:i2,km)
    #real    pl, pr, qsum, dp, esl
    #integer i, k, l, m, k0
    qsum = 0.

    for k in range(0,km):
        for i in range(0,im):
            dp1[i,k] = pe1[i,k+1] - pe1[i,k]
            q4[0,i,k] = q1[i,k]
    
    #! Compute vertical subgrid distribution
    if (kord > 7):
       #print qs, q4, dp1, km, i1, i2, iv, kord, q_min
       q4 = scalar_profile(qs, q4, dp1, km, i1, i2, iv, kord, q_min)
    else:
       q4 = ppm_profile(q4, dp1, km, i1, i2, iv, kord)
       
    for i in range(0,im):
        k0 = 0
        for k in range(0,kn):
            next_k = False
            for l in range(k0,km):  #AKA l-loop
                #! locate the top edge: pe2(i,k)
                if (pe2[i,k] >= pe1[i,l] and pe2[i,k] <= pe1[i,l+1]):
                    pl = (pe2[i,k]-pe1[i,l]) / dp1[i,l]
                    if (pe2[i,k+1] <= pe1[i,l+1]):
                        #! entire new grid is within the original grid
                        pr = (pe2[i,k+1]-pe1[i,l]) / dp1[i,l]
                        q2[i,k] = q4[1,i,l] + 0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l]) * (pr+pl)-q4[3,i,l]*r3*(pr*(pr+pl)+pl**2)
                        k0 = l
                        next_k = True
                        break
                        #goto 555 #(next iteration of "for k in range(0,kn):" loop)
                    else:
                        #! Fractional area...
                        qsum = (pe1[i,l+1]-pe2[i,k])*(q4[1,i,l]+0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(1.+pl)-q4[3,i,l]*(r3*(1.+pl*(1.+pl))))
                        for m in range(l+1,km): #AKA m-loop
                            #! locate the bottom edge: pe2(i,k+1)
                            if (pe2[i,k+1] > pe1[i,m+1]):
                                #! Whole layer
                                qsum = qsum + dp1[i,m]*q4[0,i,m]
                            else:
                                dp = pe2[i,k+1]-pe1[i,m]
                                esl = dp / dp1[i,m]
                                qsum = qsum + dp*(q4[1,i,m]+0.5*esl*(q4[2,i,m]-q4[1,i,m]+q4[3,i,m]*(1.-r23*esl)))
                                k0 = m
                                #goto 123 #(exit out of l-loop)
                                break
                        else:
                            #GJF: the following if statement is not in the fv_mapz, but it captures the case where pe2[kn] > pe1[km] where the m loop is not entered; without this, the lowest layer values are weird
                            if (l+1 == km):
                                dp = pe2[i,kn]-pe1[i,km]
                                esl = dp / dp1[i,km-1]
                                qsum = qsum + dp*(q4[1,i,km-1]+0.5*esl*(q4[2,i,km-1]-q4[1,i,km-1]+q4[3,i,km-1]*(1.-r23*esl)))
                            break #handles goto 123 statement below (exits out of l-loop even if m-loop successfully completes)
                            #continue
                        break
                        #goto 123 #(right before going to next iteration of "for k in range(1,kn):" loop)
            if not next_k:
                q2[i,k] = qsum / (pe2[i,k+1] - pe2[i,k]) #AKA label 123
                
    return q2

def map1_q2 (km, pe1, q1, kn, pe2, dp2, i1, i2, iv, kord, q_min):
    #! INPUT PARAMETERS:
    #integer, intent(in) :: i1, i2
    #integer, intent(in) :: iv                !< Mode: 0 ==  constituents 1 == ???
    #integer, intent(in) :: kord
    #integer, intent(in) :: km                !< Original vertical dimension
    #integer, intent(in) :: kn                !< Target vertical dimension
    #real, intent(in) ::  pe1(i1:i2,km+1)     !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
    #real, intent(in) ::  pe2(i1:i2,kn+1)     !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
    #real, intent(in) ::  q1(i1:i2,km) !< Field input
    #real, intent(in) ::  dp2(i1:i2,kn)
    #real, intent(in) ::  q_min
    #! INPUT/OUTPUT PARAMETERS:
    im = i2 - i1 + 1
    q2 = np.zeros([im,kn])
    #real, intent(inout):: q2(i1:i2,kn) !< Field output
    #! LOCAL VARIABLES:
    im = i2 - i1 + 1
    qs = np.zeros([im])
    dp1 = np.zeros([im,km])
    q4 = np.zeros([4,im,km])
    #real   qs(i1:i2)
    #real   dp1(i1:i2,km)
    #real   q4(4,i1:i2,km)
    #real   pl, pr, qsum, dp, esl
    #integer i, k, l, m, k0
    qsum = 0.
    for k in range(0,km):
        for i in range(0,im):
            dp1[i,k] = pe1[i,k+1] - pe1[i,k]
            q4[0,i,k] = q1[i,k]
    
    #! Compute vertical subgrid distribution
    if (kord > 7):
       q4 = scalar_profile (qs, q4, dp1, km, i1, i2, iv, kord, q_min)
    else:
       q4 = ppm_profile (q4, dp1, km, i1, i2, iv, kord)
    
    #! Mapping
    for i in range(0,im):
        k0 = 0
        for k in range(0,kn):
            next_k = False
            #print 'k new = ',k
            for l in range(k0,km):
                #print 'l old = ',l
                #! locate the top edge: pe2(i,k)
                if (pe2[i,k] >= pe1[i,l] and pe2[i,k] <= pe1[i,l+1]):
                    pl = (pe2[i,k]-pe1[i,l]) / dp1[i,l]
                    if (pe2[i,k+1] <= pe1[i,l+1]):
                        #! entire new grid is within the original grid
                        pr = (pe2[i,k+1]-pe1[i,l]) / dp1[i,l]
                        q2[i,k] = q4[1,i,l] + 0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(pr+pl)-q4[3,i,l]*r3*(pr*(pr+pl)+pl**2)
                        k0 = l
                        next_k = True
                        #print 'new grid within old; q2 = ', q2[i,k]
                        break
                        #goto 555 #next k-loop iteration
                    else:
                        #! Fractional area...
                        #print k, (pe1[i,l+1]-pe2[i,k]), (q4[1,i,l]+0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(1.+pl)-q4[3,i,l]*(r3*(1.+pl*(1.+pl)))), dp2[i,k]
                        qsum = (pe1[i,l+1]-pe2[i,k])*(q4[1,i,l]+0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(1.+pl)-q4[3,i,l]*(r3*(1.+pl*(1.+pl))))
                        for m in range(l+1,km):
                            #! locate the bottom edge: pe2(i,k+1)
                            if (pe2[i,k+1] > pe1[i,m+1]):
                                #! Whole layer..
                                qsum = qsum + dp1[i,m]*q4[0,i,m]
                                #print 'whole layer, m = ',m
                            else:
                                dp = pe2[i,k+1]-pe1[i,m]
                                esl = dp / dp1[i,m]
                                qsum = qsum + dp*(q4[1,i,m]+0.5*esl*(q4[2,i,m]-q4[1,i,m]+q4[3,i,m]*(1.-r23*esl)))
                                k0 = m
                                #print 'partial layer, m = ',m
                                #goto 123 #end l-loop
                                break
                        else:
                            #GJF: the following if statement is not in the fv_mapz, but it captures the case where pe2[kn] > pe1[km] where the m loop is not entered; without this, the lowest layer values are weird
                            if (l+1 == km):
                                dp = pe2[i,kn]-pe1[i,km]
                                esl = dp / dp1[i,km-1]
                                qsum = qsum + dp*(q4[1,i,km-1]+0.5*esl*(q4[2,i,km-1]-q4[1,i,km-1]+q4[3,i,km-1]*(1.-r23*esl)))
                            break
                        
                        break
                        #goto 123 #end l-loop
            if not next_k:
                q2[i,k] = qsum / dp2[i,k] #formerly labeled 123
                #print 'result q2 ', q2[i,k]
    #print q2
    #exit()
    return q2

def cs_profile(qs, a4, delp, km, i1, i2, iv, kord):
      #! Optimized vertical profile reconstruction:
      #! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
      #integer, intent(in):: i1, i2
      #integer, intent(in):: km      !< vertical dimension
      #integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others
      #integer, intent(in):: kord
      #real, intent(in)   ::   qs(i1:i2)
      #real, intent(in)   :: delp(i1:i2,km)     !< Layer pressure thickness
      #real, intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
      #real, intent(in):: qmin
      #!-----------------------------------------------------------------------
      im = i2 - i1 + 1
      extm = np.zeros([im,km],dtype=bool)
      ext5 = np.zeros([im,km],dtype=bool)
      ext6 = np.zeros([im,km],dtype=bool)
      
      gam = np.zeros([im,km])
      q   = np.zeros([im,km+1])
      d4  = np.zeros([im])
      
      #logical, dimension(i1:i2,km):: extm, ext5, ext6
      #real  gam(i1:i2,km)
      #real    q(i1:i2,km+1)
      #real   d4(i1:i2)
      #real   bet, a_bot, grat 
      #real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
      #integer i, k, im
      
      if (iv == -2):
          for i in range(0,im):
              gam[i,1] = 0.5
              q[i,0] = 1.5*a4[0,i,0]
          for k in range(1,km-1):
              for i in range(0,im):
                  grat = delp[i,k-1] / delp[i,k]
                  bet =  2. + grat + grat - gam[i,k]
                  q[i,k] = (3.*(a4[0,i,k-1]+a4[0,i,k]) - q[i,k-1])/bet
                  gam[i,k+1] = grat / bet
          for i in range(0,im):
             grat = delp[i,km-2] / delp[i,km-1]
             q[i,km-1] = (3.*(a4[0,i,km-2]+a4[0,i,km-1]) - grat*qs[i] - q[i,km-2]) / (2. + grat + grat - gam[i,km-1])
             q[i,km] = qs[i]
          for k in range(km-2,-1,-1):
              for i in range(0,im):
                  q[i,k] = q[i,k] - gam[i,k+1]*q[i,k+1]
      else:
          for i in range(0,im):
              grat = delp[i,1] / delp[i,0]   #! grid ratio
              bet = grat*(grat+0.5)
              q[i,0] = ((grat+grat)*(grat+1.)*a4[0,i,0] + a4[0,i,1]) / bet
              gam[i,0] = ( 1. + grat*(grat+1.5) ) / bet
          for k in range(1,km):
              for i in range(0,im):
                  d4[i] = delp[i,k-1] / delp[i,k]
                  bet =  2. + d4[i] + d4[i] - gam[i,k-1]
                  q[i,k] = ( 3.*(a4[0,i,k-1]+d4[i]*a4[0,i,k]) - q[i,k-1] )/bet
                  gam[i,k] = d4[i] / bet
          for i in range(0,im):
              a_bot = 1. + d4[i]*(d4[i]+1.5)
              q[i,km] = (2.*d4[i]*(d4[i]+1.)*a4[0,i,km-1]+a4[0,i,km-2]-a_bot*q[i,km-1]) / ( d4[i]*(d4[i]+0.5) - a_bot*gam[i,km-1])
          for k in range(km-1,-1,-1):
              for i in range(0,im):
                  q[i,k] = q[i,k] - gam[i,k]*q[i,k+1]
      

      #!----- Perfectly linear scheme --------------------------------
      if (np.abs(kord) > 16):
          for k in range(0,km):
              for i in range(0,im):
                  a4[1,i,k] = q[i,k]
                  a4[2,i,k] = q[i,k+1]
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          return a4

      #!----- Perfectly linear scheme --------------------------------
      #!------------------
      #! Apply constraints
      #!------------------
      
      #! Apply *large-scale* constraints 
      for i in range(0,im):
          q[i,1] = np.min([q[i,1], np.max([a4[0,i,0], a4[0,i,1]])])
          q[i,1] = np.max([q[i,1], np.min([a4[0,i,0], a4[0,i,1]])])
      
      for k in range(1,km):
          for i in range(0,im):
              gam[i,k] = a4[0,i,k] - a4[0,i,k-1]
      
      #! Interior:
      for k in range(2,km-1):
          for i in range(0,im):
              if (gam[i,k-1]*gam[i,k+1] > 0.):
                  #! Apply large-scale constraint to ALL fields if not local max/min
                  q[i,k] = np.min([q[i,k], np.max([a4[0,i,k-1],a4[0,i,k]])])
                  q[i,k] = np.max([q[i,k], np.min([a4[0,i,k-1],a4[0,i,k]])])
              else:
                  if (gam[i,k-1] > 0):
                      #! There exists a local max
                      q[i,k] = np.max([q[i,k], np.min([a4[0,i,k-1],a4[0,i,k]])])
                  else:
                      #! There exists a local min
                      q[i,k] = np.min([q[i,k], np.max([a4[0,i,k-1],a4[0,i,k]])])
                      if (iv == 0):
                          q[i,k] = np.max([0., q[i,k]])

      #! Bottom:
      for i in range(0,im):
          q[i,km-1] = np.min([q[i,km-1], np.max([a4[0,i,km-2], a4[0,i,km-1]])])
          q[i,km-1] = np.max([q[i,km-1], np.min([a4[0,i,km-2], a4[0,i,km-1]])])

      for k in range(0,km):
          for i in range(0,im):
              a4[1,i,k] = q[i,k  ]
              a4[2,i,k] = q[i,k+1]
      
      for k in range(0,km):
          if (k == 0 or k == km-1):
              for i in range(0,im):
                  extm[i,k] = (a4[1,i,k]-a4[0,i,k]) * (a4[2,i,k]-a4[0,i,k]) > 0.
          else:
              for i in range(0,im):
                  extm[i,k] = gam[i,k]*gam[i,k+1] < 0.
          if ( np.abs(kord) > 9 ):
              for i in range(0,im):
                  x0 = 2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k])
                  x1 = np.abs(a4[1,i,k]-a4[2,i,k])
                  a4[3,i,k] = 3.*x0
                  ext5[i,k] = np.abs(x0) > x1
                  ext6[i,k] = np.abs(a4[3,i,k]) > x1

      #!---------------------------
      #! Apply subgrid constraints:
      #!---------------------------
      #! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
      #! Top 2 and bottom 2 layers always use monotonic mapping

      if (iv == 0):
          for i in range(0,im):
              a4[1,i,0] = np.max([0., a4[1,i,0]])
      elif (iv == -1):
          for i in range(0,im):
              if ( a4[1,i,0]*a4[0,i,0] <= 0. ):
                  a4[1,i,0] = 0.
      elif (iv == 2):
          for i in range(0,im):
              a4[1,i,0] = a4[0,i,0]
              a4[2,i,0] = a4[0,i,0]
              a4[3,i,0] = 0.
              
      if (iv != 2):
          for i in range(0,im):
              a4[3,i,0] = 3.*(2.*a4[0,i,0] - (a4[1,i,0]+a4[2,i,0]))
          a4[:,:,0] = cs_limiters(im, extm[:,0], a4[:,:,0], 1)
      
      #! k=1
      for i in range(0,im):
          a4[3,i,1] = 3.*(2.*a4[0,i,1] - (a4[1,i,1]+a4[2,i,1]))
      a4[:,:,1] = cs_limiters(im, extm[:,1], a4[:,:,1], 2)
      
      #!-------------------------------------
      #! Huynh's 2nd constraint for interior:
      #!-------------------------------------
      for k in range(2,km-2):
          if (np.abs(kord) < 9):
              for i in range(0,im):
                  #! Left  edges
                  pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                  lac_1 = pmp_1 + 1.5*gam[i,k+2]
                  a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                  #! Right edges
                  pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                  lac_2 = pmp_2 - 1.5*gam[i,k-1]
                  a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])

                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          elif (np.abs(kord) == 9):
              for i in range(0,im):
                  if (extm[i,k] and extm[i,k-1]):
                      #! grid-scale 2-delta-z wave detected
                      a4[1,i,k] = a4[0,i,k]
                      a4[2,i,k] = a4[0,i,k]
                      a4[3,i,k] = 0.
                  elif (extm[i,k] and extm[i,k+1]):
                      #! grid-scale 2-delta-z wave detected
                      a4[1,i,k] = a4[0,i,k]
                      a4[2,i,k] = a4[0,i,k]
                      a4[3,i,k] = 0.
                  else:
                      a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])
                      #! Check within the smooth region if subgrid profile is non-monotonic
                      if(np.abs(a4[3,i,k]) > np.abs(a4[1,i,k]-a4[2,i,k])):
                          pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                          lac_1 = pmp_1 + 1.5*gam[i,k+2]
                          a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                          pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                          lac_2 = pmp_2 - 1.5*gam[i,k-1]
                          a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                          a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])

          elif (np.abs(kord) == 10):
              for i in range(0,im):
                  if (ext5[i,k]):
                      if (ext5[i,k-1] or ext5[i,k+1]):
                          a4[1,i,k] = a4[0,i,k]
                          a4[2,i,k] = a4[0,i,k]
                      elif (ext6[i,k-1] or ext6[i,k+1]):
                          pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                          lac_1 = pmp_1 + 1.5*gam[i,k+2]
                          a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                          pmp_2 = a4[1,i,k] + 2.*gam[i,k]
                          lac_2 = pmp_2 - 1.5*gam[i,k-1]
                          a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                  elif (ext6[i,k]):
                      if (ext5[i,k-1] or ext5[i,k+1]):
                          pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                          lac_1 = pmp_1 + 1.5*gam[i,k+2]
                          a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                          pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                          lac_2 = pmp_2 - 1.5*gam[i,k-1]
                          a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
              for i in range(0,im):
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          elif (np.abs(kord) == 12):
              for i in range(0,im):
                  if (extm[i,k]):
                      a4[1,i,k] = a4[0,i,k]
                      a4[2,i,k] = a4[0,i,k]
                      a4[3,i,k] = 0.
                  else:        #! not a local extremum
                      a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])
                      #! Check within the smooth region if subgrid profile is non-monotonic
                      if (np.abs(a4[3,i,k]) > np.abs(a4[1,i,k]-a4[2,i,k])):
                          pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                          lac_1 = pmp_1 + 1.5*gam[i,k+2]
                          a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                          pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                          lac_2 = pmp_2 - 1.5*gam[i,k-1]
                          a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
                          a4[3,i,k] = 6.*a4[0,i,k] - 3.*(a4[1,i,k]+a4[2,i,k])
          elif (np.abs(kord) == 13):
              for i in range(0,im):
                  if (ext6[i,k]):
                      if (ext6[i,k-1] and ext6[i,k+1]):
                          #! grid-scale 2-delta-z wave detected
                          a4[1,i,k] = a4[0,i,k]
                          a4[2,i,k] = a4[0,i,k]
              for i in range(0,im):
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          elif (np.abs(kord) == 14):
              for i in range(0,im):
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          elif (np.abs(kord) == 15):   #! Revised abs(kord)=9 scheme
              for i in range(0,im):
                  if (ext5[i,k] ):
                      if (ext5[i,k-1] or ext5[i,k+1]):
                          a4[1,i,k] = a4[0,i,k]
                          a4[2,i,k] = a4[0,i,k]
                  elif (ext6[i,k]):
                      #! Check within the smooth region if subgrid profile is non-monotonic
                      pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                      lac_1 = pmp_1 + 1.5*gam[i,k+2]
                      a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                      pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                      lac_2 = pmp_2 - 1.5*gam[i,k-1]
                      a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
              for i in range(0,im):
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          elif (np.abs(kord) == 16):
              for i in range(0,im):
                  if (ext5[i,k]):
                      if (ext5[i,k-1] or ext5[i,k+1]):
                          a4[1,i,k] = a4[0,i,k]
                          a4[2,i,k] = a4[0,i,k]
                      elif (ext6[i,k-1] or ext6[i,k+1]):
                          #! Left  edges
                          pmp_1 = a4[0,i,k] - 2.*gam[i,k+1]
                          lac_1 = pmp_1 + 1.5*gam[i,k+2]
                          a4[1,i,k] = np.min([np.max([a4[1,i,k], np.min([a4[0,i,k], pmp_1, lac_1])]), np.max([a4[0,i,k], pmp_1, lac_1])])
                          #! Right edges
                          pmp_2 = a4[0,i,k] + 2.*gam[i,k]
                          lac_2 = pmp_2 - 1.5*gam[i,k-1]
                          a4[2,i,k] = np.min([np.max([a4[2,i,k], np.min([a4[0,i,k], pmp_2, lac_2])]), np.max([a4[0,i,k], pmp_2, lac_2])])
              for i in range(0,im):
                  a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          else:      #! kord = 11
              for i in range(0,im):
                  if (ext5[i,k] and (ext5[i,k-1] or ext5[i,k+1])):
                      #! Noisy region:
                      a4[1,i,k] = a4[0,i,k]
                      a4[2,i,k] = a4[0,i,k]
                      a4[3,i,k] = 0.
                  else:
                      a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
                      
          #! Additional constraint to ensure positivity
          if (iv == 0):
              a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 0)

      ####end for k in range(3,km-2)

      #!----------------------------------
      #! Bottom layer subgrid constraints:
      #!----------------------------------
      if (iv == 0):
          for i in range(0,im):
              a4[2,i,km-1] = np.max([0., a4[2,i,km-1]])
      elif (iv == -1):
          for i in range(0,im): 
              if (a4[2,i,km-1]*a4[0,i,km-1] <= 0.):
                  a4[2,i,km-1] = 0.

      for k in range(km-2,km):
          for i in range(0,im):
              a4[3,i,k] = 3.*(2.*a4[0,i,k] - (a4[1,i,k]+a4[2,i,k]))
          if (k == (km-2)):
              a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 2)
          if (k == km-1):
              a4[:,:,k] = cs_limiters(im, extm[:,k], a4[:,:,k], 1)

      return a4

def mappm (km, pe1, q1, kn, pe2, i1, i2, iv, kord, ptop):
  #! IV = 0: constituents
  #! IV = 1: potential temp
  #! IV =-1: winds
   
  #! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
   
  #! pe1: pressure at layer edges (from model top to bottom surface)
  #!      in the original vertical coordinate
  #! pe2: pressure at layer edges (from model top to bottom surface)
  #!      in the new vertical coordinate

  # integer, intent(in):: i1, i2, km, kn, kord, iv
  # real, intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1) !< pe1: pressure at layer edges from model top to bottom
  #                                                        !!      surface in the ORIGINAL vertical coordinate 
  #                                                    !< pe2: pressure at layer edges from model top to bottom 
  #                                                        !!      surface in the NEW vertical coordinate
  #! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
  # real, intent(in )::  q1(i1:i2,km)
  # real, intent(out)::  q2(i1:i2,kn)
  # real, intent(IN) :: ptop
  #! local
  #        real  qs(i1:i2)
  #        real dp1(i1:i2,km)
  #    real a4(4,i1:i2,km)
  #         integer i, k, l
  #        integer k0, k1
  #        real pl, pr, tt, delp, qsum, dpsum, esl
  im = i2 - i1 + 1
  dp1 = np.zeros([im,km]) 
  a4 = np.zeros([4,im,km])
  q2 = np.zeros([im,kn])
  qs = np.zeros(im)
  
  for k in range(0,km):
      for i in range(0,im):
        dp1[i,k] = pe1[i,k+1] - pe1[i,k]
        a4[0,i,k] = q1[i,k]
  
  if ( kord > 7 ):
      a4 =  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
  else:
      a4 =  ppm_profile( a4, dp1, km, i1, i2, iv, kord )

  #!------------------------------------
  #! Lowest layer: constant distribution
  #!------------------------------------
  ##ifdef NGGPS_SUBMITTED
  #        do i=i1,i2
  #           a4(2,i,km) = q1(i,km)
  #           a4(3,i,km) = q1(i,km)
  #           a4(4,i,km) = 0.
  #    enddo
  #endif
  qsum = 0.
  for i in range(0,im):
      k0 = 0
      for k in range(0,kn):
          next_k = False
          if (pe2[i,k] <= pe1[i,0]):
  #! above old ptop
              q2[i,k] = q1[i,0]
          elif (pe2[i,k] >= pe1[i,km]):
  #! Entire grid below old ps
  ##ifdef NGGPS_SUBMITTED
  #          q2(i,k) = a4(3,i,km)   ! this is not good.
  #else
              q2[i,k] = q1[i,km-1]
  #endif
          else:
              for l in range(k0,km):
                  #! locate the top edge at pe2(i,k)
                  if ( pe2[i,k] >= pe1[i,l] and pe2[i,k] <= pe1[i,l+1]):
                      k0 = l
                      pl = (pe2[i,k]-pe1[i,l]) / dp1[i,l]
                      if (pe2[i,k+1] <= pe1[i,l+1]):
                          #! entire new grid is within the original grid
                          pr = (pe2[i,k+1]-pe1[i,l]) / dp1[i,l]
                          tt = r3*(pr*(pr+pl)+pl**2)
                          q2[i,k] = a4[1,i,l] + 0.5*(a4[3,i,l]+a4[2,i,l]-a4[1,i,l])*(pr+pl)-a4[3,i,l]*tt
                          next_k = True
                          break
                          #goto 555
                      else:
                          #! Fractional area...
                          delp = pe1[i,l+1] - pe2[i,k]
                          tt   = r3*(1.+pl*(1.+pl))
                          qsum = delp*(a4[1,i,l]+0.5*(a4[3,i,l]+a4[2,i,l]-a4[1,i,l])*(1.+pl)-a4[3,i,l]*tt)
                          dpsum = delp
                          k1 = l + 1
                          break
                          #goto 111
              if not next_k:
                  #labeled 111
                  for l in range(k1,km):
                      if( pe2[i,k+1] > pe1[i,l+1] ):
                          #! Whole layer..
                          qsum  =  qsum + dp1[i,l]*q1[i,l]
                          dpsum = dpsum + dp1[i,l]
                      else:
                          delp = pe2[i,k+1]-pe1[i,l]
                          esl  = delp / dp1[i,l]
                          qsum = qsum + delp * (a4[1,i,l]+0.5*esl*(a4[2,i,l]-a4[1,i,l]+a4[3,i,l]*(1.-r23*esl)))
                          dpsum = dpsum + delp
                          k0 = l
                          break #goto 123
                  else:  #when l-loop completes without breaking
                      delp = pe2[i,k+1] - pe1[i,km] #should this be km?
                      if (delp > 0.):
                          #! Extended below old ps
                          ##ifdef NGGPS_SUBMITTED
                          #qsum = qsum + delp * a4(3,i,km)    ! not good.
                          ##else
                          qsum = qsum + delp * q1[i,km-1] # should this be km-1?
                          ##endif
                          dpsum = dpsum + delp
                  q2[i,k] = qsum / dpsum #formerly labeled 123
    
  return q2

def map1_ppm(km, pe1, q1, qs, kn, pe2, i1, i2, iv, kord):
# subroutine map1_ppm( km,   pe1,    q1,   qs,           &
#                       kn,   pe2,    q2,   i1, i2,       &
#                       j,    ibeg, iend, jbeg, jend, iv,  kord)
#  integer, intent(in) :: i1                !< Starting longitude
#  integer, intent(in) :: i2                !< Finishing longitude
#  integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == ??? 2 == remap temp with cs scheme
#  integer, intent(in) :: kord              !< Method order
#  integer, intent(in) :: j                 !< Current latitude
#  integer, intent(in) :: ibeg, iend, jbeg, jend
#  integer, intent(in) :: km                !< Original vertical dimension
#  integer, intent(in) :: kn                !< Target vertical dimension
#  real, intent(in) ::   qs(i1:i2)       !< bottom BC
#  real, intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges from model top to bottom surface in the original vertical coordinate
#  real, intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges from model top to bottom surface in the new vertical coordinate
#  real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
# ! INPUT/OUTPUT PARAMETERS:
#  real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
# 
# ! DESCRIPTION:
# ! IV = 0: constituents
# ! pe1: pressure at layer edges (from model top to bottom surface)
# !      in the original vertical coordinate
# ! pe2: pressure at layer edges (from model top to bottom surface)
# !      in the new vertical coordinate
# 
# ! LOCAL VARIABLES:
#    real    dp1(i1:i2,km)
#    real   q4(4,i1:i2,km)
#    real    pl, pr, qsum, dp, esl
#    integer i, k, l, m, k0
# 
  im = i2 - i1 + 1
  q2 = np.zeros([im,kn])
  
  qs = np.zeros([im])
  dp1 = np.zeros([im,km])
  q4 = np.zeros([4,im,km])
  qsum = 0.
  for k in range(0,km):
      for i in range(i1-1,i2):
          dp1[i,k] = pe1[i,k+1] - pe1[i,k]
          q4[0,i,k] = q1[i,k]
# ! Compute vertical subgrid distribution
  if (kord > 7):
      q4 = cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
  else:
      q4 = ppm_profile( q4, dp1, km, i1, i2, iv, kord )

  for i in range(i1-1,i2):
      k0 = 0
      for k in range(0,kn):
          next_k = False
          for l in range(k0,km):
              # ! locate the top edge: pe2(i,k)
              if( pe2[i,k] >= pe1[i,l] and pe2[i,k] <= pe1[i,l+1] ):
                  pl = (pe2[i,k]-pe1[i,l]) / dp1[i,l]
                  if( pe2[i,k+1] <= pe1[i,l+1] ):
                      # ! entire new grid is within the original grid
                      pr = (pe2[i,k+1]-pe1[i,l]) / dp1[i,l]
                      q2[i,k] = q4[1,i,l] + 0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(pr+pl)-q4[3,i,l]*r3*(pr*(pr+pl)+pl**2)
                      k0 = l
                      next_k = True
                      #print 'new grid within old; q2 = ', q2[i,k]
                      break
                      #goto 555 #next k-loop iteration
                  else:
                      # ! Fractional area...
                      qsum = (pe1[i,l+1]-pe2[i,k])*(q4[1,i,l]+0.5*(q4[3,i,l]+q4[2,i,l]-q4[1,i,l])*(1.+pl)-q4[3,i,l]*(r3*(1.+pl*(1.+pl))))
                      for m in range(l+1,km): # was do m = l+1,km
                          # ! locate the bottom edge: pe2(i,k+1)
                          if( pe2[i,k+1] > pe1[i,m+1] ):
                              # ! Whole layer
                              qsum = qsum + dp1[i,m]*q4[0,i,m]
                          else:
                              dp = pe2[i,k+1]-pe1[i,m]
                              esl = dp / dp1[i,m]
                              qsum = qsum + dp*(q4[1,i,m]+0.5*esl*(q4[2,i,m]-q4[1,i,m]+q4[3,i,m]*(1.-r23*esl)))
                              k0 = m
                          #   goto 123
                              break
                      else:
                          #GJF: the following if statement is not in the fv_mapz, but it captures the case where pe2[kn] > pe1[km] where the m loop is not entered; without this, the lowest layer values are weird
                          if (l+1 == km):
                               dp = pe2[i,kn]-pe1[i,km]
                               esl = dp / dp1[i,km-1]
                               qsum = qsum + dp*(q4[1,i,km-1]+0.5*esl*(q4[2,i,km-1]-q4[1,i,km-1]+q4[3,i,km-1]*(1.-r23*esl)))
                          break
                      
                      break
                      #goto 123 #end l-loop
          if not next_k:
              q2[i,k] = qsum / ( pe2[i,k+1] - pe2[i,k] ) #formerly labeled 123
              
  return q2           


def fillq(im, km, nq, q, dp):
  
  for ic in range(0,nq):
      for k in range(km-1,0,-1):
          k1 = k-1
          for i in range(0,im):
              if( q[i,k,ic] < 0. ):
                  q[i,k1,ic] = q[i,k1,ic] + q[i,k,ic]*dp[i,k]/dp[i,k1]
                  q[i,k ,ic] = 0.
       
    #! Top down:
      for k in range(0,km-1):
          k1 = k+1
          for i in range(0,im):
              if( q[i,k,ic] < 0. ):
                  q[i,k1,ic] = q[i,k1,ic] + q[i,k,ic]*dp[i,k]/dp[i,k1]
                  q[i,k ,ic] = 0.

  return q

def fillz(im, km, nq, q, dp):
 #integer,  intent(in):: im                !< No. of longitudes
 #integer,  intent(in):: km                !< No. of levels
 #integer,  intent(in):: nq                !< Total number of tracers
 #real , intent(in)::  dp(im,km)           !< pressure thickness
 #real , intent(inout) :: q(im,km,nq)      !< tracer mixing ratio
 #! LOCAL VARIABLES:
 #logical:: zfix(im)
 #real ::  dm(km)
 #integer i, k, ic, k1
 #real  qup, qly, dup, dq, sum0, sum1, fac
 
 dm = np.zeros([km])
 
 #print ('orig q')
 #print q
 
 for ic in range(0,nq):
     for i in range(0,im):
         #top layer
         if( q[i,0,ic] < 0. ):
             q[i,1,ic] = q[i,1,ic] + q[i,0,ic]*dp[i,0]/dp[i,1]
             q[i,0,ic] = 0.
     #! Interior
     zfix = [False] * im
     for k in range(1,km-1):
        for i in range(0,im):
            if( q[i,k,ic] < 0. ):
                #print('neg in layer',k,q[i,k,ic])
                zfix[i] = True
                if ( q[i,k-1,ic] > 0. ):
                    #print('borrow from above')
                    #! Borrow from above
                    dq = np.min( [q[i,k-1,ic]*dp[i,k-1], -q[i,k,ic]*dp[i,k]] ) 
                    q[i,k-1,ic] = q[i,k-1,ic] - dq/dp[i,k-1]
                    q[i,k  ,ic] = q[i,k  ,ic] + dq/dp[i,k  ]
                if ( q[i,k,ic] < 0.0 and q[i,k+1,ic] > 0. ):
                    #! Borrow from below:
                    #print('borrow from below')
                    dq = np.min ( [q[i,k+1,ic]*dp[i,k+1], -q[i,k,ic]*dp[i,k]] ) 
                    q[i,k+1,ic] = q[i,k+1,ic] - dq/dp[i,k+1]
                    q[i,k  ,ic] = q[i,k  ,ic] + dq/dp[i,k  ]
                #print ('new q',q[i,k  ,ic])
     #! Bottom layer
     k = km-1
     for i in range(0,im):
        if( q[i,k,ic] < 0. and q[i,k-1,ic] > 0.):
            zfix[i] = True
            #! Borrow from above
            qup =  q[i,k-1,ic]*dp[i,k-1]
            qly = -q[i,k  ,ic]*dp[i,k  ]
            dup =  np.min([qly, qup])
            q[i,k-1,ic] = q[i,k-1,ic] - dup/dp[i,k-1] 
            q[i,k,  ic] = q[i,k,  ic] + dup/dp[i,k  ]

     #! Perform final check and non-local fix if needed
     for i in range(0,im):
        if ( zfix[i] ):
            sum0 = 0.
            for k in range(1,km):
                dm[k] = q[i,k,ic]*dp[i,k]
                sum0 = sum0 + dm[k]
            #print('sum0',sum0)
            if ( sum0 > 0. ):
                sum1 = 0.
                for k in range(1,km):
                   sum1 = sum1 + np.max([0., dm[k]])
                fac = sum0 / sum1
                #print('fac',fac)
                for k in range(1,km):
                   q[i,k,ic] = np.max([0., fac*dm[k]/dp[i,k]])
 
 return q

def mp_auto_conversion(ql, qi):
  qi0_max = 2.0E-3
  ql0_max = 2.5E-3
  qr = 0.0
  qs = 0.0  
  
  #! Convert excess cloud water into rain:
  if ( ql > ql0_max ):
      qr = ql - ql0_max
      ql = ql0_max
  #! Convert excess cloud ice into snow:
  if ( qi > qi0_max ):
      qs = qi - qi0_max
      qi = qi0_max
  
  return (ql, qr, qi, qs)

def latlon2xyz(p):

#real(kind=R_GRID), intent(in) :: p(2)
#real(kind=R_GRID), intent(out):: e(3)

#integer n
#real (f_p):: q(2)
#real (f_p):: e1, e2, e3

  e = np.zeros(3)

  e1 = math.cos(p[1]) * math.cos(p[0])
  e2 = math.cos(p[1]) * math.sin(p[0])
  e3 = math.sin(p[1])
#!-----------------------------------
#! Truncate to the desired precision:
#!-----------------------------------
  e = [e1, e2, e3]

  return e

def mid_pt3_cart(p1, p2): 
#     real(kind=R_GRID), intent(IN)  :: p1(3), p2(3)
#     real(kind=R_GRID), intent(OUT) :: e(3)
#!
#     real (f_p):: q1(3), q2(3)
#     real (f_p):: dd, e1, e2, e3
#     integer k
     
     e = np.zeros(3)
     
     # do k=1,3
     #    q1(k) = p1(k)
     #    q2(k) = p2(k)
     # enddo

     e1 = p1[0] + p2[0]
     e2 = p1[1] + p2[1]
     e3 = p1[2] + p2[2]

     dd = math.sqrt( e1**2 + e2**2 + e3**2 )
     e1 = e1 / dd
     e2 = e2 / dd
     e3 = e3 / dd

     e = [e1, e2, e3]
     
     return e

def cart_to_latlon(q):
#! vector version of cart_to_latlon1
#integer, intent(in):: np
#real(kind=R_GRID), intent(inout):: q(3,np)
#real(kind=R_GRID), intent(inout):: xs(np), ys(np)
#! local
#real(kind=R_GRID), parameter:: esl=1.d-10
#real (f_p):: p(3)
#real (f_p):: dist, lat, lon
#integer i,k
  esl = 1.0E-10
  
  dist = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)
  q = np.divide(q,dist)

  if ( (abs(q[0])+abs(q[1]))  < esl ):
    lon = 0.0
  else:
    lon = math.atan2( q[1], q[0] )   #! range [-pi,pi]
   

  if ( lon < 0.):
      lon = 2*math.pi + lon
  
  #! RIGHT_HAND system:
  lat = math.asin(q[2])

  return (lon, lat)

def mid_pt_sphere(p1, p2):
#    real(kind=R_GRID) , intent(IN)  :: p1(2), p2(2)
#    real(kind=R_GRID) , intent(OUT) :: pm(2)
#!------------------------------------------
#    real(kind=R_GRID) e1(3), e2(3), e3(3)
    
    pm = np.zeros(2)
    
    e1 = latlon2xyz(p1)
    e2 = latlon2xyz(p2)
    e3 = mid_pt3_cart(e1, e2)
    (pm[0], pm[1]) = cart_to_latlon(e3)
    
    return pm

def vect_cross(p1, p2):
    #real(kind=R_GRID), intent(in) :: p1(3), p2(3)
    #real(kind=R_GRID), intent(out):: e(3)
    e = np.zeros(3)
    
    e[0] = p1[1]*p2[2] - p1[2]*p2[1]
    e[1] = p1[2]*p2[0] - p1[0]*p2[2]
    e[2] = p1[0]*p2[1] - p1[1]*p2[0]

    return e

def normalize_vect(e):

    #real(kind=R_GRID), intent(inout):: e(3)
    #real(f_p):: pdot
    #integer k

  pdot = e[0]**2 + e[1]**2 + e[2]**2
  pdot = math.sqrt( pdot ) 
  e = e/pdot

  return e

def get_unit_vect2( e1, e2):
 #real(kind=R_GRID), intent(in) :: e1(2), e2(2)
 #real(kind=R_GRID), intent(out):: uc(3) !< unit vector e1--->e2
#! Local:
 #real(kind=R_GRID), dimension(3):: pc, p1, p2, p3
  uc = np.zeros(3)
  
#! RIGHT_HAND system:
  p1 = latlon2xyz(e1)
  p2 = latlon2xyz(e2)

  pc = mid_pt3_cart(p1, p2)
  p3 = vect_cross(p2, p1)
  uc = vect_cross(pc, p3)
  uc = normalize_vect( uc )
  
  return uc

def get_latlon_vector(pp):
    #real(kind=R_GRID), intent(IN)  :: pp(2)
    #real(kind=R_GRID), intent(OUT) :: elon(3), elat(3)
    elon = np.zeros(3)
    elat = np.zeros(3)
    
    elon[0] = -math.sin(pp[0])
    elon[1] = math.cos(pp[0])
    elon[2] =  0.0
    elat[0] = -math.sin(pp[1])*math.cos(pp[0])
    elat[1] = -math.sin(pp[1])*math.sin(pp[0])
#!!! RIGHT_HAND
    elat[2] =  math.cos(pp[1])
#! Left-hand system needed to be consistent with rest of the codes
#!  elat[2] = -math.cos(pp[1])

    return (elon, elat)

def inner_prod(v1, v2):
   #real(kind=R_GRID),intent(in):: v1(3), v2(3)
   #real (f_p) :: vp1(3), vp2(3), prod16
   #integer k

    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
