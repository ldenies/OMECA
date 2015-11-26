# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:50:25 2015

@author: luka
"""
from __future__ import division
import numpy as np
from scipy.interpolate import interp1d

def bartz(T0,Tw,p0,M,Dt,area,visc,cp,Pr,gamma,cstar):
    s = (0.5*Tw/T0*(1 + (gamma-1)/2*M**2) + 0.5)**(-0.68) * (1 + (gamma-1)/2*M**2)**(-0.12)
    h = 0.026 / Dt**0.2 * visc**0.2 * cp / Pr**0.6 * (p0/cstar)**0.8 * (1.0/area)**0.9 * s
    return h
    
def bartzFilm(T0,Tw,p0,M,Dt,area,muFilm,cpFilm,PrFilm,gamma,cstar):
    s = (0.5*Tw/T0*(1 + (gamma-1)/2*M**2) + 0.5)**(-0.68)
    h = 0.026 / Dt**0.2 * muFilm**0.2 * cpFilm / PrFilm**0.6 * (p0/cstar)**0.8 * (1.0/area)**0.9 * s
    return h    

def adiabatic_wall(T_free,gamma,M,Pr):
    r = Pr**(1.0/3.0)
    return  T_free * (1+r* (gamma-1)/2 * M**2)

def interpol(Tv,param,Tf):
    intpFunction = interp1d(Tv[::-1],param[::-1])
    return intpFunction(Tf)

def dittusBoelter(Re,Pr):
    return 0.023 * Re**0.8 * Pr**0.4

def Taylor(Re,Pr,T,Tw,Dh,x):
    return 0.023 * Re**0.8 * Pr**0.4 * (T/Tw)**(0.57-1.59*Dh/x)
    
def Ruan(Re,Pr,rho,rhow,Dh,x):
    return 0.0069*Re**0.9*Pr**0.66*(rhow/rho)**0.43*(1+2.4*Dh/x)

def Dh_shell(Do,Di):
    return Do-Di
    
def Dh_rect(w,h):
    return 2*w*h/(w+h)
    
def colebrook(Dh, roughness, Re, f):
    return 1/(-2*np.log10(roughness/(3.7*Dh) + 2.51 / (Re*np.sqrt(f))))**2

def frictionFactor(Dh,roughness,Re):
    guess = 1e-5
    for i in range(5):
        guess = colebrook(Dh,roughness,Re,guess)
    return guess

class thermo:
    # Define constants       
    
    # General
    
    # critical temperature
    Tcrit = 190.564000000
    # critical molar density
    rho_M_crit = 10.139342719
    # Molar mass of methane (from GERG)    
    M = 16.042460
 
    # Equation of state (GERG)
    # Current and previous gas constants
    R = 8.314472
    Rstar = 8.314510   
    
    doik = [0,1,1,2,2,4,4,1,1,1,2,3,6,2,3,3,4,4,2,3,4,5,6,6,7]
    toik = [0,0.125,1.125,0.375,1.125,0.625,1.5,0.625,2.625,2.750,2.125,2,1.750,4.5,4.75,5,4,4.5,7.5,14,11.5,26,28,30,16]
    coik = [0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,6,6,6,6]
    noik = [0,0.57335704239162,-0.16760687523730e1,0.23405291834916,-0.21947376343441,0.16369201404128e-1,
    0.15004406389280e-1,0.98990489492918e-1,0.58382770929055,-0.74786867560390,0.30033302857974,0.20985543806568,
    -0.18590151133061e-1,-0.15782558339049,0.12716735220791,-0.32019743894346e-1,-0.68049729364536e-1,
    0.24291412853736e-1,0.51440451639444e-2,-0.19084949733532e-1,0.55229677241291e-2,-0.44197392976085e-2,
    0.40061416708429e-1,-0.33752085907575e-1,-0.25127658213357e-2]
    
    noik0 = [0,19.597538587,-83.959667892,3.000880000,0.763150000,0.004600000,8.744320000,-4.469210000]
    thoik = [0,0,0,0,4.306474465,0.936220902,5.577233895,5.722644361]
    
    Kpoli = 6
    Kexpi = 18  
    
    # Define viscosity constants
    d = [0.00260536,-0.0185247,0.0234216,0]

    a = [-3.12118e-5,1.99422e-7,0]
    b = [5.98858e-5,-4.91143e-5,0]
    c = [3.49668e-5,-1.73176e-5,0]
    A = [-8.52992e-10,-3.58009e-10,0]
    B = [1.60099e-8,8.50221e-10,0]
    C = [-3.55631e-7,2.80326e-7,0]    
    
    
    # Thermal conductivity constants
    a_t = [-0.5154269e-1,-0.51986885,-0.45638189e-1,0,0.14228623,0,-0.68641484e-1,0.20844530e-1]
    eps_kb = 163.558
    sigma = 0.3709
    M_t = 16.043
    f = [0,0.269386063023e1,-0.213367335398e1,0.799266772081,-0.160476700181,0.170796753887e-1,\
    -0.785578368269e-3,0.855930945735e-7,-0.145719286035e-10,0.918971501008,2000] 
    l = np.array([[0,0,0],
    [0.3365456978e1,-0.1655638071e1,0.2665846261],\
    [0.1429803995,0.1033504449,-0.173411081e-1],\
    [-0.2131339934e-1,0,0],[0.9568157622e-3,0,0],[-0.7929420721e-5,0,0]])
    
    Gamma = 0.0609
    ksi0 = 0.18e-9
    nu = 0.63
    gamma = 1.2415
    qdb = 1/(0.38e-9)
    
    k_B = 1.380658e-23
    R_cond = 1.01
    M_cond = 16.043
    
    e_over_k = 174     
    
    ###
    # Equation of state    
    ###
    
    def rhocrit(self):
        return self.rho_M_crit*self.M
    
    def tau(self,T):
        # Tcrit = 190.564000000
        return self.Tcrit/T
        
    def delta(self,rho):
        return rho/self.M/self.rho_M_crit

    def ar(self,delta,tau):
        alpha_r = sum( [ self.noik[i] * delta**self.doik[i] * tau**self.toik[i] for i in range(1,self.Kpoli+1)   ] ) + \
        sum( [ self.noik[i] * delta**self.doik[i] * tau**self.toik[i] * np.exp(-delta**self.coik[i]) for i in range(self.Kpoli+1,self.Kpoli+self.Kexpi+1)   ] )
        return alpha_r

    def a0(self,delta,tau):      
        alpha_0 = self.Rstar/self.R * (np.log(delta) + self.noik0[1] + self.noik0[2]*tau + self.noik0[3]*np.log(tau) + 
        sum ( [ self.noik0[i]*np.log(abs(np.sinh(self.thoik[i]*tau))) for i in [4,6] ] ) -
        sum ( [ self.noik0[i]*np.log(abs(np.cosh(self.thoik[i]*tau))) for i in [5,7] ] )
        )
        return alpha_0
        
    def darddelta(self,delta,tau,ddelta=0.001):
        return (self.ar(delta+ddelta,tau) - self.ar(delta-ddelta,tau))/(2*ddelta)

    def da0dtau(self,delta,tau,dtau=0.001):
        return (self.a0(delta,tau+dtau) - self.a0(delta,tau-dtau))/(2*dtau)
        
    def dardtau(self,delta,tau,dtau=0.001):
        return (self.ar(delta,tau+dtau) - self.ar(delta,tau-dtau))/(2*dtau)

    # Second derivatives
    def dar2dtau2(self,delta,tau,dtau=0.001):
        return (self.ar(delta,tau+dtau) - 2*self.ar(delta,tau) + self.ar(delta,tau-dtau))/(dtau**2)
        
    def da02dtau2(self,delta,tau,dtau=0.001):
        return (self.a0(delta,tau+dtau) - 2*self.a0(delta,tau) + self.a0(delta,tau-dtau))/(dtau**2)
        
    def dar2ddelta2(self,delta,tau,ddelta=0.001):
        return (self.ar(delta+ddelta,tau) - 2*self.ar(delta,tau) + self.ar(delta-ddelta,tau))/(ddelta**2)        

    def dar2ddeltadtau(self,delta,tau,ddelta=0.001,dtau=0.001):
        return (self.ar(delta+ddelta,tau+dtau) + self.ar(delta-ddelta,tau-dtau) \
        - self.ar(delta-ddelta,tau+dtau) - self.ar(delta+ddelta,tau-dtau)   )/(4*ddelta*dtau) 

    def pressure(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho)
        # Model is based on mol/l density, so pressure must be adjusted by factor 1000
        return (1 + delta * self.darddelta(delta,tau)) * rho/self.M*self.R*T*1000

    def cp(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho)        
        return 1000/self.M*self.R*(-tau**2*(self.da02dtau2(delta,tau)+self.dar2dtau2(delta,tau)) + \
        (1 + delta*self.darddelta(delta,tau) - delta*tau*self.dar2ddeltadtau(delta,tau))**2 / \
        (1 + 2*delta*self.darddelta(delta,tau) + delta**2 * self.dar2ddelta2(delta,tau)) )

    def cv(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho) 
        return 1000/self.M*self.R*(-tau**2*(self.da02dtau2(delta,tau)+self.dar2dtau2(delta,tau)))

    def h(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho) 
        return 1000/self.M*self.R*T*(1 + tau*(self.da0dtau(delta,tau) + self.dardtau(delta,tau)) + delta*self.darddelta(delta,tau))        

    def u(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho) 
        return 1000/self.M*self.R*T*(tau*(self.da0dtau(delta,tau) + self.dardtau(delta,tau)))     

    def dpdt(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho)
        # Pressure must be adjusted by factor 1000
        return 1000*rho/self.M*self.R*(1+delta*self.darddelta(delta,tau) - delta*tau*self.dar2ddeltadtau(delta,tau))

    def drhodp(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho)
        return self.M/(1000*self.R*T*(1 + 2*delta*self.darddelta(delta,tau) + delta**2 * self.dar2ddelta2(delta,tau)  ))

    def sonic(self,rho,T):
        tau = self.tau(T)
        delta = self.delta(rho)
        return np.sqrt(T*self.R/self.M*1000*(1 + 2*delta*self.darddelta(delta,tau) + delta**2 * self.dar2ddelta2(delta,tau)
        - (1+delta*self.darddelta(delta,tau) - delta*tau*self.dar2ddeltadtau(delta,tau))**2 / \
        (tau**2*(self.da02dtau2(delta,tau)+self.dar2dtau2(delta,tau)) ) ))
                    
            
    def eqState(self,p,T,rhoGuess=500):
        pGuess = self.pressure(rhoGuess,T)
        while abs(pGuess-p)>1:
            drdp = self.drhodp(rhoGuess,T)
            rhoGuess = rhoGuess + drdp*(p-pGuess)
            pGuess = self.pressure(rhoGuess,T)
        return rhoGuess

        
    # Friction theory viscosity  model
    # Described in Generalization of the friction theory for viscosity modeling, Quiñones-Cisneros, Sergio E. and Deiters, Ulrich K.
    # Doi: 10.1021/jp0618577   
    
    def psi1(self,tau):
        return np.exp(tau) - 1
        
    def psi2(self,tau):
        return np.exp(tau**2) - 1
        
    def Tr(self,T):
        return T/self.Tcrit

    def viscosity(self,rho,T):
        Tr = self.Tr(T)
        Gamma = 1/Tr
        
        mu0 = self.d[0]+self.d[1]*Tr**0.25+self.d[2]*Tr**0.5+self.d[3]*Tr**0.75
        psi1 = np.exp(Gamma)-1
        psi2 = np.exp(Gamma**2)-1    
        
        ka = (self.a[0] + self.a[1]*psi1+self.a[2]*psi2)*Gamma
        kr = (self.b[0] + self.b[1]*psi1+self.b[2]*psi2)*Gamma
        ki = (self.c[0] + self.c[1]*psi1+self.c[2]*psi2)*Gamma
        
        kaa = (self.A[0] + self.A[1]*psi1+self.A[2]*psi2)*Gamma**3
        krr = (self.B[0] + self.B[1]*psi1+self.B[2]*psi2)*Gamma**3
        kii = (self.C[0] + self.C[1]*psi1+self.C[2]*psi2)*Gamma**3
        
        # Find attractive, ideal and additional repulsive pressure
        (pa,pid,dpr) = self.pressures(rho,T)
        # Rewrite from mPa*s to Pa*s
        return (mu0+pa*ka+dpr*kr+pid*ki+pa**2*kaa+dpr**2*krr+pid**2*kii)/1000
        
    def pressures(self,rho,T):
        p = self.pressure(rho,T) 
        # Attractive pressure pa equals "internal pressure"
        pa = p - self.dpdt(rho,T)*T
        # Ideal pressure pid found from ideal gas theory
        pid = rho*T*self.R/self.M*1000
        dpr = p-pa-pid
        # Return pressures in bar, Quinones method is based on bar
        return pa/1e5,pid/1e5,dpr/1e5
    
#    ###
#    # Thermal conductivity
#    ###
#  
#     Model from  "The thermal conductivity of methane in the critical region"
#     By Sakonidou, E. P., van den Berg, H. R., ten Seldam, C. A. and Sengers, J. V.
#     Doi: 10.1063/1.472943
#     Critical enhancement by Olchowy and Sengers    
    
    
    def conductivity(self,rho,T):
        return self.lambda0(rho,T) + self.lambda_exc(rho,T) + self.lambda_cr(rho,T)
    
    def lambda0(self,rho,T):
        Tst = T/self.eps_kb
        return (1.39463*(T/self.M_t)**(1.0/2.0) * (self.cint_block(T)) )/(np.pi*self.sigma**2*self.G(Tst))/1e3
         
    def G(self,Tst):
        lnG = sum( [self.a_t[i]*np.log(Tst)**i for i in range(8)] )
        return np.exp(lnG)

    def cint_block(self,T):
        f = self.f        
        u =f[10]/T
        return 2.0/5.0*4.0* (sum([f[i]*T**(i/3.0) for i in range(1,7)]) + \
        sum([f[i]*T**(i-4) for i in range(7,9)]) + f[9]*(u**2*np.exp(u))/(np.exp(u)-1)**2 )
        
    def lambda_exc(self,rho,T):
        rhoM = rho/self.M_t
        Tr = T/self.Tcrit
        l = self.l
        s = 0
        for i in range(1,6):
            for j in range(3):
                s = s+l[i][j]*rhoM**i*Tr**j
        return s/1000
   
    def lambda_cr(self,rho,T):
        Rd = self.R_cond
        k_B = self.k_B
        eta = self.viscosity(rho,T)
        ksi = self.ksi(rho,T)
        cp = self.cp(rho,T)
        omega = self.omega_bar(rho,T,ksi)
        omega0 = self.omega_bar0(rho,T,ksi)
        return rho*cp*Rd*k_B*T*(omega-omega0) / (6*np.pi*eta*ksi)
       
    def ksi(self,rho,T):
        return self.ksi0 * (self.deltachistar(rho,T)/self.Gamma)**(self.nu/self.gamma)
        
    def chistar(self,rho,T):  
        pc = self.pressure(self.rhocrit(),self.Tcrit)
        drhodp = self.drhodp(rho,T)
        return rho * drhodp * pc/self.rhocrit()**2
        
    def deltachistar(self,rho,T):
        Tr = 2 * self.Tcrit
        deltachi = self.chistar(rho,T) - self.chistar(rho,Tr)*Tr/T
        if deltachi>0:
            return deltachi
        else:
            # Give finite but non-zero value if temperature is above reference temperature
            return 1e-15
    
    def omega_bar(self,rho,T,ksi):
        cp = self.cp(rho,T)
        cv = self.cv(rho,T)
        return 2.0/np.pi * (((cp-cv)/cp) * np.arctan(self.qdb*ksi) + cv/cp*self.qdb*ksi)
        
    def omega_bar0(self,rho,T,ksi):
        return 2.0/np.pi*(1-np.exp(-1/((self.qdb*ksi)**(-1) + \
        (self.qdb*ksi*self.rhocrit()/rho)**2/3)))