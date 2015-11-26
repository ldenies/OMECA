# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:59:24 2015

@author: luka
"""

# Initialize and import packages
import numpy as np
import CEA
import thermoClass as th
import scipy.interpolate as intp
import matplotlib.pyplot as plt

# Calculate fuel flow from engine parameters
mflow = 3.176 # kg/s
OFratio = 3.16 
fFlow = 1/(1+OFratio)*mflow

# Define nozzle material and thickness
tChamber = 4.2e-3 # wall thickness
kChamber = 295 # W/(m2 K)
rhoChamber = 9134
mu = 0.34 # Poisson's ratio
E = 85e9  # Modulus of elasticity
s_yield = 120932000 # Yield strength

# Define channel geometry
NChannels = 64
tRib = 1e-3
channelHeight = 1e-3
roughness = 6e-6    

# Initialize coolant pressure and temperature
p = pin = 60e5
T = Tin = 110

# Read nozzle coordinates
cont = np.genfromtxt("nozzleContour.csv",delimiter=",")

# Define function for radius of curvature based on coordinates of 3 points
def radiusCurvature(x1,y1,x2,y2,x3,y3):
    num = np.sqrt(((x2-x1)**2+(y2-y1)**2)*((x2-x3)**2+(y2-y3)**2)*((x3-x1)**2+(y3-y1)**2))
    den = 2*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2)
    return num/den

xVals = cont[0,::]
yVals = cont[1,::]
# Define engine size (throat radius and area)
rt = 0.0216
At = rt**2*np.pi
aRatioMinm = min(yVals**2/rt**2)


def interpol(x,y,xNew,how="linear"):
    f = intp.interp1d(x,y,kind=how)
    return f(xNew)

xHeight = np.array([0,  9,  11, 13, 15, 16, 18, 20, 30])*1e-2
Height = np.array([ 0.8,0.8,0.6,1.0,3.0,1.0,0.4,1.1,2])*1e-3

# Check for inward buckling (due to coolant pressure)
l = max(xVals)
r = max(yVals)
t = tChamber

gamma = 1
pcrit = 0.855 * E * np.sqrt(gamma) / ( (1-mu**2)**(3./4.) * (r/t)**(5./2.) * (l/r))
if pcrit>pin:
    print "Buckling pressure okay:",pcrit/1e5,"bar"
else:
    print "Buckling pressure exceeded:",pcrit/1e5,"bar"
    
# Check for hoop stress (due to chamber pressure)    
s_h = pin*r/t
if s_h<s_yield:
    print "Hoop stress okay",s_h/1e6,"MPa"
else:
    print "Hoop stress exceeded",s_h/1e6,"MPa"

# Read CEA file to find adiabatic wall temperature and convective coefficient
CEAfile = "dataCea/methalox.out"
AreaCEA,pCEA,TCEA,rhoCEA,MCEA,muCEA,cpCEA,PrCEA,gCEA,pH2O,pCO2,cstar = CEA.read(CEAfile)
T0 = TCEA[0]
p0 = pCEA[0]
# Create class with methane thermophysical model
methane = th.thermo()

# Initialize empty lists for parameter storage
pvals = []
Tvals = []
p0vals = []
T0vals = []
rhovals = []
Twvals = []
hcvals = []
hgvals = []
wvals = []
Revals = []
Nuvals = []
qvals = []
qradvals = []
Tawvals = []
Civals = []
Dhvals = []
Prvals = []
muvals = []
kapvals = []
cpvals = []
channelHeightvals = []
# Set channel wall temperature to coolant inlet temperature for 1st channel
TwChannel = T
# Pointer to indicate what the current CEA station is, start at nozzle end
CEAval_curr = len(AreaCEA)-1
# Start channel length at nonzero value to make sure Taylor equation does not crash
x = 0.01
Q = 0
Atot = 0
V = 0
rho = methane.eqState(p,T)
cp = methane.cp(rho,T)
Tw = 400
mTot = 0
# Start calculation loop from end of nozzle towards combustion chamber
for i in range(1,len(xVals)):

    # Calculate lenght of channel part and geometry of chamber
    l = np.sqrt( (xVals[-i-1]-xVals[-i])**2 + (yVals[-i-1]-yVals[-i])**2   )
    Rnozzle = yVals[-i]
    aRatio = yVals[-i]**2/rt**2
    
    channelHeight = interpol(xHeight,Height,xVals[-i])    
    channelHeightvals.append(channelHeight)
    # Calculate channel cross-sectional dimensions at this nozzle station
    if NChannels==1:
        A = np.pi*( (Rnozzle+channelHeight)**2 - Rnozzle**2 )      
        Dh = th.Dh_shell(Rnozzle+channelHeight,Rnozzle)
    else:
        channelWidth = Rnozzle*2*np.pi/NChannels - tRib
        if channelWidth<0:
            print "Error: channel width smaller than 0"
        A = NChannels*channelWidth*channelHeight
        Dh = th.Dh_rect(channelWidth,channelHeight)

    # Calculate dynamic pressure and temperature at previous station
    dynPres1 = 0.5*rho*V**2
    dynTemp1 = 0.5*V**2/cp

    # Calculate density and flow velocity
    rho = methane.eqState(p,T,rho)
    V = fFlow / (A * rho)
    
    # Calculate/update static pressure and temperature
    dynPres2 = 0.5*rho*V**2
    p = p - (dynPres2-dynPres1)    
    
    dynTemp2 = 0.5*V**2/cp
    T = T - (dynTemp2 - dynTemp1)
    
    # Calculate thermodynamic properties of methane at current (rho,T)
    mu = methane.viscosity(rho,T)
    cp = methane.cp(rho,T)
    gam = cp/methane.cv(rho,T)
    kap = methane.conductivity(rho,T)
    # Calculate bulk flow properties of coolant

    Re = V*rho*Dh/mu
    Pr = mu*cp/kap
    
    # Correct for curvature of channel alongside nozzle    
    if i>1 and i<len(xVals):
        (x1,y1) = (xVals[-i-1],yVals[-i-1])
        (x2,y2) = (xVals[-i],yVals[-i])
        (x3,y3) = (xVals[-i+1],yVals[-i+1])
        Rc = radiusCurvature(x1,y1,x2,y2,x3,y3)
        # Use Niino's formula
        Ci = (Re*(Dh/4/abs(Rc))**2)**(np.sign(Rc)*0.05)
    # If radius is too high, set correction to 1 (no correction)
        if abs(Rc) > 1:
            Ci = 1
            Rc = 1e9
    else:
        Ci = 1
        Rc = 1e9
    ksi = th.frictionFactor(Dh,roughness,Re)/th.frictionFactor(Dh,0,Re)      
    Cksi = (1+1.5*Pr**(-1./6.)*Re**(-1./8.)*(Pr-1))*ksi/(1+1.5*Pr**(-1./6.)*Re**(-1./8.)*(Pr*ksi-1)) 
    
    # Check if CEA station should be shifted, depending on current area ratio
    if aRatio >= AreaCEA[CEAval_curr] and aRatio >= AreaCEA[CEAval_curr-1] and CEAval_curr>1:    
        CEAval_curr = CEAval_curr - 1
    elif aRatio <= AreaCEA[CEAval_curr] and aRatio <= AreaCEA[CEAval_curr-1] and CEAval_curr>1:     
        CEAval_curr = CEAval_curr - 1 
    elif abs(aRatio-aRatioMinm)<1e-6:
        CEAval_curr = CEAval_curr - 1
    
    # Calculate hot gas parameters depending on CEA values
    pWater = CEA.interpol(aRatio,AreaCEA,CEAval_curr,pH2O)
    pCarbDiox = CEA.interpol(aRatio,AreaCEA,CEAval_curr,pCO2)
    Tg = CEA.interpol(aRatio,AreaCEA,CEAval_curr,TCEA)
    Mg = CEA.interpol(aRatio,AreaCEA,CEAval_curr,MCEA)
    gg = CEA.interpol(aRatio,AreaCEA,CEAval_curr,gCEA)
    Prg = CEA.interpol(aRatio,AreaCEA,CEAval_curr,PrCEA)    
    cpg = CEA.interpol(aRatio,AreaCEA,CEAval_curr,cpCEA)
    mug = CEA.interpol(aRatio,AreaCEA,CEAval_curr,muCEA)    
    Taw = th.adiabatic_wall(Tg,gg,Mg,Prg)
    
    # Increase TwNew to avoid missing loop
    TwNew = Tw+10
    TwChannelNew = TwChannel+10

    while (abs(TwNew-Tw)>0.1) and (abs(TwChannel-TwChannelNew)>0.1):
        Tw = TwNew  
        TwChannel = TwChannelNew             
        # Calculate convective coefficient using Bartz
        hg = th.bartz(T0,Tw,p0,Mg,rt*2,aRatio,mug,cpg,Prg,gg,cstar)
        hg = hg/0.026*0.0195

        # Calculate Nusselt number
        Nu = th.Taylor(Re,Pr,T,TwChannel,Dh,x)
        #Nu = th.dittusBoelter(Re,Pr)
        #rhow = methane.eqState(p,TwChannel)
        #Nu = th.Ruan(Re,Pr,rho,rhow,Dh,x)
        # Apply correction to Nusselt number
        Nu = Nu*Ci*Cksi
        # Calculate coolant convective coefficient
        hc = Nu*kap/Dh
        
        # Calculate fin effectiveness        
        m = np.sqrt(2*hc*tRib/kChamber)
        finEffectiveness = np.tanh(m/tRib*channelHeight)/(m/tRib*channelHeight)
        hc = hc*(channelWidth + finEffectiveness*2*channelHeight) / (channelWidth+tRib)
        
        # Calculate radiative heat transfer
        qW = 5.74 * (pWater/1e5*Rnozzle)**0.3 * (Taw/100)**3.5
        qC = 4 * (pCarbDiox/1e5*Rnozzle)**0.3 * (Taw/100)**3.5
        qRad = qW+qC   
        
        # Calculate heat flux
        q = (Taw-T+qRad/hg) / (1/hg + tChamber/kChamber + 1/hc)
        # Calculate hot gas wall temperature and channel wall temperature    
        TwNew = Taw-(q-qRad)/hg
        TwChannelNew = T+q/hc
        
    Tw = TwNew
    TwChannel = TwChannelNew
    

    # Calculate change in temperature and pressure
    A_heat = 2*np.pi*Rnozzle*l
    deltaT = q*A_heat/(fFlow*cp)
    deltap = th.frictionFactor(Dh,roughness,Re) * l/Dh * rho*V**2 / 2.0    
    Q = Q+q*A_heat   
    Atot = Atot+A_heat
    mCur = (2*np.pi*Rnozzle*l*tChamber+l*tRib*channelHeight*NChannels)*rhoChamber
    mTot = mTot + mCur
    # Update pressure, temperature and channel length
    p = p - deltap
    T = T + deltaT    
    x = x + l
    
    p0vals.append(p+0.5*rho*V**2)
    T0vals.append(T+0.5*V**2/cp) 

    # Store parameters in lists
    pvals.append(p)
    Tvals.append(T)
    rhovals.append(rho)
    Twvals.append(Tw)
    Tawvals.append(Taw)
    hcvals.append(hc)
    hgvals.append(hg)
    wvals.append(channelWidth)
    Revals.append(Re)
    Nuvals.append(Nu)
    qvals.append(q)
    Civals.append(Ci)
    Dhvals.append(Dh)
    Prvals.append(Pr)
    muvals.append(mu)
    kapvals.append(kap)
    cpvals.append(cp)
    
# Print output for user
print min(wvals)*1e3, "mm minimum channel width"
print max(Twvals), "K maximum wall temperature"
print (p0vals[0]-p0vals[-1])/1e5,"bar pressure loss"
print T-Tvals[0], "K temperature rise"
print Q, "Total heat input"
print mTot, "kg chamber mass"
# Plot results


# Create figure
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(15/2.54,6/2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = range(4)

# Wall temperature
lins[0] = ax.plot(xVals[1:]*100, Twvals[::-1], 'g--', lw=2, label = r'$T_w$')
ax.set_ylim([0,round(max(Twvals)+100,2)])

# Heat flux
ax2 = ax.twinx()
lins[1] = ax2.plot(xVals[1:]*100, np.array(qvals[::-1])/1e7, 'r-.',lw=2,  label = r'$q$')

# Geometry
heights = interpol(xHeight,Height,xVals)   
lins[2] = ax2.plot(xVals*100,heights*1e3, 'b:',lw=2,label = r'$d_c$')
lins[3] = ax2.plot(xVals*100, yVals*10, 'k-', label = r'Contour')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines,labs,loc=6,labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Temperature [K]")
ax2.set_ylabel(r"$d_c$ [mm]; $q$ [$\mathrm{10^7 W/m^2}$]; Radius [10 cm]")
ax.set_ylim([400,800])
ax2.set_ylim([0,4])
ax.grid()
plt.show()

# Create figure
fig = plt.figure(2)
fig.clf()
fig.set_size_inches(15/2.54,6/2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = range(3)

# Reynolds number
lins[0] = ax.plot(xVals[1:]*100, np.array(Revals[::-1])/1e4, 'r--', lw=2, label = r'Re')
#ax.set_ylim([0,round(max(Revals)+100,2)])

# Nusselt number
ax2 = ax.twinx()
lins[1] = ax2.plot(xVals[1:]*100, Nuvals[::-1], 'b-.',lw=2,  label = r'Nu')

# Nozzle contour
lins[2] = ax.plot(xVals*100, yVals*100, 'k-', label = r'Contour')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines,labs,loc=0,labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Radius [cm]; Re [$\mathrm{10^4}$]")
ax2.set_ylabel(r"Nu [-]")
ax.set_ylim([0,80])
ax2.set_ylim([0,4000])
ax.grid()
plt.show()


# Create figure
fig = plt.figure(3)
fig.clf()
fig.set_size_inches(15/2.54,6/2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = range(2)

# Pressure
lins[0] = ax.plot(xVals[1:]*100, np.array(p0vals[::-1])/1e5, 'b--', lw=2, label = r'$p_0$')

# Temperature
ax2 = ax.twinx()
lins[1] = ax2.plot(xVals[1:]*100, T0vals[::-1], 'r-.',lw=2,  label = r'$T_0$')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines,labs,loc=7,labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"$p_{0,c,b}$ [bar]")
ax2.set_ylabel(r"$T_{0,c,b}$ [K]")
ax.set_ylim([53,61])
ax2.set_ylim([100,500])
ax.grid()
plt.show()

# Create figure
fig = plt.figure(4)
fig.clf()
fig.set_size_inches(15/2.54,6/2.54)
ax = fig.add_subplot(111)

# Create four plots
lins = range(1)

# Pressure
lins[0] = ax.plot(xVals*100, yVals*100, 'k-', label = r'Contour')
ax.plot(xVals*100, -yVals*100, 'k-')

# Create legend
labs = [line[0].get_label() for line in lins]
lines = [line[0] for line in lins]
ax.legend(lines,labs,loc=7,labelspacing=0)

# Create labels, show and save
ax.set_xlabel(r"$x$ coordinate [cm]")
ax.set_ylabel(r"Radius [cm]")
ax.grid()
plt.show()
