import numpy as np

def read(fName):
    data = open(fName)
    content = data.readlines()
    data.close()
    
    Aratio_comb = 12
    
    area = np.array([])
    for line in content:
        if line[0:6] == " Ae/At":
            area = np.append(area,Aratio_comb)
            for val in line[17:].split():
                area = np.append(area,float(val))
                 
    for line in content:
        if line[0:6] == " CSTAR":
            cstar = float(line[17:].split()[0])
            break
                        
    gamma = np.array([])
    for line in content:
        if line[0:7] == " GAMMAs":
            for val in line[17:].split():
                gamma = np.append(gamma,float(val))            
    
    temper = np.array([])
    for line in content:
        if line[0:5] == " T, K":
            for val in line[17:].split():
                temper = np.append(temper,float(val))
                
    visc = np.array([])
    for line in content:
        if line[0:5] == " VISC":
            for val in line[17:].split():
                visc = np.append(visc,float(val)/10000)
    
    pressure = np.array([])
    for line in content:
        if line[0:7] == " P, BAR":
            for val in line[17:].split():
                pressure = np.append(pressure,float(val)*100000)
                
    rho = np.array([])
    for line in content:
        if line[0:4] == " RHO":
            mag = "0"
            for val in line[17:].split():
                if val=="0":
                    pass
                elif len(val)==1:
                    mag = val
                else:
                    if val[-2]=="-":
                        val = val[:-2] + "e" + val[-2:]
                    else:
                        val = val + "e" + mag
                        mag = "0"
                    rho = np.append(rho,float(val))                
                
    nH2O = np.array([])
    for line in content:
        if line[0:7] == " H2O   ":
            for val in line[17:].split():
                nH2O = np.append(nH2O,float(val))
                
    nCO2 = np.array([])
    for line in content:
        if line[0:8] == " *CO2   ":
            for val in line[17:].split():
                nCO2 = np.append(nCO2,float(val))
                                
    mach = np.array([])
    switch_frozen = 0
    for line in content:
        if line[0:69] == "           THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION" and switch_frozen == 0:
            switch_frozen = len(mach)-1
        if line[0:5] == " MACH":
            for val in line[17:].split():
                mach = np.append(mach,float(val))            
    
    cp = np.array([])
    switch_write = 1
    for line in content:
        if line[:28] == "  WITH EQUILIBRIUM REACTIONS":
            switch_write = 0
        if line[:23] == "  WITH FROZEN REACTIONS":
            switch_write = 1
        if line[:7] == " Pinf/P":
            switch_write = 0
        if line[0:15] == " Cp, KJ/(KG)(K)" and switch_write == 1:
            for val in line[17:].split():
                cp = np.append(cp,float(val)*1000)
    
    prandtl = np.array([])
    freeze = 0
    switch_write  = 0
    for line in content:
        if line[0:69] == "           THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION":
            freeze = 1
            switch_write = 1
        if line[:28] == "  WITH EQUILIBRIUM REACTIONS" and freeze == 0:
            switch_write = 0
        if line[:23] == "  WITH FROZEN REACTIONS" and freeze == 0:   
            switch_write = 1
        if line[0:8] == " PRANDTL" and switch_write==1:
            for val in line[17:].split():
                prandtl = np.append(prandtl,float(val))                  
    
    sub = mach<=1
    sup = mach>1
    frozen = np.array(range(len(area))) > switch_frozen
    
    nH2O = np.concatenate((nH2O,nH2O[1]*np.ones(len(pressure)-len(nH2O))))
    nCO2 = np.concatenate((nCO2,nCO2[1]*np.ones(len(pressure)-len(nCO2))))
    
    area = np.concatenate((area[sub], area[sup&frozen]))
    pressure = np.concatenate((pressure[sub], pressure[sup&frozen]))
    rho = np.concatenate((rho[sub], rho[sup&frozen]))
    cp = np.concatenate((cp[sub], cp[sup&frozen]))
    gamma = np.concatenate((gamma[sub], gamma[sup&frozen]))
    temper = np.concatenate((temper[sub], temper[sup&frozen]))
    visc = np.concatenate((visc[sub], visc[sup&frozen]))
    mach = np.concatenate((mach[sub], mach[sup&frozen]))
    prandtl = np.concatenate((prandtl[sub], prandtl[sup&frozen]))
    nH2O = np.concatenate((nH2O[sub], nH2O[sup&frozen]))
    nCO2 = np.concatenate((nCO2[sub], nCO2[sup&frozen]))

    pressure,ind_select = np.unique(pressure,return_index=True)
    pressure = pressure[::-1]
    rho = rho[ind_select[::-1]]
    area = area[ind_select[::-1]]
    cp = cp[ind_select[::-1]]
    gamma = gamma[ind_select[::-1]]
    temper = temper[ind_select[::-1]]
    visc = visc[ind_select[::-1]]
    mach = mach[ind_select[::-1]]
    prandtl = prandtl[ind_select[::-1]]
    nH2O = nH2O[ind_select[::-1]]
    nCO2 = nCO2[ind_select[::-1]]
    
    pH2O = nH2O*pressure
    pCO2 = nCO2*pressure

    return area,pressure,temper,rho,mach,visc,cp,prandtl,gamma,pH2O,pCO2,cstar
    
def interpol(aRatio,area_CEA,CEAval_curr,param_CEA):
    dist = abs(aRatio-area_CEA[CEAval_curr])/abs(area_CEA[CEAval_curr]-area_CEA[CEAval_curr-1])
    return (1-dist)*param_CEA[CEAval_curr]+dist*param_CEA[CEAval_curr-1]