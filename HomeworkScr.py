
# coding: utf-8

# In[138]:


import numpy as np
import scipy.special as sp
import xlwt as xw
import matplotlib.pyplot as plt

numSensori = 9
detectionSensoriTeor = []
detectionSensoriSim = []
detectionSensoriSemi = []

pfa = pow(10,-2)

snrDbMin = -25
snrDbMax = -5
snrDbArray = np.arange(snrDbMin, snrDbMax + 1, 1)
snrL = pow(10,snrDbArray/10)

numProveH0 = 1000
numProveH1 = 1000
numCampioni = 1000

ps = 1
pn = ps/snrL
stdRumore = np.sqrt(pn)

segnalePu = np.sign(0.5 - np.random.rand(1, numCampioni)) + 1j*np.sign(0.5 - np.random.rand(1, numCampioni))
segnalePu = segnalePu/np.std(segnalePu)

# In[2]:


def calcoloSoglia():
    energiaH0 = []
    for i in range(numCampioni):
        energiaH0.append(0)

    energiaH0Sorted = []
    sogliaTeor = []
    sogliaSim = []
    
    for i in range(len(snrDbArray)):
        for j in range(numProveH0):
            rumore = np.random.randn(1, numCampioni) + 1j*np.random.randn(1, numCampioni)
            rumore = rumore/np.std(rumore)
            rumore = stdRumore[i] * rumore
            energiaH0[j] = np.sum(pow(np.absolute(rumore),2)) / numCampioni
        
        #Calcolo soglia via teorica
        mediaH0 = np.mean(energiaH0)
        varianzaH0 = np.var(energiaH0)
        sogliaTeor.append(mediaH0 + (pow((2*varianzaH0), 0.5) * sp.erfinv(1 - 2*pfa))) 

        #Calcolo soglia via simulata
        energiaH0Sorted = np.sort(energiaH0)
        indiceSoglia = int(numProveH0 * pfa)
        sogliaSim.append(energiaH0Sorted[ int(len(energiaH0Sorted) - indiceSoglia )])
        
    return sogliaTeor, sogliaSim


# In[83]:


def calcoloProbDetection():
    energiaH1 = []
    pdTeor = []
    pdSim = []
    pdSemi = []
    
    for i in range(numCampioni):
        energiaH1.append(0)
    
    for i in range(len(snrDbArray)):
        for j in range(numProveH1):
            rumore = np.random.randn(1, numCampioni) + 1j*np.random.randn(1, numCampioni)
            rumore = rumore/np.std(rumore)
            rumore = stdRumore[i] * rumore
            segnaleRicevuto = segnalePu + rumore
            energiaH1[j] = np.sum(pow(np.absolute(segnaleRicevuto),2)) / numCampioni
        
        #Calcolo probabilità di detection via teoria
        pdTeor.append(0.5 + 0.5 * sp.erf((-sogliaTeor[i] + np.mean(energiaH1)) * (pow((2*np.var(energiaH1)), -0.5))))
        
        #Calcolo probabilità di detection via simulata
        pdSim.append(np.sum(energiaH1 > sogliaSim[i]) / numProveH1)
        
        #Calcolo probabilità di detection via semi analitica
        pdSemi.append(np.sum(energiaH1 > sogliaTeor[i]) / numProveH1)
    
    return pdTeor, pdSim, pdSemi
        


# In[80]:


def calcoloPdSensori(sensori):
    probDiAvvenutaDetection = 0.7
    binarySensori = []
    pdAnd = []
    pdOr = []
    pdMagg = []
    for i in range(len(sensori)):
        binarySensori.append([])

    for i in range(len(sensori)):
        for j in range(len(snrDbArray)):
            if(sensori[i][j] >= probDiAvvenutaDetection):
                binarySensori[i].append(1)
            else:
                binarySensori[i].append(0)
    
    binarySensori = np.transpose(binarySensori)
    for i in range(len(snrDbArray)):
        pdAnd.append(recursiveAnd(binarySensori[i], 0, numSensori))
        pdOr.append(recursiveOr(binarySensori[i], 0, numSensori))
        pdMagg.append(iterativeMagg(binarySensori[i]))
        
    return pdAnd, pdOr, pdMagg


# In[52]:


def recursiveAnd(binarySensoriRow, i, lun):
    if(i == lun):
        return 1
    if(i < lun):
        return binarySensoriRow[i] and recursiveAnd(binarySensoriRow, i+1, lun)


# In[53]:


def recursiveOr(binarySensoriRow, i, lun):
    if i == lun:
        return 0
    if(i < lun):
        return binarySensoriRow[i] or recursiveOr(binarySensoriRow, i+1, lun)
    return


# In[54]:


def iterativeMagg(binarySensoriRow):
    magg = 0
    for i in range(numSensori):
        magg += binarySensoriRow[i]
    if magg >= len(binarySensoriRow):
        return 1
    return 0


# In[137]:


def createXls(pdTeor, pdSim, pdSemi):
    book = xw.Workbook(encoding = "utf-8")
    
    #Teorico
    f = book.add_sheet("Detection via teorica")
    f.write(0, 0, "Snr")
    f.write(0, 1, "AndTeorico")
    f.write(0, 2, "OrTeorico")
    f.write(0, 3, "MaggioranzaTeorico")
    
    for i in range(numSensori):
        f.write(0, (i+4), "PdSensoreTeorico: " + str(i+1))
    
    for i in range(len(snrDbArray)):
        f.write((i+1), 0, str(snrDbArray[i]))
        f.write((i+1), 1, str(pdAndTeor[i]))
        f.write((i+1), 2, str(pdOrTeor[i]))
        f.write((i+1), 3, str(pdMaggTeor[i]))
        for j in range(numSensori):
            f.write((i+1),(j+4), str(pdTeor[j][i]))
            
    #Simulata
    f = book.add_sheet("Detection via simulata")
    f.write(0, 0, "Snr")
    f.write(0, 1, "AndSimulato")
    f.write(0, 2, "OrSimulato")
    f.write(0, 3, "MaggioranzaSimulato")
    
    for i in range(numSensori):
        f.write(0, (i+4), "PdSensoreSimulato: " + str(i+1))
    
    for i in range(len(snrDbArray)):
        f.write((i+1), 0, str(snrDbArray[i]))
        f.write((i+1), 1, str(pdAndSim[i]))
        f.write((i+1), 2, str(pdOrSim[i]))
        f.write((i+1), 3, str(pdMaggSim[i]))
        for j in range(numSensori):
            f.write((i+1),(j+4), str(pdSim[j][i]))

    #Semianalitica
    f = book.add_sheet("Detection via semianalitica")
    f.write(0, 0, "Snr")
    f.write(0, 1, "AndSemianalitico")
    f.write(0, 2, "OrSemianalitico")
    f.write(0, 3, "MaggioranzaSemianalitico")
    
    for i in range(numSensori):
        f.write(0, (i+4), "PdSensoreSemianalitco: " + str(i+1))
    
    for i in range(len(snrDbArray)):
        f.write((i+1), 0, str(snrDbArray[i]))
        f.write((i+1), 1, str(pdAndSemi[i]))
        f.write((i+1), 2, str(pdOrSemi[i]))
        f.write((i+1), 3, str(pdMaggSemi[i]))
        for j in range(numSensori):
            f.write((i+1),(j+4), str(pdSemi[j][i]))
            
    book.save("EnergyDetection.xls")


# In[104]:


def showPlot(pdTeor, pdSim, pdSemi, num):
    #fig, axes = plt.subplots(nrows=2, ncols=2)
    plt.title('Sensore ' + str(num +1))
    plt.ylabel('Snr')
    plt.plot(snrDbArray, pdTeor, 'b', label = 'pdTeorica')
    plt.plot(snrDbArray, pdSim, 'g', label = 'pdSimulata')
    plt.plot(snrDbArray, pdSemi, 'r', label = 'pdSemianalitica')
    
    legend = plt.legend(loc='lower right', shadow=True)
    
    plt.show()


# In[134]
#Main
sogliaTeor, sogliaSim = calcoloSoglia()

for i in range(numSensori):
    pdTeor, pdSim, pdSemi = calcoloProbDetection()
    detectionSensoriTeor.append(pdTeor)
    detectionSensoriSim.append(pdSim)
    detectionSensoriSemi.append(pdSemi)
    
    showPlot(pdTeor, pdSim, pdSemi, i)
    
pdAndTeor, pdOrTeor, pdMaggTeor = calcoloPdSensori(detectionSensoriTeor)
pdAndSim, pdOrSim, pdMaggSim = calcoloPdSensori(detectionSensoriSim)
pdAndSemi, pdOrSemi, pdMaggSemi = calcoloPdSensori(detectionSensoriSemi)
createXls(detectionSensoriTeor, detectionSensoriSim, detectionSensoriSemi)
