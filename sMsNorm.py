import os
import numpy as np
from scipy.interpolate import interp1d
import sys
sys.path.append('./baseTools')
from plotClass import plotCLASS


class sMsNormCalc():

  def __init__(self, deBrogWL, scatAmpFolder, NpadBins):

    self.plc          = plotCLASS()
    self.atomTypes    = ["carbon", "nitrogen", "oxygen", "iodine", "flourine", "hydrogen"]
    self.NpadBins     = NpadBins
    self.scatAmpIntrp = {}

    for atm in self.atomTypes:      
      angStr = []
      sctStr = []
      with open(scatAmpFolder + atm + "_dcs.dat", 'r') as inpFile:
        ind=0
        for line in inpFile:
          if ind < 31:
            ind += 1
            continue

          angStr.append(line[2:11])
          sctStr.append(line[39:50])

      angs = np.array(angStr).astype(np.float64)
      S    = 4*np.pi*np.sin((angs/2)*(np.pi/180))/deBrogWL
      scts = np.sqrt(np.array(sctStr).astype(np.float64))

      self.scatAmpIntrp[atm] = interp1d(S, scts, 'cubic')


  def analyze_sMs(self, molParams, Nbins, maxS):

    if "atoms" not in molParams or "Natoms" not in molParams:
      raise RuntimeError("Cannot find key 'atoms' or 'Natoms'!!!")

    QperPix       = float(maxS)/float(Nbins)
    atomicScat    = np.zeros(Nbins)
    self.Seval    = np.linspace(0, maxS, Nbins)
    scatAmps      = []
    bondScatAmps  = []
    bondSMSscale  = []
    bondTypes     = []
    for i,atm in enumerate(molParams["atoms"]):
      scatAmps.append(self.scatAmpIntrp[atm](self.Seval))
      atomicScat += molParams["Natoms"][i]\
                      *(self.scatAmpIntrp[atm](self.Seval))**2
    
      if atm == "hydrogen":
        continue

      for j in range(i, len(molParams["atoms"])):
        if (i == j) and (molParams["Natoms"][i] == 1):
          continue
        if molParams["atoms"][j] == "hydrogen":
          continue

        bondTypes.append(atm + "-" + molParams["atoms"][j])
        bondScatAmps.append(self.scatAmpIntrp[atm](self.Seval)*
            self.scatAmpIntrp[molParams["atoms"][j]](self.Seval))

    #Normalize by atomic scattering
    for i in range(len(bondScatAmps)):
      bondSMSscale.append(bondScatAmps[i]/atomicScat)


    # Plotting
    opts = {
        "labels" : bondTypes,
        "xTitle" : r"Q [$\AA^{-1}$]"}

    opts["legOpts"] = {
        "loc"             : 'upper left'}
    self.plc.print1d(
        np.array(bondSMSscale), 
        "./plots/" + molParams["molName"] + "_bondsMsScale",
        xRange=[0,maxS],
        isFile=False,
        options=opts)
    del opts["legOpts"] 

    opts["yLog"] = True
    self.plc.print1d(
        np.array(bondScatAmps), 
        "./plots/" + molParams["molName"] + "_bondScatAmps",
        xRange=[0,maxS],
        isFile=False,
        options=opts)

    opts["labels"] = molParams["atoms"]
    self.plc.print1d(
        np.array(scatAmps), 
        "./plots/" + molParams["molName"] + "_inputScatAmps",
        xRange=[0,maxS],
        isFile=False,
        options=opts)


    ##### Expected Diffraction Signal  #####
    sinusiod = np.sin(8*self.Seval*2*np.pi/maxS)

    sMsSignal = [sinusiod/sum(molParams["Natoms"])]
    for i in range(len(bondTypes)):
      sMsSignal.append(bondSMSscale[i]*sinusiod)

    del opts["yLog"]
    opts["labels"] = ["Single Atom Molecule"] + bondTypes
    opts["legOpts"] = {
        "loc"             : 'upper center',
        "bbox_to_anchor"  : (0.5, 1.14),
        "ncol"            : 2}

    self.plc.print1d(
        np.array(sMsSignal), 
        "./plots/" + molParams["molName"] + "_sMsDiffSignal",
        xRange=[0,maxS],
        isFile=False,
        options=opts)
    del opts["legOpts"]


    #####  Expected Pair Correlation Signal  #####
    
    # Sanity Check
    sinusiod.astype(np.double).tofile(
        "./results/sanityCheck_" 
        + molParams["molName"]
        + "[" + str(Nbins) + "].dat")

    os.system("./pairCorr.exe "
        "./results/sanityCheck_" 
        + molParams["molName"]
        + "[" + str(Nbins) + "].dat "
        + " -QpPix " + str(QperPix)
        + " -FiltSTD " + str(molParams["filterSTD"])
        + " -Npad 0"
        + " -doFilt false")

    # Signal for each bond type
    sinCoeff = []
    cosCoeff = []
    curList = ["homoNuclear"] + bondTypes
    for i in range(len(sMsSignal)):
      sMsSignal[i].astype(np.double).tofile(
          "./results/" + molParams["molName"]
          + "_" + curList[i] 
          + "[" + str(Nbins) + "].dat")

      os.system("./pairCorr.exe "
          "./results/" + molParams["molName"]
          + "_" + curList[i] + "[" + str(Nbins) + "].dat "
          + " -QpPix " + str(QperPix)
          + " -FiltSTD " + str(molParams["filterSTD"])
          + " -Npad " + str(self.NpadBins))

      sinCoeff.append(np.fromfile("./results/pairCorrOdd_"
          + molParams["molName"]
          + "_" + curList[i] + ".dat"))

      cosCoeff.append(np.fromfile("./results/pairCorrEven_"
          + molParams["molName"]
          + "_" + curList[i] + ".dat"))

    opts["xTitle"] = r"R [$\AA$]"
    maxR = ((Nbins+self.NpadBins)/2 + 1)\
              *(2*np.pi/((float(maxS)/float(Nbins))\
              *(Nbins+self.NpadBins)))*0.007
    self.plc.print1d(
        np.array(sinCoeff), 
        "./plots/" + molParams["molName"] + "_sinPairCorrSignal",
        xRange=[0,maxR],
        isFile=False,
        options=opts)

    self.plc.print1d(
        np.array(cosCoeff), 
        "./plots/" + molParams["molName"] + "_cosPairCorrSignal",
        xRange=[0,maxR],
        isFile=False,
        options=opts)











if __name__ == "__main__":

  scatAmpFile = "./scatteringAmplitudes/"
  deBrog      = 0.002966 #Ang
  NpadBins    = 0

  sms = sMsNormCalc(deBrog, scatAmpFile, NpadBins)

  nbzParams = {
    "molName"   : "nitrobenzene",
    "atoms"     : ["carbon", "nitrogen", "oxygen", "hydrogen"],
    "Natoms"    : [6, 1, 2, 5],
    "filterSTD" : 3.25}

  sms.analyze_sMs(nbzParams, 5000, 12)


  CF3IParams = {
    "molName" : "CF3I",
    "atoms"   : ["carbon", "flourine", "iodine"],
    "Natoms"  : [1, 3, 1],
    "filterSTD" : 3.25}

  sms.analyze_sMs(CF3IParams, 5000, 12)

  CHDParams = {
    "molName" : "CHD",
    "atoms"   : ["carbon", "hydrogen"],
    "Natoms"  : [6, 8],
    "filterSTD" : 3.25}

  sms.analyze_sMs(CHDParams, 5000, 12)
