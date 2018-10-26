

This package will calculate the atomic scattering and plot the Q dependent 
variations in the modified molecular scattering caused by dividing the 
scattering amplitude contribution from a single bond by the total atomic 
scattering of the molecule. This package will additionally plot the expected
signal for each bond both in recipricol and real space. In addition to the 
sine transformation, the cosine transformation is shown as well to show 
how much power is transferred from the sine projection to the cosine 
projection through the modified molecular scattering Q dependence and by 
applying a smoothing filter to send the high Q contributions to 0 by 
weighting with a gaussian centered at Q=0.


Setup:
  mkdir plots results

  If executable (pairCorr.exe) does not work when called by sMsNorm.py
  then modify Makefile with the correct library and include addresses 
  and remake pairCorr.exe via the command below:
    make clean; make
    

How to run:

  python sMsNorm.py

    code that runs:

      sms = sMsNormCalc(deBrog, scatAmpFolder, NpadBins)

      nbzParams = {
        "molName"   : "nitrobenzene",
        "atoms"     : ["carbon", "nitrogen", "oxygen", "hydrogen"],
        "Natoms"    : [6, 1, 2, 5],
        "filterSTD" : 3.25}

      sms.analyze_sMs(nbzParams, Nbins, maxS)


        Variables:
        
          deBrog:         DeBroglie wavelength of the electron in Angstroms
          scatAmpFolder:  Folder path to the scattering amplitudes
          NpadBins:       Number of bins to pad the FFT to increase granularity in the pair correlation
          Nbins:          Number of bins in the diffraction line out
          maxS:           How far in Q the diffraction pattern will go.

          molName:        Molecule name
          atoms:          List of all atom species in the molecule, scatttering amplitudes 
                              for these atoms must exist in ./scatteringAmplitudes
          Natoms:         The multiplicity of each atom in 'atoms'.
          filterSTD:      The standard deviation of the smoothing gaussian is defined as
                              Nbins/filterSTD


pairCorr.cpp

  This program claculates the sine and cosine transform of the guassian filtered
        input. pairCorr.exe is called by sMsNorm.py to do the sine and cosine 
        tranformation of the modified molecular diffraction line out saved in 
        results. The coefficients from the transforms are saved in results as 
        well and then plotted by sMsNorm.py.


baseTools

  Directory that contains various classes and namespaces used to save, plot,
        and analyze data. 
