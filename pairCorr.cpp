//Home Grown
#include "baseTools/tools.h"
#include "baseTools/plotClass.h"
#include "baseTools/saving.h"

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'runName' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }


  cout.setf(ios::scientific);
  
  //////////////////////////////////
  /////  Setting up variables  /////
  //////////////////////////////////

  int NQbins    = -1;
  int NPadBins  = -1;
  float QperPix = -1;
  float filtSTD = -1;
  float rMaxRat = 0.007;
  bool doFiltering = true;
  bool verbose  = true;
  bool pltVerbose  = true;
  PLOTclass plt;

  /////  Importing variables from command line  /////
  if (verbose)
    std::cout << "Importing command line variables.\n";

  if (argc < 2) {
    cerr << "Must give name of file.\n";
    exit(0);
  }

  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg], "-QpPix") == 0) {
      QperPix = atof(argv[iarg+1]);
    }
    else if (strcmp(argv[iarg], "-doFilt") == 0) {
      string str(argv[iarg+1]);
      if (str.compare("false") == 0) doFiltering = false;
    }
    else if (strcmp(argv[iarg], "-FiltSTD") == 0) {
      filtSTD = atof(argv[iarg+1]);
    }
    else if (strcmp(argv[iarg], "-Npad") == 0) {
      NPadBins = atoi(argv[iarg+1]);
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }


  /////  Getting Parameters/names from File Names  ///// 
  std::string inputFileName(argv[1]);
  auto sInd = inputFileName.find("results");
  std::string folderName  = inputFileName.substr(0, sInd+7);
  cout<<"folder: "<<folderName<<endl;
  std::string fileName    = inputFileName.substr(
                            sInd+8, 
                            inputFileName.length()-sInd-8);
  cout<<"filename: "<<fileName<<endl;
  std::string molName     = inputFileName.substr(sInd+8, 
                              inputFileName.find("_")-sInd-8);
  cout<<"molecule: "<<molName<<endl;
  sInd = inputFileName.find("_");
  std::string bondType    = inputFileName.substr(sInd+1, 
                              inputFileName.find("[")-sInd-1);

  cout<<"bond: "<<bondType<<endl;


  /////  Plotting variables  /////
  std::vector<PLOToptions> opts(2);
  std::vector<std::string> vals(2);
  std::vector<PLOToptions> oppts(3);
  std::vector<std::string> vaals(3);
  oppts[0] = yLabel;   vaals[0] = "Time [ps]";
  oppts[1] = xLabel;   vaals[1] = "Scattering Q [arb units]";
  oppts[2] = draw;     vaals[2] = "COLZ";
  opts[0] = xLabel;   vals[0] = "R [#AA]";
  opts[1] = xSpan;    


  ////////////////////////////
  /////  Importing data  /////
  ////////////////////////////

  if (verbose)
    std::cout << "Importing data.\n";

  std::vector<int> shape = save::getShape(folderName, fileName);
  std::vector<double> sMsDiff;
  sMsDiff.resize(shape[0]);

  if (verbose)
    std::cout << "Retriving data from " << fileName << endl; 

  save::  importDat<double>(sMsDiff, folderName + "/" + fileName);
    
  if (pltVerbose) {
    std::vector<PLOToptions> pOpts(2);
    std::vector<std::string> pVals(2);
    pOpts[0] = minimum;  pVals[0] = "-2e-2";
    pOpts[1] = maximum;  pVals[1] = "2e-2";
    plt.print1d(sMsDiff, "./plots/importedTDdiff", pOpts, pVals);
  }

  NQbins = shape[0];
  float filterVar = std::pow(NQbins/filtSTD, 2);


  /////  Check Varaiables  /////

  if ((NQbins == -1)
        || (NPadBins  == -1)
        || (QperPix == -1)
        || (filtSTD == -1)) {
 
    cerr << "Check Parameters NQbins ("
        << NQbins << ")  NPadBins ("
        << NPadBins <<")  QperPix ("
        << QperPix <<")  filtSTD ("
        << filtSTD << ")\n";
    exit(0);
  } 


  ///////////////////////////////////////
  /////  Pair correlation function  /////
  ///////////////////////////////////////
  
  if (verbose) {
    std::cout << "Begin calculating pair correlation functions\n";
  }

  bool mirrorDiffPattern = false;
  int inpFFTsize = NQbins + NPadBins;
  int outFFTsize = inpFFTsize/2 + 1;
  int outSize = rMaxRat*outFFTsize;
  cout<<"outsize: "<<outSize<<endl;

  double* qSpace = (double*) fftw_malloc(sizeof(double)*inpFFTsize);
  fftw_complex* rSpace =
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*outFFTsize);
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(inpFFTsize, qSpace, rSpace, FFTW_MEASURE);

  int FFTcent = inpFFTsize/2;
  std::vector<double> fftInp(inpFFTsize);
  std::vector<double> smoothed(NQbins, 0);
  std::vector<double> filter(NQbins, 0);
  std::vector<double> pairCorrEven(outSize);
  std::vector<double> pairCorrOdd(outSize);
  std::vector< std::vector<double> > results;


  double filterScale = 1;

  int iq = 0;
  /////  Filling FFT Input  /////
  for (iq=0; iq<NQbins; iq++) {
    //  Apply Filter
    if (doFiltering) {
      filterScale = exp(-1*std::pow(iq, 2)
                        /(2*filterVar));
    }
    smoothed[iq] = sMsDiff[iq]*filterScale;
    filter[iq] = filterScale;
  }

  if (true || pltVerbose) {
    plt.print1d(sMsDiff, "./plots/inpDiff_" + molName);
    plt.print1d(smoothed, "./plots/smoothed_" + molName);
    plt.print1d(filter, "./plots/filter");
  }


  /////  Apply FFT  /////
  for (iq=0; iq<NQbins; iq++) {
    fftInp[iq] = smoothed[iq];
  }
  results = tools::fft1dRtoC(fftInp, fftB, qSpace, rSpace, 
                false, false, NULL, "fftInp_" + molName);

  if (pltVerbose)
    plt.print1d(fftInp, "./plots/fftFuncInp_" + molName);

  for (int ir=0; ir<outSize; ir++) {
    pairCorrEven[ir] = std::pow(results[0][ir], 1);
    pairCorrOdd[ir] = std::pow(results[1][ir], 1);
  }
  

  if (pltVerbose) {
    vals[1] = "0," + to_string(rMaxRat*outFFTsize*(2*PI/(QperPix*inpFFTsize)));
    if (doFiltering) {
      plt.print1d(pairCorrEven, "./plots/pairCorrEven", opts, vals);
      plt.print1d(pairCorrOdd, "./plots/pairCorrOdd", opts, vals);
    }
    else {
      plt.print1d(pairCorrEven, "./plots/pairCorrEven_sanityCheck", opts, vals);
      plt.print1d(pairCorrOdd, "./plots/pairCorrOdd_sanityCheck", opts, vals);
    }
  }


  ////////////////////
  /////  Saving  /////
  ////////////////////

  if (verbose)
    std::cout << "Saving.\n";

  save::saveDat<double>(pairCorrEven, 
      "./results/pairCorrEven_"   
      + molName + "_" 
      + bondType + ".dat");
  save::saveDat<double>(pairCorrOdd, 
      "./results/pairCorrOdd_"   
      + molName + "_"
      + bondType + ".dat");


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (verbose) 
    std::cout << "Cleaning up" << endl;

  fftw_free(rSpace);
  fftw_free(qSpace);
  fftw_destroy_plan(fftB);
  return 1;
}
