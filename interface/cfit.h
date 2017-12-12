// CFIT class definition and includes

#ifndef CFIT_H
#define CFIT_H

#include "TROOT.h"
#include "THStack.h"
#include "Riostream.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompLU.h"
#include "TDecompQRH.h"
#include "TDecompChol.h"
#include "TDecompBK.h"

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <memory>

const double PI = 3.141592653589793238462;

const int NBINMAX = 1000; // the max number of bins
const int NTEMPLATEMAX = 100; // the max number of templates
const int SYSMAX = 100; // the max number of systematics
const int HMAX = 10000; // the total max number of histograms

const bool doRescaleSys = 1;
const bool doNormalise = 1;

extern double funcDEFM(double vdata1,
		       double vdataErr1,
		       std::vector<double> vmc1,
		       double vdata2,
		       double vdataErr2,
		       std::vector<double> vmc2,
		       double covEL,
		       double norm1,
		       double norm2,
		       double meanBias1,
		       double meanBias2,
		       double *par);

extern void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag);

typedef enum {
   OPT_NONE,
   OPT_MORPH,
   OPT_SGN,
   OPT_SIGMA,
   OPT_MORPH_SGN,
   OPT_MORPH_SIGMA,
   OPT_SGN_SIGMA,
   OPT_MORPH_SGN_SIGMA,
   OPT_NOCORR
}
OPTMODE;

typedef enum {
   OPTMORPH_CUTOFF,
   OPTMORPH_GEOMETRIC,
   OPTMORPH_ARITHMETIC,
   OPTMORPH_FUNC
}
OPTMORPHMODE;

typedef enum {
   COV_MAX,
   COV_BARLOW
}
COVMODE;

namespace CFIT
{
   
   class cfit
     {
	
      public:
	cfit(std::string name = "");
	virtual ~cfit();

	// user methods
	
	// run the fit procedure
	// option="": perform the fit with pre-tag templates
	// option="tag": perform the fit with tag templates 
	void Run(std::string option = "");
	
	// switch on verbose
	// 0: no verbose
	// 1: verbose on
	void SetVerbose(int v);

	// set optimization option for systematics correlation treatment
	// OPT_NONE: no optimization (DEFAULT)
	// OPT_MORPH: perform the morphing of systematic correlations according to OPTMORPHMODE
	// OPT_SGN: ignore correlations in bins where systematics flips sign or is of the same sign
 	// OPT_SIGMA: reset relative bin uncertainties to sig = sqrt[ sum_i(sig_i*sig_i*b_i)/sum_i(b_i) ]
	//            if sig_i/sig < 0.1
	// combination of various options is also possible:
	// OPT_MORPH_SGN, OPT_MORPH_SIGMA, OPT_SGN_SIGMA, OPT_MORPH_SGN_SIGMA
	void SetOptimization(OPTMODE optMode=OPT_NONE);
	
	int GetOptimization() {return optMode;};
	
	// set morphing correlation factor for systematic variations
	// OPTMORPH_CUTOFF: f=0 if b_i/b_max < frac
	// OPTMORPH_GEOMETRIC: f=sqrt(b_i*b_j)/b_max (DEFAULT)
	// OPTMORPH_ARITHMETIC: f=(b_i+b_j)/(2*b_max)
	// OPTMORPH_FUNC: f=b_i*b_j/(b_max*b_max)
	void SetMorphing(OPTMORPHMODE optMorphMode=OPTMORPH_GEOMETRIC,float frac = 1.);
	
	int GetMorphing() {return optMorphMode;};

	// set variance computation mode
	// COV_MAX: symmetrise systematic uncertainty by choosing the maximum deviation (DEFAULT)
	// COV_BARLOW: asymmetric systematic uncertainty approximation by two half-gaussians (Model #1)
	void SetCovarianceMode(COVMODE covMode=COV_MAX);
	
	int GetCovarianceMode() {return *covMode;};
	
	// produce plots and store results in file
	void ProducePlots(bool produce = 1);
	
	// return chi2 of the fit
	double GetChisq();

	// return nDOF
	int GetNDOF();
	
	// get the number of templates
	int GetNPar();
	
	// get the i fitted parameter
	double GetPar(int i);
	
	// get the error of the i fitted parameter
	double GetParErr(int i);
	
	// set the name of the input file with histograms
	void SetInputFile(std::string fin);

	// set the name of the nominal data histogram
	void SetData(std::string name);
	
	// set the name of the post-tag nominal data histogram
	void SetDataTag(std::string name);

	// set the name of the un-tag nominal data histogram
	void SetDataUntag(std::string name);
	
	// add the title and histogram name of the nominal template
	void AddTemplate(std::string name,std::string nominalName,int colour);
	
	// add the title and histogram name of the tag template
	void AddTemplateTag(std::string name,std::string tagName,int colour);

	// add the title and histogram name of the un-tag template
	void AddTemplateUntag(std::string name,std::string untagName,int colour);
	
	// vary templates simultaneously in the fit
	void GlueTemplates(std::vector<std::string> nameV,std::string nameMerged = "merged",int colour = 46);
	
	// vary templates simultaneously in the fit
	void GlueTemplatesTag(std::vector<std::string> nameV,std::string nameMerged = "merged",int colour = 46);
	
	// set the write/read option for correlation matrix
	// option="WRITE": write the result matrix to file
	// option="READ": read the matrix from file
	void SetMatrixOption(std::string option);

	// set the output name for root file with matrix
	void SetMatrixName(std::string name = "matrix");

	// set histogram legend header
	void SetLegendHeader(std::string name = "");
	
	// add systematics template
	// name = name of the systematics, could be any name
	// up = histogram prefix name for the up variation with the full histogram name to be nominalName+up
	// down = histogram prefix name for the down variation with the full histogram name to be nominalName+down
	void AddSys(std::string name,std::string upName,std::string downName);
	
	// replace nominal template by a systematic variation
	// name = either upName or downName
	void SetSysVariation(std::string name);
	
	std::string GetSysVariation() {return nameSYSVAR;};
	
	// vary nominal template within the statistical uncertainty
	// rnd = seed (>= 666)
	void SetStatVariation(int rnd);
	
	int GetStatVariation() {return sysIdx;};

	// get total sum of events in the input data template
	float GetNData();
	
	// get total sum of weighted events for the given template
	float GetNTemplate(std::string templateName);

	// get total error for the sum of weighted events for the given template
	float GetErrTemplate(std::string templateName);

	// get total error for the sum of weighted events for the given template
	void DrawPrePostFit() {drawPPF = 1;};
	
      protected:

	// functions
	
	void reset();
	
	void resetNames();

	void setMatrixIO(std::string option);
	
	void processInput(std::string option);
	
	void processGlueTemplates(std::string option);

	void sortTemplateNames();
	
	int getTemplateId(std::string templateName);
	
	void SetPlotStyle();
	
	TStyle* PlotStyle();
	
	void totSys(TH1D *h_nom,TH1D *h_down,TH1D *h_up);
	
	void combSysLinear(std::vector<TH1D*> h_nom,
			   std::vector<TH1D*> h_sys_up,
			   std::vector<TH1D*> h_sys_down,
			   TH1D *h_comb,
			   TH1D *h_sys_down_comb,
			   TH1D *h_sys_up_comb);
	
	void combSys(TH1D *h_nom,TH1D** h_sys_down,TH1D** h_sys_up,TH1D *h_sys_down_comb,TH1D *h_sys_up_comb);

	void applySys(std::vector<double> sf,
		      std::vector<double> sferr,
		      bool write,
		      std::string option);

	void adjustSys(TH1D *hnom, TH1D *hsys_down, TH1D *hsys_up, int isys);
	
	void doFit(double *chis,std::string option);
	
	void doFracSys(TH1D *hnom,TH1D *hsysDown,TH1D *hsysUp,int isys,std::string option,std::string postfix);
	
	double bbCorr(TH1D *hsysDown,TH1D *hsysUp);
	
	void applySF(std::string option);
	
	void addbin(TH1D *h);

	void removeMCNeg(TH1D *h1);
	
	void adjust(TH1D *h);
	
	void normalise(TH1D *h1);
	
	void rescale(TH1D *hnom,TH1D *hsysUp,TH1D *hsysDown);

	double* mergeBins(TH1D *h,int ibmin,int ibmax,double xnew[]);
	
	TGraphAsymmErrors* makeErrorBand(TH1D* tot, TH1D* plus, TH1D* minus);
	
	void removeMCSys(TH1D *hnom,TH1D *hsysUp,TH1D *hsysDown);
	
	void drawPrePostFit(std::string option,bool postfit);
	
	// variables
	
	double chis[10];
	
	int sysIdx;
	int optMode;
	int optMorphMode;
	float optMorphFrac;
	std::vector<double> sf;
	std::vector<double> sferr;
	double covFactor[SYSMAX];
	double binFactor[SYSMAX][NBINMAX];
	double binFactor2D[SYSMAX][NBINMAX][NBINMAX];
	TFile *fcov;
	std::string covOption;
	std::string covName;
	std::string legendName;
	std::string runName;
	bool producePlots;
	float _NDATA;
	float _NMC[NTEMPLATEMAX];
	float _NMCERR[NTEMPLATEMAX];
	int ibinMax;
	int ibinMin;
	int nTAG;
	int nUNTAG;
	int nSYS;
	bool rcov;
	bool drawPPF;
	std::string nameDATA_CURRENT;
	std::string nameDATA;
	std::string nameDATATAG;
	std::string nameDATAUNTAG;
	std::vector<std::string> nameT_CURRENT;
	std::vector<std::string> nameT;
	std::vector<std::string> nameTAG;
	std::vector<std::string> nameUNTAG;
	std::vector<std::string> titleT;
	std::vector<std::string> titleTAG;
	std::vector<std::string> titleUNTAG;
	std::vector<std::vector<std::string> > titleGLUE;
	std::vector<std::vector<std::string> > titleGLUETAG;
	std::vector<int> colourT;
	std::vector<int> colourTAG;
	std::vector<int> colourUNTAG;
	std::vector<std::string> labelGLUE;
	std::vector<std::string> labelGLUETAG;
	std::vector<int> colourGLUE;
	std::vector<int> colourGLUETAG;
	std::string nameSYSVAR;
	std::vector<std::string> nameSYS;
	std::vector<std::string> nameSYSUP;
	std::vector<std::string> nameSYSDOWN;
	int ISMERGEDFIT[NTEMPLATEMAX];
	double cov[NBINMAX][NBINMAX];
	double norm1D[NBINMAX];
	
	double PAR[NTEMPLATEMAX];
	double ERR[NTEMPLATEMAX];
	
	TFile *dfile;
	TH1D *hist[HMAX];
	
	TH1D *h_comb;
	TH1D *h_sys_down_comb[SYSMAX];
	TH1D *h_sys_up_comb[SYSMAX];
	TH1D *h_combNATURAL;
	TH1D *h_sys_down_combNATURAL[SYSMAX];
	TH1D *h_sys_up_combNATURAL[SYSMAX];
	TH1D *h_sys_down_combSum;
	TH1D *h_sys_up_combSum;

      public:
	
	static std::unique_ptr<int> covMode; 
	
	static std::unique_ptr<std::vector<int> > vecM;
	static std::unique_ptr<int> nT;
	static std::unique_ptr<int> nBINS;
	static std::unique_ptr<double> CHISQ;
	static std::unique_ptr<int> NDOF;
	static std::unique_ptr<bool> VALID;

	static std::vector<std::unique_ptr<TH1D> > histNOM;
	static std::vector<std::unique_ptr<TH1D> > histNOMTRUE;
	static std::vector<std::unique_ptr<TH1D> > histTAGNOM;
	static std::vector<std::unique_ptr<TH1D> > histUNTAGNOM;
	static std::vector<std::unique_ptr<TH1D> > histDOWN;
	static std::vector<std::unique_ptr<TH1D> > histUP;
	
	static std::unique_ptr<TH1D> h_data;
	static std::unique_ptr<TH1D> h_data_tag;
	static std::unique_ptr<TH1D> h_data_untag;

	static std::unique_ptr<TMatrixD> covMIp;
	static std::unique_ptr<TMatrixD> covMp;
	static std::unique_ptr<TVectorD> norm1Dp;

	static std::unique_ptr<bool> verb;
	
	static std::unique_ptr<std::vector<int> > PARIDXFIT;
	
	static std::unique_ptr<int> nTFIT;
	static std::unique_ptr<int> nTTAGFIT;
	
      ClassDef(CFIT::cfit,1)
     };
}

#endif
