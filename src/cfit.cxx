// CFIT class implementation
#include "RecoBTag/CFIT/interface/cfit.h"

ClassImp(CFIT::cfit)

double funcDEFM(double vdata1,
		double vdataErr1,
		std::vector<double> vmc1,
		double vdata2,
		double vdataErr2,
		std::vector<double> vmc2,
		double covEL,
		double norm1,
		double norm2,
		double *par);

void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag);

std::unique_ptr<std::vector<int> > CFIT::cfit::vecM;

std::unique_ptr<std::vector<int> > CFIT::cfit::PARIDXFIT;

std::unique_ptr<int> CFIT::cfit::covMode;
std::unique_ptr<int> CFIT::cfit::nT;
std::unique_ptr<int> CFIT::cfit::nTFIT;
std::unique_ptr<int> CFIT::cfit::nBINS;
std::unique_ptr<double> CFIT::cfit::CHISQ;
std::unique_ptr<int> CFIT::cfit::NDOF;
std::unique_ptr<bool> CFIT::cfit::VALID;

std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histNOM;
std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histNOMTRUE;
std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histTAGNOM;
std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histUNTAGNOM;
std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histDOWN;
std::vector<std::unique_ptr<TH1D> > CFIT::cfit::histUP;
std::unique_ptr<TH1D> CFIT::cfit::h_data;
std::unique_ptr<TH1D> CFIT::cfit::h_data_tag;
std::unique_ptr<TH1D> CFIT::cfit::h_data_untag;

std::unique_ptr<TMatrixD> CFIT::cfit::covMIp;
std::unique_ptr<TMatrixD> CFIT::cfit::covMp;
std::unique_ptr<TVectorD> CFIT::cfit::norm1Dp;

std::unique_ptr<bool> CFIT::cfit::verb;

CFIT::cfit::cfit(std::string name)
{
   runName = name;
   
   covMIp.reset(); covMIp = std::unique_ptr<TMatrixD>(new TMatrixD(NBINMAX,NBINMAX));
   covMp.reset(); covMp = std::unique_ptr<TMatrixD>(new TMatrixD(NBINMAX,NBINMAX));
   norm1Dp.reset(); norm1Dp = std::unique_ptr<TVectorD>(new TVectorD(NBINMAX));
   
   covMode.reset(); covMode = std::unique_ptr<int>(new int(0));
   nT.reset(); nT = std::unique_ptr<int>(new int(0));
   nTFIT.reset(); nTFIT = std::unique_ptr<int>(new int(0));
   nBINS.reset(); nBINS = std::unique_ptr<int>(new int(0));
   CHISQ.reset(); CHISQ = std::unique_ptr<double>(new double(0));
   NDOF.reset(); NDOF = std::unique_ptr<int>(new int(0));
   VALID.reset(); VALID = std::unique_ptr<bool>(new bool(0));
   
   verb.reset(); verb = std::unique_ptr<bool>(new bool(0));
   
   vecM.reset(); vecM = std::unique_ptr<std::vector<int> >(new std::vector<int>);

   PARIDXFIT.reset(); PARIDXFIT = std::unique_ptr<std::vector<int> >(new std::vector<int>);
   
   optMode = OPT_NONE;
   optMorphMode = OPTMORPH_GEOMETRIC;
   optMorphFrac = 1.;
   *covMode = COV_MAX;

   covOption = "";
   covName = "matrix";
   legendName = "";
   
   producePlots = 0;

   *verb = 0;
   
   drawPPF = 0;
}

void CFIT::cfit::SetInputFile(std::string fin)
{
   std::ifstream infile(fin.c_str());
   
   if( infile.good() )
     { 	
	if( dfile && dfile->IsOpen() ) dfile->Close();
	
	dfile = TFile::Open(fin.c_str(),"READ");
     }
   else
     {
	std::cout << "Input file does not exist: " << fin << std::endl;
	reset();
	exit(1);
     }      
}

void CFIT::cfit::ProducePlots(bool produce)
{
   producePlots = produce;
}

void CFIT::cfit::SetOptimization(OPTMODE mode)
{
   bool optFound = 0;
   
   for(int i=OPT_NONE;i!=OPT_NOCORR+1;i++)
     {
	OPTMODE cur = static_cast<OPTMODE>(i);
	if( mode == cur )
	  {
	     optFound = 1;
	     break;
	  }	
     }   
   if( !optFound )
     {
	std::cout << "Optimization mode " << mode << " is not defined" << std::endl;
	reset();
	exit(1);
     }
   else
     optMode = mode;
}

void CFIT::cfit::SetMorphing(OPTMORPHMODE mode,float frac)
{
   bool optFound = 0;
   
   for(int i=OPTMORPH_CUTOFF;i!=OPTMORPH_FUNC+1;i++)
     {
	OPTMORPHMODE cur = static_cast<OPTMORPHMODE>(i);
	if( mode == cur )
	  {
	     optFound = 1;
	     break;
	  }	
     }   
   if( !optFound )
     {
	std::cout << "Morphing optimization mode " << mode << " is not defined" << std::endl;
	reset();
	exit(1);
     }
   else
     {	
	optMorphMode = mode;
	optMorphFrac = frac;
     }   
}

void CFIT::cfit::SetCovarianceMode(COVMODE mode)
{
   bool optFound = 0;
   
   for(int i=COV_MAX;i!=COV_BARLOW+1;i++)
     {
	COVMODE cur = static_cast<COVMODE>(i);
	if( mode == cur )
	  {
	     optFound = 1;
	     break;
	  }	
     }   
   if( !optFound )
     {
	std::cout << "Covariance computation mode " << mode << " is not defined" << std::endl;
	reset();
	exit(1);
     }
   else
     {	
	*covMode = mode;
     }   
}

void CFIT::cfit::SetData(std::string name)
{
   nameDATA = name;
   nameDATA_CURRENT = name;
}

void CFIT::cfit::SetDataTag(std::string name)
{
   nameDATATAG = name;
}

void CFIT::cfit::SetDataUntag(std::string name)
{
   nameDATAUNTAG = name;
}

void CFIT::cfit::AddTemplate(std::string name,std::string nameNominal,int colour)
{
   titleT.push_back(name);
   nameT.push_back(nameNominal);
   nameT_CURRENT.push_back(nameNominal);
   colourT.push_back(colour);
}

void CFIT::cfit::AddTemplateTag(std::string name,std::string nameTag,int colour)
{
   titleTAG.push_back(name);
   nameTAG.push_back(nameTag);
   colourTAG.push_back(colour);
}

void CFIT::cfit::AddTemplateUntag(std::string name,std::string nameUntag,int colour)
{
   titleUNTAG.push_back(name);
   nameUNTAG.push_back(nameUntag);
   colourUNTAG.push_back(colour);
}

void CFIT::cfit::GlueTemplates(std::vector<std::string> nameV,std::string nameMerged,int colour)
{
   titleGLUE.push_back(nameV);
   labelGLUE.push_back(nameMerged);
   colourGLUE.push_back(colour);
}

void CFIT::cfit::GlueTemplatesTag(std::vector<std::string> nameV,std::string nameMerged,int colour)
{
   titleGLUETAG.push_back(nameV);
   labelGLUETAG.push_back(nameMerged);
   colourGLUETAG.push_back(colour);
}

void CFIT::cfit::processGlueTemplates(std::string option)
{
   int nReq = titleGLUE.size();
   if( strcmp(option.c_str(),"tag") == 0 ) nReq = titleGLUETAG.size();
   
   *nTFIT = *nT+nReq;
   
   int nGlueTotal = 0;
   
   for(int ir=0;ir<nReq;ir++)
     {
	std::vector<int> glueIndex;
	
	int nGlue = 0;
	if( strcmp(option.c_str(),"tag") == 0 ) nGlue = titleGLUETAG[ir].size();
	else nGlue = titleGLUE[ir].size();
	nGlueTotal += nGlue;
	
	for(int it=0;it<*nT;it++)
	  {
	     std::string titName = "";
	     if( strcmp(option.c_str(),"tag") == 0 ) titName = titleTAG[it];
	     else titName = titleT[it];
	     
	     for(int ig=0;ig<nGlue;ig++)
	       {
		  std::string gluName = "";
		  if( strcmp(option.c_str(),"tag") == 0 ) gluName = titleGLUETAG[ir][ig];
		  else gluName = titleGLUE[ir][ig];
		  
		  if( strcmp(gluName.c_str(),titName.c_str()) == 0 )
		    {
		       ISMERGEDFIT[it] = ir+1;
		       break;
		    }		  
	       }
	  }	
     }   
   
   int OUT[*nT];
   for(int it=0;it<*nT;it++)
     {
	OUT[it] = 0;	
     }   
   int k = 0;
   
   for(int it=1;it<*nT;it++)
     {
	if( ISMERGEDFIT[it] == 0 )
	  {
	     k++;
	     OUT[it] = k;
	  }	
	else
	  {
	     for(int j=0;j<it;j++)
	       {
		  if( ISMERGEDFIT[it] == ISMERGEDFIT[j] ) OUT[it] = OUT[j];
	       }	     
	     
	     if( OUT[it] == 0 )
	       {
		  k++;
		  OUT[it] = k;
	       }	     
	  }	
     }   
   
   if( nReq != 0 ) *nTFIT -= nGlueTotal;
   else *nTFIT = *nT;
   
   for(int it=0;it<*nT;it++)
     {
	PARIDXFIT->at(it) = OUT[it];
//	std::cout << PARIDXFIT->at(it) << " " << ISMERGEDFIT[it] << std::endl;
     }   
}

void CFIT::cfit::sortTemplateNames()
{
   // tag
   std::vector<int> nameId_tag;
   
   for(int it=0;it<*nT;it++)
     {
	std::string tit = titleT[it];
	for(int itt=0;itt<nTAG;itt++)
	  {
	     if( tit == titleTAG[itt] )
	       {
		  nameId_tag.push_back(itt);
		  break;
	       }	     
	  }	
     }   
   
   if( nameId_tag.size() != (unsigned)*nT && nameId_tag.size() > 0 )
     {
	std::cout << "Templates in pre-tag and tag collections do not match" << std::endl;
	reset();
	exit(1);
     }

   if( nameId_tag.size() > 0 )
     {	
	std::vector<std::string> titleTAG_sorted = titleTAG;
	std::vector<std::string> nameTAG_sorted = nameTAG;
	std::vector<int> colourTAG_sorted = colourTAG;
	nameTAG.clear();
	titleTAG.clear();
	colourTAG.clear();
	
	for(int it=0;it<*nT;it++)
	  {
	     int idx = nameId_tag[it];
	     nameTAG.push_back(nameTAG_sorted[idx]);
	     titleTAG.push_back(titleTAG_sorted[idx]);
	     colourTAG.push_back(colourTAG_sorted[idx]);
	  }   
     }   

   // untag
   std::vector<int> nameId_untag;
   
   for(int it=0;it<*nT;it++)
     {
	std::string tit = titleT[it];
	for(int itt=0;itt<nUNTAG;itt++)
	  {
	     if( tit == titleUNTAG[itt] )
	       {
		  nameId_untag.push_back(itt);
		  break;
	       }	     
	  }	
     }   
   
   if( nameId_untag.size() != (unsigned)*nT && nameId_untag.size() > 0 )
     {
	std::cout << "Templates in pre-tag and un-tag collections do not match" << std::endl;
	reset();
	exit(1);
     }

   if( nameId_untag.size() > 0 )
     {	   
	std::vector<std::string> titleUNTAG_sorted = titleUNTAG;
	std::vector<std::string> nameUNTAG_sorted = nameUNTAG;
	std::vector<int> colourUNTAG_sorted = colourUNTAG;
	nameUNTAG.clear();
	titleUNTAG.clear();
	colourUNTAG.clear();
	
	for(int it=0;it<*nT;it++)
	  {
	     int idx = nameId_untag[it];
	     nameUNTAG.push_back(nameUNTAG_sorted[idx]);
	     titleUNTAG.push_back(titleUNTAG_sorted[idx]);
	     colourUNTAG.push_back(colourUNTAG_sorted[idx]);
	  }   
     }   
}

void CFIT::cfit::setMatrixIO(std::string option)
{
   if( covOption == "READ" )
     {	
	if( option == "" )
	  {		  
	     std::string name = covName+".root";
	     if( fcov ) fcov->Close();
	     delete gROOT->GetListOfFiles()->FindObject(name.c_str());
	     fcov = TFile::Open(name.c_str(),"READ");
	     rcov = 1;
	     if( *verb ) std::cout << "Read "<< name << std::endl;
	  }
	else if( option == "tag" )
	  {		  
	     std::string name = covName+"_tag.root";
	     if( fcov ) fcov->Close();
	     delete gROOT->GetListOfFiles()->FindObject(name.c_str());
	     fcov = TFile::Open(name.c_str(),"READ");
	     rcov = 1;
	     if( *verb ) std::cout << "Read "<< name << std::endl;
	  }
	else
	  {
	     std::cout << "Option " << option << " is not known" << std::endl;
	     reset();
	     exit(1);
	  }      
     }   
   else
     {
	if( *verb ) std::cout << "Create new correlation matrix" << std::endl;
	if( fcov ) fcov->Close();
	fcov = NULL;
     }
}

void CFIT::cfit::SetMatrixOption(std::string option)
{
   covOption = option;
   if( covOption != "READ" && covOption != "WRITE" )
     {
	std::cout << "Matrix options READ and WRITE are only allowed" << std::endl;
	reset();
	exit(1);
     }
}

void CFIT::cfit::SetMatrixName(std::string name)
{
   covName = name;
}

void CFIT::cfit::SetLegendHeader(std::string name)
{
   legendName = name;
}

void CFIT::cfit::AddSys(std::string name,std::string up,std::string down)
{   
   nameSYS.push_back(name);
   nameSYSDOWN.push_back(down);
   nameSYSUP.push_back(up);
}

void CFIT::cfit::SetSysVariation(std::string name)
{   
   nameSYSVAR = name;
}

void CFIT::cfit::SetStatVariation(int rnd)
{   
   sysIdx = rnd;
}

void CFIT::cfit::reset()
{
   histDOWN.clear();
   histUP.clear();
   
   if( h_data.get() ) h_data.reset();
   if( h_data_tag.get() ) h_data_tag.reset();
   if( h_data_untag.get() ) h_data_untag.reset();
   
   histNOM.clear();
   histNOMTRUE.clear();
   histTAGNOM.clear();
   histUNTAGNOM.clear();
   
   if( h_comb ) delete h_comb;
   if( h_combNATURAL ) delete h_combNATURAL;
   if( h_sys_down_combSum ) delete h_sys_down_combSum;
   if( h_sys_up_combSum ) delete h_sys_up_combSum;
   
   for(int is=0;is<nSYS;is++)
     {	
	if( h_sys_down_comb[is] ) delete h_sys_down_comb[is];
	if( h_sys_up_comb[is] ) delete h_sys_up_comb[is];
	if( h_sys_down_combNATURAL[is] ) delete h_sys_down_combNATURAL[is];
	if( h_sys_up_combNATURAL[is] ) delete h_sys_up_combNATURAL[is];
     }   
   
   ibinMax = NBINMAX;
   ibinMin = 1;

   if( vecM.get() ) {vecM.reset(); vecM = std::unique_ptr<std::vector<int> >(new std::vector<int>);}

   *VALID = 0;

   covMIp.reset(); covMIp = std::unique_ptr<TMatrixD>(new TMatrixD(NBINMAX,NBINMAX));
   covMp.reset(); covMp = std::unique_ptr<TMatrixD>(new TMatrixD(NBINMAX,NBINMAX));
   norm1Dp.reset(); norm1Dp = std::unique_ptr<TVectorD>(new TVectorD(NBINMAX));
   
   PARIDXFIT.reset(); PARIDXFIT = std::unique_ptr<std::vector<int> >(new std::vector<int>);
   
   for(int i=0;i<NTEMPLATEMAX;i++)
     {
	PAR[i] = 1.;
	ERR[i] = 0.;
	PARIDXFIT->push_back(i);
	ISMERGEDFIT[i] = 0;
     }   
}

void CFIT::cfit::resetNames()
{
   nameDATA = "";
   nameDATA_CURRENT = "";
   nameDATATAG = "";
   nameDATAUNTAG = "";
   
   nameT.clear();
   nameT_CURRENT.clear();
   nameTAG.clear();
   nameUNTAG.clear();
   
   titleT.clear();
   titleTAG.clear();
   titleUNTAG.clear();

   colourT.clear();
   colourTAG.clear();
   colourUNTAG.clear();
   
   nameSYS.clear();
   nameSYSDOWN.clear();
   nameSYSUP.clear();
}

void CFIT::cfit::processInput(std::string option)
{
   reset();
   if( *verb ) std::cout << "Reset done" << std::endl;

   if( option != "tag" && producePlots )
     system("if [[ ! -d pics ]]; then mkdir pics; else rm -r pics; mkdir pics; fi");
   
   gErrorIgnoreLevel = 2000;

   SetPlotStyle();

   rcov = 0;
   
   setMatrixIO(option);
   if( *verb ) std::cout << "Update matrix option done" << std::endl;

   nSYS = nameSYS.size();

   for(int is=0;is<nSYS;is++)
     {
	covFactor[is] = 1.;
	for(int ib=0;ib<NBINMAX;ib++)
	  {
	     binFactor[is][ib] = 1.;
	     for(int ibb=0;ibb<NBINMAX;ibb++)
	       {
		  binFactor2D[is][ib][ibb] = 1.;
	       }	     
	  }	
     }   
   
   *nT = nameT.size();
   if( *nT == 0 )
     {
	std::cout << "Please specify at least one template" << std::endl;
	reset();
	exit(1);
     }   

   nTAG = nameTAG.size();
   if( nTAG != *nT && nTAG > 0 )
     {
	std::cout << "Please specify the same number of tagged templates" << std::endl;
	reset();
	exit(1);
     }   

   nUNTAG = nameUNTAG.size();
   if( nUNTAG != *nT && nTAG > 0 )
     {
	std::cout << "Please specify the same number of un-tagged templates" << std::endl;
	reset();
	exit(1);
     }   

   if( option == "tag" )
     {
	nameDATA_CURRENT = nameDATATAG;
	
	nameT_CURRENT.clear();
	for(int it=0;it<*nT;it++)
	  {	     
	     nameT_CURRENT.push_back(nameTAG[it]);
	  }	
     }   
   else
     {
	nameDATA_CURRENT = nameDATA;
	
	nameT_CURRENT.clear();
	for(int it=0;it<*nT;it++)
	  {
	     nameT_CURRENT.push_back(nameT[it]);
	  }	
     }
   
   if( nameDATA_CURRENT != "" )
     {
	TH1D *hdNOM = (TH1D*)dfile->Get(nameDATA_CURRENT.c_str());
	if( ! hdNOM )
	  {	     
	     std::cout << "Data histogram does not exist: " << nameDATA_CURRENT << std::endl;
	     reset();
	     exit(1);
	  }	
	std::string hdcopyNOM = nameDATA_CURRENT+"_NOMcopy";
	h_data = std::unique_ptr<TH1D>( (static_cast<TH1D*>(hdNOM->Clone(hdcopyNOM.c_str()))) );
	delete hdNOM;

	if( h_data->Integral() == 0. )
	  {
	     std::cout << "Data histogram is empty" << std::endl;
	     reset();
	     exit(1);
	  }

	// adjust empty bins in data
	adjust(h_data.get());
	
	_NDATA = h_data->Integral(0,h_data->GetXaxis()->GetNbins()+1);
     }
   else
     {
	std::cout << "Please specify the nominal data histogram" << std::endl;
	reset();
	exit(1);
     }

   if( nameDATATAG != "" ) 
     {
	TH1D *hdTAG = (TH1D*)dfile->Get(nameDATATAG.c_str());

	if( ! hdTAG )
	  {
	     std::cout << "Data tag histogram does not exist: " << nameDATATAG << std::endl;
	     reset();
	     exit(1);
	  }	

	std::string hdcopyTAG = nameDATATAG+"_TAGcopy";
	h_data_tag = std::unique_ptr<TH1D>( (static_cast<TH1D*>(hdTAG->Clone(hdcopyTAG.c_str()))) );
	delete hdTAG;
	
	if( h_data_tag->Integral() == 0. )
	  {
	     std::cout << "Data tag histogram is empty" << std::endl;
	     reset();
	     exit(1);
	  }	
	TH1D *hdUNTAG = (TH1D*)dfile->Get(nameDATAUNTAG.c_str());

	if( ! hdUNTAG )
	  {
	     if( nameDATAUNTAG == "" ) std::cout << "Please also specify the un-tag data histogram" << std::endl;
	     else std::cout << "Histogram does not exist: " << nameDATAUNTAG << std::endl;
	     reset();
	     exit(1);
	  }	

	std::string hdcopyUNTAG = nameDATAUNTAG+"_UNTAGcopy";
	h_data_untag = std::unique_ptr<TH1D>( (static_cast<TH1D*>(hdUNTAG->Clone(hdcopyUNTAG.c_str()))) );
	delete hdUNTAG;
	
	if( h_data_untag->Integral() == 0. )
	  {
	     std::cout << "Data untag histogram is empty" << std::endl;
	     reset();
	     exit(1);
	  }	
	int nbNom = h_data->GetXaxis()->GetNbins();
	int nbTag = h_data_tag->GetXaxis()->GetNbins();
	int nbUntag = h_data_untag->GetXaxis()->GetNbins();
	float xMinNom = h_data->GetXaxis()->GetBinLowEdge(1);
	float xMaxNom = h_data->GetXaxis()->GetBinUpEdge(nbNom);
	float xMinTag = h_data_tag->GetXaxis()->GetBinLowEdge(1);
	float xMaxTag = h_data_tag->GetXaxis()->GetBinUpEdge(nbTag);
	float xMinUntag = h_data_untag->GetXaxis()->GetBinLowEdge(1);
	float xMaxUntag = h_data_untag->GetXaxis()->GetBinUpEdge(nbUntag);
	if( nbNom != nbTag )
	  {
	     std::cout << "Pre-tag and tag histograms should have the same number of bins" << std::endl;
	     reset();
	     exit(1);
	  }	
	if( nbNom != nbUntag )
	  {
	     std::cout << "Pre-tag and tag histograms should have the same number of bins" << std::endl;
	     reset();
	     exit(1);
	  }	
	if( xMinNom != xMinTag || xMaxNom != xMaxTag )
	  {
	     std::cout << "Pre-tag and tag histograms should have the same x-axis ranges" << std::endl;
	     reset();
	     exit(1);
	  }
	if( xMinNom != xMinUntag || xMaxNom != xMaxUntag )
	  {
	     std::cout << "Pre-tag and un-tag histograms should have the same x-axis ranges" << std::endl;
	     reset();
	     exit(1);
	  }
     }

   if( nameSYSVAR != "" && sysIdx < 666. )
     {	
	bool foundSYSVAR = 0;
	for(int is=0;is<nSYS;is++)
	  {
	     bool isDown = (strcmp(nameSYSVAR.c_str(),nameSYSDOWN[is].c_str()) == 0);
	     bool isUp = (strcmp(nameSYSVAR.c_str(),nameSYSUP[is].c_str()) == 0);
	     if( isDown || isUp )
	       {
		  sysIdx = is;
		  foundSYSVAR = 1;
		  break;
	       }	     
	  }
	if( !foundSYSVAR )
	  {
	     std::cout << "Systematic variation should correspond to one of the existing systematics" << std::endl;
	     reset();
	     exit(1);
	  }	
     }   
   
   if( *nT > NTEMPLATEMAX )
     {
	std::cout << "Max number of templates is " << NTEMPLATEMAX << std::endl;
	reset();
	exit(1);
     }   
   if( nSYS > SYSMAX )
     {
	std::cout << "Max number of systematic variations is " << SYSMAX << std::endl;
	reset();
	exit(1);
     }
   
   int nbins = h_data->GetXaxis()->GetNbins();
   if( nbins > NBINMAX )
     {
	std::cout << "Max number of bins in template is " << NBINMAX << std::endl;
	reset();
	exit(1);
     }   

   sortTemplateNames();

   for(int it=0;it<*nT;it++)
     {	
	std::string hnameNOM = nameT_CURRENT[it];
	if( sysIdx < 666 && nameSYSVAR != "" )
	  {
	     hnameNOM = nameT_CURRENT[it]+nameSYSVAR;
	  }
	TH1D *htNOM = (TH1D*)dfile->Get(hnameNOM.c_str());
	if( !htNOM )
	  {
	     std::cout << "Histogram does not exist: " << hnameNOM << std::endl;
	     reset();
	     exit(1);
	  }	
	if( htNOM->Integral() == 0. )
	  {
	     std::cout << "Template histogram is empty" << std::endl;
	     reset();
	     exit(1);
	  }	
	std::string hcopyNOMTRUE = nameT_CURRENT[it]+"_NOMTRUEcopy";
	if( hnameNOM != nameT_CURRENT[it] )
	  {	     
	     TH1D *htNOMTRUE = (TH1D*)dfile->Get(nameT_CURRENT[it].c_str());
	     histNOMTRUE.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htNOMTRUE->Clone(hcopyNOMTRUE.c_str()))) ));
	     delete htNOMTRUE;
	  }
	else
	  {
	     histNOMTRUE.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htNOM->Clone(hcopyNOMTRUE.c_str()))) ));
	  }	
	std::string hcopyNOM = hnameNOM+"_NOMcopy";
	histNOM.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htNOM->Clone(hcopyNOM.c_str()))) ));
	delete htNOM;
	     
	if( nTAG > 0 )
	  {	     
	     std::string hnameTAGNOM = nameTAG[it];
	     if( sysIdx < 666 && nameSYSVAR != "" )
	       {
		  hnameTAGNOM = nameTAG[it]+nameSYSVAR;
	       }
	     
	     TH1D *htTAGNOM = (TH1D*)dfile->Get(hnameTAGNOM.c_str());
	     if( !htTAGNOM )
	       {
		  std::cout << "Histogram does not exist: " << hnameTAGNOM << std::endl;
		  reset();
		  exit(1);
	       }	
	     if( htTAGNOM->Integral() == 0. )
	       {
		  std::cout << "Template tag histogram is empty" << std::endl;
		  reset();
		  exit(1);
	       }	
	     std::string hcopyTAGNOM = hnameTAGNOM+"_TAGNOMcopy";
	     histTAGNOM.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htTAGNOM->Clone(hcopyTAGNOM.c_str()))) ));
	     delete htTAGNOM;

	     std::string hnameUNTAGNOM = nameUNTAG[it];
	     if( sysIdx < 666 && nameSYSVAR != "" )
	       {
		  hnameUNTAGNOM = nameUNTAG[it]+nameSYSVAR;
	       }
	     TH1D *htUNTAGNOM = (TH1D*)dfile->Get(hnameUNTAGNOM.c_str());
	     if( !htUNTAGNOM )
	       {
		  if( hnameUNTAGNOM == "" ) std::cout << "Please also specify the un-tag histogram" << std::endl;
		  else std::cout << "Histogram does not exist: " << hnameUNTAGNOM << std::endl;
		  reset();
		  exit(1);
	       }	
	     if( htUNTAGNOM->Integral() == 0. )
	       {
		  std::cout << "Template untag histogram is empty" << std::endl;
		  reset();
		  exit(1);
	       }
	     std::string hcopyUNTAGNOM = hnameUNTAGNOM+"_UNTAGNOMcopy";
	     histUNTAGNOM.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htUNTAGNOM->Clone(hcopyUNTAGNOM.c_str()))) ));
	     delete htUNTAGNOM;

	     int nbNom = histNOM[it]->GetXaxis()->GetNbins();
	     int nbTag = histTAGNOM[it]->GetXaxis()->GetNbins();
	     int nbUntag = histUNTAGNOM[it]->GetXaxis()->GetNbins();
	     float xMinNom = histNOM[it]->GetXaxis()->GetBinLowEdge(1);
	     float xMaxNom = histNOM[it]->GetXaxis()->GetBinUpEdge(nbNom);
	     float xMinTag = histTAGNOM[it]->GetXaxis()->GetBinLowEdge(1);
	     float xMaxTag = histTAGNOM[it]->GetXaxis()->GetBinUpEdge(nbTag);
	     float xMinUntag = histUNTAGNOM[it]->GetXaxis()->GetBinLowEdge(1);
	     float xMaxUntag = histUNTAGNOM[it]->GetXaxis()->GetBinUpEdge(nbUntag);
	     if( nbNom != nbTag )
	       {
		  std::cout << "Pre-tag and tag histograms should have the same number of bins" << std::endl;
		  reset();
		  exit(1);
	       }	
	     if( nbNom != nbUntag )
	       {
		  std::cout << "Pre-tag and un-tag histograms should have the same number of bins" << std::endl;
		  reset();
		  exit(1);
	       }	
	     if( xMinNom != xMinTag || xMaxNom != xMaxTag )
	       {
		  std::cout << "Pre-tag and tag histograms should have the same x-axis ranges" << std::endl;
		  reset();
		  exit(1);
	       }		     
	     if( xMinNom != xMinUntag || xMaxNom != xMaxUntag )
	       {
		  std::cout << "Pre-tag and un-tag histograms should have the same x-axis ranges" << std::endl;
		  reset();
		  exit(1);
	       }		     
	  }	

	double err = 0.;
	_NMC[it] = histNOMTRUE[it]->IntegralAndError(0,histNOM[it]->GetXaxis()->GetNbins()+1,err);
	_NMCERR[it] = err;
	
	for(int is=0;is<nSYS;is++)
	  {
	     std::string hnameDOWN = nameT_CURRENT[it]+nameSYSDOWN[is];
	     if( nameSYSDOWN[is] == nameSYSVAR ) hnameDOWN = nameT_CURRENT[it];
	     
	     TH1D *htDOWN = (TH1D*)dfile->Get(hnameDOWN.c_str());
	     if( !htDOWN )
	       {
		  std::cout << "Histogram does not exist: " << hnameDOWN << std::endl;
		  reset();
		  exit(1);
	       }	
	     std::string hcopyDOWN = hnameDOWN+"_DOWNcopy";
	     histDOWN.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htDOWN->Clone(hcopyDOWN.c_str()))) ));
	     delete htDOWN;	     

	     std::string hnameUP = nameT_CURRENT[it]+nameSYSUP[is];
	     if( nameSYSUP[is] == nameSYSVAR ) hnameUP = nameT_CURRENT[it];
	     
	     TH1D *htUP = (TH1D*)dfile->Get(hnameUP.c_str());
	     if( !htUP )
	       {
		  std::cout << "Histogram does not exist: " << hnameUP << std::endl;
		  reset();
		  exit(1);
	       }	
	     std::string hcopyUP = hnameUP+"_UPcopy";
	     histUP.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htUP->Clone(hcopyUP.c_str()))) ));
	     delete htUP;
	  }
	
	if( nSYS == 0 )
	  {
	     std::string hnameDOWN = nameT_CURRENT[it];
	     TH1D *htDOWN = (TH1D*)dfile->Get(hnameDOWN.c_str());
	     if( !htDOWN )
	       {
		  std::cout << "Histogram does not exist: " << hnameDOWN << std::endl;
		  reset();
		  exit(1);
	       }	
	     std::string hcopyDOWN = hnameDOWN+"_DOWNcopy";
	     histDOWN.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htDOWN->Clone(hcopyDOWN.c_str()))) ));
	     delete htDOWN;	     

	     std::string hnameUP = nameT_CURRENT[it];
	     TH1D *htUP = (TH1D*)dfile->Get(hnameUP.c_str());
	     if( !htUP )
	       {
		  std::cout << "Histogram does not exist: " << hnameUP << std::endl;
		  reset();
		  exit(1);
	       }	
	     std::string hcopyUP = hnameUP+"_UPcopy";
	     histUP.push_back(std::unique_ptr<TH1D>( (static_cast<TH1D*>(htUP->Clone(hcopyUP.c_str()))) ));
	     delete htUP;
	  }
     }
   
   for(int is=0;is<nSYS;is++)
     {
	for(int i=0;i<NBINMAX;i++)
	  {
	     binFactor[is][i] = 1.;
	     
	     for(int j=0;j<NBINMAX;j++)
	       {
		  binFactor2D[is][i][j] = 1.;
	       }	     
	  }	
     }

   // statistical variation of templates for SF study

   if( sysIdx >= 666 )
     {
	for(int it=0;it<*nT;it++)
	  {
	     TRandom3 *rt = new TRandom3(sysIdx+it);
	     int nbins = histNOM[it]->GetXaxis()->GetNbins();
	   
	     if( nTAG > 0 )
	       {		  
		  double integ_nom = histNOM[it]->Integral(0,nbins+1);
		  double integ_tag = histTAGNOM[it]->Integral(0,nbins+1);
		  double integ_untag = histUNTAGNOM[it]->Integral(0,nbins+1);
	     
		  for(int ibb=0;ibb<=nbins+1;ibb++)
		    {
		       double v = histNOM[it]->GetBinContent(ibb);
		       double e = histNOM[it]->GetBinError(ibb);
		       double v_t = histTAGNOM[it]->GetBinContent(ibb);
		       double e_t = histTAGNOM[it]->GetBinError(ibb);
		       double v_nt = histUNTAGNOM[it]->GetBinContent(ibb);
		       double e_nt = histUNTAGNOM[it]->GetBinError(ibb);
		       double rnd = rt->Gaus(0,1);
		       double rnd2 = rt->Gaus(0,1);
		       double vNew = v;
		       
		       if( integ_nom/integ_tag-1 < 10E-5 )
			 {		       
			    vNew += e_t*rnd;
			 }		  
		       else
			 {		       
			    vNew += e_t*rnd+e_nt*rnd2;
			 }		  

		       histNOM[it]->SetBinContent(ibb,vNew);
		    }
	       }	     
	     else
	       {
		  double integ_nom = histNOM[it]->Integral(0,nbins+1);
		  
		  for(int ibb=0;ibb<=nbins+1;ibb++)
		    {
		       double v = histNOM[it]->GetBinContent(ibb);
		       double e = histNOM[it]->GetBinError(ibb);
		       double rnd = rt->Gaus(0,1);
		       double vNew = v;
		       
		       vNew += e*rnd;
		       
		       histNOM[it]->SetBinContent(ibb,vNew);		       
		    }		  
	       }	
	     delete rt;
	  }
	
	TRandom3 *rt = new TRandom3(sysIdx+10000);
	int nbins = h_data->GetXaxis()->GetNbins();

	if( nTAG > 0 )
	  {	     
	     double integ_nom = h_data->Integral(0,nbins+1);
	     double integ_tag = h_data_tag->Integral(0,nbins+1);
	     double integ_untag = h_data_untag->Integral(0,nbins+1);
	     
	     for(int ibb=0;ibb<=nbins+1;ibb++)
	       {
		  double v = h_data->GetBinContent(ibb);
		  double e = h_data->GetBinError(ibb);
		  double v_t = h_data_tag->GetBinContent(ibb);
		  double e_t = h_data_tag->GetBinError(ibb);
		  double v_nt = h_data_untag->GetBinContent(ibb);
		  double e_nt = h_data_untag->GetBinError(ibb);
		  double rnd = rt->Gaus(0,1);
		  double rnd2 = rt->Gaus(0,1);
		  double vNew = v;
		  
		  if( integ_nom/integ_tag-1 < 10E-5 )
		    {
		       vNew += e_t*rnd;
		    }	     
		  else
		    {		  
		       vNew += e_t*rnd+e_nt*rnd2;
		    }	     
		  
		  h_data->SetBinContent(ibb,vNew);
	       }
	  }	
	else
	  {
	     double integ_nom = h_data->Integral(0,nbins+1);
	     
	     for(int ibb=0;ibb<=nbins+1;ibb++)
	       {
		  double v = h_data->GetBinContent(ibb);
		  double e = h_data->GetBinError(ibb);
		  double rnd = rt->Gaus(0,1);
		  double vNew = v;
		  
		  vNew += e*rnd;
		  
		  h_data->SetBinContent(ibb,vNew);
	       }	     
	  }	     
	delete rt;	
     }
   
   // normalise all templates to give unity when summed up   
   if( doNormalise )
     {	
	float totSum = 0.;
	for(int i=0;i<*nT;i++)
	  {
	     totSum += histNOMTRUE[i]->Integral(0,histNOMTRUE[i]->GetXaxis()->GetNbins()+1);
	  }   
	float k = (totSum > 0.) ? 1./totSum : 1.;
	for(int i=0;i<*nT;i++)
	  {
	     float kSys = histNOMTRUE[i]->Integral(0,histNOMTRUE[i]->GetXaxis()->GetNbins()+1)/
	       histNOM[i]->Integral(0,histNOMTRUE[i]->GetXaxis()->GetNbins()+1);
	     histNOM[i]->Scale(k*kSys);
	     histNOMTRUE[i]->Scale(k);
	  }
     }
   
   // bin optimisation
   
   for(int it=0;it<*nT;it++)
     {	
	addbin(histNOM[it].get());
	addbin(histNOMTRUE[it].get());
	
	for(int is=0;is<nSYS;is++)
	  {
	     addbin(histDOWN[is+it*nSYS].get());
	     addbin(histUP[is+it*nSYS].get());
	
	     if( doRescaleSys ) rescale(histNOM[it].get(),histUP[is+it*nSYS].get(),histDOWN[is+it*nSYS].get());
	  }
	
	if( nSYS == 0 )
	  {
	     addbin(histDOWN[it].get());
	     addbin(histUP[it].get());	     

	     if( doRescaleSys ) rescale(histNOM[it].get(),histUP[it].get(),histDOWN[it].get());
	  }	
     }  

   addbin(h_data.get());
      
   // normalise data histogram to unity
   if( doNormalise ) h_data->Scale(1./h_data->Integral());

   if( *nT > 0 )
     {	
	TH1D *hadj = (TH1D*)histNOMTRUE[0]->Clone("hadj");
	for(int it=1;it<*nT;it++)
	  {	
	     hadj->Add(histNOMTRUE[it].get(),1.);
	  }
//	TH1D *hadj = (TH1D*)histNOM[0]->Clone("hadj");
//	for(int it=1;it<*nT;it++)
//	  {	
//	     hadj->Add(histNOM[it].get(),1.);
//	  }
	
	double vIntData = h_data->Integral();
	double ibinMaxData = ibinMax;
    std::cout << ibinMaxData <<std::endl;
	if( vIntData > 0. )
	  {
	     int nb = h_data->GetXaxis()->GetNbins();
	     double vIntCur = 0.;
	     for(int ib=0;ib<=nb+1;ib++)
	       {
		  double v = h_data->GetBinContent(ib);
		  vIntCur += v;
		  if( vIntCur/vIntData >= (1.-0.02) )
		    {
		       ibinMaxData = ib;
		       break;
		    }		  
	       }	     
	  }	
	
	double vIntMC = hadj->Integral();
	double ibinMaxMC = ibinMax;
	if( vIntMC > 0. )
	  {
	     int nb = hadj->GetXaxis()->GetNbins();
	     double vIntCur = 0.;
	     for(int ib=0;ib<=nb+1;ib++)
	       {
		  double v = hadj->GetBinContent(ib);
		  vIntCur += v;
		  if( vIntCur/vIntMC >= (1.-0.02) )
		    {
		       ibinMaxMC = ib;
		       break;
		    }		  
	       }	     
	  }
	
	double ibinMinMC = ibinMin;
	if( vIntMC > 0. )
	  {
	     int nb = hadj->GetXaxis()->GetNbins();
	     double vIntCur = 0.;
	     for(int ib=nb+1;ib>=0;ib--)
	       {
		  double v = hadj->GetBinContent(ib);
		  vIntCur += v;
		  if( (vIntMC-vIntCur)/vIntMC <= 0.01 )
		    {
		       ibinMinMC = ib;
		       break;
		    }		  
	       }	     
	  }
	
	ibinMax = ibinMaxMC;
	ibinMin = ibinMinMC;
	
	if( *verb ) std::cout << "Fit is done in the range of bins #" << ibinMin << "-" << ibinMax << std::endl;
	
	// adjust empty bins in data
//	adjust(h_data.get());

	for(int it=0;it<*nT;it++)
	  {	
	     removeMCNeg(histNOM[it].get());

	     for(int is=0;is<nSYS;is++)
	       {
		  removeMCNeg(histDOWN[is+it*nSYS].get());
		  removeMCNeg(histUP[is+it*nSYS].get());
	       }
	     
	     if( nSYS == 0 )
	       {
		  removeMCNeg(histDOWN[it].get());
		  removeMCNeg(histUP[it].get());	     
	       }	
	  }  	
	
	// merge last bins
	double *xn;
	int nb0 = histNOM[0]->GetXaxis()->GetNbins();
	bool doRebin = ( ibinMax < nb0 || ibinMin > 1 );
	const int arrs = ibinMax-(ibinMin-1)+1;
	double xnew[arrs];

	xn = mergeBins(h_data.get(),ibinMin,ibinMax,xnew);
	if( doRebin ) h_data = std::unique_ptr<TH1D>( (static_cast<TH1D*>(h_data->Rebin(arrs-1,"",xn))) );
	
	for(int it=0;it<*nT;it++)
	  {		
	     xn = mergeBins(histNOM[it].get(),ibinMin,ibinMax,xnew);
	     if( doRebin ) histNOM[it] = std::unique_ptr<TH1D>( (static_cast<TH1D*>(histNOM[it]->Rebin(arrs-1,"",xn))) );

	     for(int is=0;is<nSYS;is++)
	       {
		  xn = mergeBins(histDOWN[is+it*nSYS].get(),ibinMin,ibinMax,xnew);
		  if( doRebin ) histDOWN[is+it*nSYS] = std::unique_ptr<TH1D>( (static_cast<TH1D*>(histDOWN[is+it*nSYS]->Rebin(arrs-1,"",xn))) );

		  xn = mergeBins(histUP[is+it*nSYS].get(),ibinMin,ibinMax,xnew);
		  if( doRebin ) histUP[is+it*nSYS] = std::unique_ptr<TH1D>( (static_cast<TH1D*>(histUP[is+it*nSYS]->Rebin(arrs-1,"",xn))) );
	       }	     
	     
	     if( nSYS == 0 )
	       {
		  xn = mergeBins(histDOWN[it].get(),ibinMin,ibinMax,xnew);
		  if( doRebin ) histDOWN[it] = std::unique_ptr<TH1D>( (static_cast<TH1D*>(histDOWN[it]->Rebin(arrs-1,"",xn))) );

		  xn = mergeBins(histUP[it].get(),ibinMin,ibinMax,xnew);
		  if( doRebin ) histUP[it] = std::unique_ptr<TH1D>( (static_cast<TH1D*>(histUP[it]->Rebin(arrs-1,"",xn))) );
	       }
	  }

	*nBINS = h_data->GetXaxis()->GetNbins();
     }
   
   if( *verb ) std::cout << "Initialisation finished" << std::endl;   
}

// Destructor
CFIT::cfit::~cfit()
{
   reset();
   resetNames();
   if( dfile ) dfile->Close();
   delete dfile;
   if( fcov ) delete fcov;
}

// chi2 term calculation
double funcDEFM(double vdata1,
		double vdataErr1,
		std::vector<double> vmc1,
		double vdata2,
		double vdataErr2,
		std::vector<double> vmc2,
		double covEL,
		double norm1,
		double norm2,
		double *par)
{
   double term1 = 0.;
   double term2 = 0.;
   
   for(unsigned int i=0;i<vmc1.size();i++)
     {
	term1 += par[CFIT::cfit::PARIDXFIT->at(i)]*vmc1[i];
	term2 += par[CFIT::cfit::PARIDXFIT->at(i)]*vmc2[i];
     }
   
   double val = ((vdata1-term1)/norm1)*((vdata2-term2)/norm2)*covEL;
   
   return val;
}

// chi2 definition
void fcnSysM(int &npar, double *gin, double &f, double *par, int iflag)
{
   double chisq = 0;
   double delta;

   double nb = CFIT::cfit::vecM->size();
   
   for( int i1=0;i1<nb;i1++ )
     {
	int idx1 = CFIT::cfit::vecM->at(i1);
	
	std::vector<double> vmc1;
	vmc1.clear();
	for(int it1=0;it1<*CFIT::cfit::nT;it1++)
	  {
	     vmc1.push_back(CFIT::cfit::histNOM[it1]->GetBinContent(idx1));
	  }	
	
	double data_1 = CFIT::cfit::h_data->GetBinContent(idx1);
	double data_err_1 = CFIT::cfit::h_data->GetBinError(idx1);
	
	for( int i2=0;i2<nb;i2++ )
	  {
	     int idx2 = CFIT::cfit::vecM->at(i2);

	     std::vector<double> vmc2;
	     vmc2.clear();
	     for(int it2=0;it2<*CFIT::cfit::nT;it2++)
	       {
		  vmc2.push_back(CFIT::cfit::histNOM[it2]->GetBinContent(idx2));
	       }		     
	     
	     double data_2 = CFIT::cfit::h_data->GetBinContent(idx2);
	     double data_err_2 = CFIT::cfit::h_data->GetBinError(idx2);
	     
	     delta = funcDEFM(data_1,
			      data_err_1,
			      vmc1,
			      data_2,
			      data_err_2,
			      vmc2,
			      (*CFIT::cfit::covMIp)[i1][i2],
			      (*CFIT::cfit::norm1Dp)[i1],
			      (*CFIT::cfit::norm1Dp)[i2],
			      par);

	     chisq += delta;
	  }
     }
   
//   if( iflag == 1 )
//     {		
//	if( *verb ) std::cout << "CHISQ/NDOF (NO FIT) = "<< chisq/(nBINS-nT) << std::endl;
//     }
   if( iflag == 3 )
     {
	*CFIT::cfit::VALID = 1;
	*CFIT::cfit::CHISQ = chisq/(*CFIT::cfit::nBINS-*CFIT::cfit::nT);
	*CFIT::cfit::NDOF = *CFIT::cfit::nBINS-*CFIT::cfit::nT;
	if( *CFIT::cfit::verb ) std::cout << "CHISQ/NDOF = " << *CFIT::cfit::CHISQ << std::endl;
     }
   
   f = chisq;
}

// measure max and min systematic variations
void CFIT::cfit::totSys(TH1D *h_nom,TH1D *h_down,TH1D *h_up)
{   
   int nbins = h_nom->GetXaxis()->GetNbins();
   
   for(int ib=1;ib<=nbins;ib++)
     {
	double b_nom = h_nom->GetBinContent(ib);
	double b_cur_down = h_down->GetBinContent(ib);
	double b_cur_up = h_up->GetBinContent(ib);
	
	double del_down = b_cur_down - b_nom;
	double del_up = b_cur_up - b_nom;
	double del_down_res = del_down;
	double del_up_res = del_up;
	
	double del_min = (del_down < del_up) ? del_down : del_up;
	double del_max = (del_down > del_up) ? del_down : del_up;

	del_down_res = (del_min > 0.) ? 0. : del_min;
	del_up_res = (del_max < 0.) ? 0. : del_max;
	
	double sys_down = b_nom + del_down_res;
	double sys_up = b_nom + del_up_res;
	
	h_down->SetBinContent(ib,sys_down);
	h_up->SetBinContent(ib,sys_up);
     }
}

// quadratic combination of systematics 
void CFIT::cfit::combSys(TH1D *h_nom,TH1D** h_sys_down,TH1D** h_sys_up,TH1D *h_sys_down_comb,TH1D *h_sys_up_comb)
{
   double sys_up[NBINMAX];
   double sys_down[NBINMAX];
   
   for(int is=1;is<=*nBINS;is++)
     {
	sys_up[is-1] = 0.;
	sys_down[is-1] = 0.;
     }      
   
   for(int j=0;j<nSYS;j++)
     {
	for(int ib=1;ib<=*nBINS;ib++)
	  {
	     double nom = h_nom->GetBinContent(ib);
	     double down = h_sys_down[j]->GetBinContent(ib);
	     double up = h_sys_up[j]->GetBinContent(ib);
	     double delta_down = nom - down;
	     double delta_up = up - nom;
	     sys_up[ib-1] = sqrt(pow(sys_up[ib-1],2) + pow(delta_up,2));
	     sys_down[ib-1] = sqrt(pow(sys_down[ib-1],2) + pow(delta_down,2));
	  }
     }      
   
   for(int ib=1;ib<=*nBINS;ib++)
     {
	double nom = h_nom->GetBinContent(ib);
	double err = h_nom->GetBinError(ib);
	
	sys_up[ib-1] = sqrt(pow(sys_up[ib-1],2) + pow(err,2));
	sys_down[ib-1] = sqrt(pow(sys_down[ib-1],2) + pow(err,2));
     }   
   
   for(int ib=1;ib<=*nBINS;ib++)
     {		
	double nom = h_nom->GetBinContent(ib);
	h_sys_down_comb->SetBinContent(ib,nom-sys_down[ib-1]);
	h_sys_up_comb->SetBinContent(ib,nom+sys_up[ib-1]);
     }
}

// linear combination of systematics
void CFIT::cfit::combSysLinear(std::vector<TH1D*> h_nom,
			       std::vector<TH1D*> h_sys_up,
			       std::vector<TH1D*> h_sys_down,
			       TH1D *h_comb,
			       TH1D *h_sys_down_comb,
			       TH1D *h_sys_up_comb)
{
   double sys_up[NBINMAX];
   double sys_down[NBINMAX];
   double stat[NBINMAX];
   
   for(int is=1;is<=*nBINS;is++)
     {
	sys_up[is-1] = 0.;
	sys_down[is-1] = 0.;
	stat[is-1] = 0.;
     }

   for(int ib=1;ib<=*nBINS;ib++)
     {
	double delta_down = 0.;
	double delta_up = 0.;
	
	for(unsigned int i=0;i<h_nom.size();i++)
	  {
	     double nom = h_nom[i]->GetBinContent(ib);
	     double up = h_sys_up[i]->GetBinContent(ib);
	     double down = h_sys_down[i]->GetBinContent(ib);
	     double err = h_nom[i]->GetBinError(ib);
	     
	     delta_down += (nom-down);
	     delta_up += (up-nom);
	     
	     stat[ib-1] += pow(err,2);
	  }	
	
	sys_up[ib-1] = sys_up[ib-1] + delta_up;
	sys_down[ib-1] = sys_down[ib-1] + delta_down;

	stat[ib-1] = sqrt(stat[ib-1]);
     }      
   
   for(int ib=1;ib<=*nBINS;ib++)
     {
	double nom = 0.;
	for(unsigned int i=0;i<h_nom.size();i++)
	  {
	     nom += h_nom[i]->GetBinContent(ib);
	  }	

	h_comb->SetBinContent(ib,nom);
	h_comb->SetBinError(ib,stat[ib-1]);
	h_sys_down_comb->SetBinContent(ib,nom-sys_down[ib-1]);
	h_sys_up_comb->SetBinContent(ib,nom+sys_up[ib-1]);
     }
}

// analysis of systematic uncertainties and correlation matrix construction
void CFIT::cfit::applySys(std::vector<double> sfl,
			  std::vector<double> sferrl,
			  bool write,
			  std::string option)
{
   std::vector<TH1D*> v_fit;
   
   for(int i=0;i<*nT;i++)
     {
	std::string c_fit = std::string(histNOM[i]->GetName()) + "_fit";
	TH1D *h_fit = (TH1D*)histNOM[i]->Clone(c_fit.c_str());
	h_fit->Scale(sfl[i]);
	v_fit.push_back(h_fit);
     }   

   std::vector<std::vector<TH1D*> > v_sys_down_fit;
   std::vector<std::vector<TH1D*> > v_sys_up_fit;

   for(int i=0;i<nSYS;i++)
     {
	std::vector<TH1D*> v_sys_down;
	std::vector<TH1D*> v_sys_up;
	
	for(int j=0;j<*nT;j++)
	  {
	     std::string hname_down = histDOWN[i+nSYS*j]->GetName();
	     std::string hname_up = histUP[i+nSYS*j]->GetName();
	     
	     TH1D *h_down = (TH1D*)histDOWN[i+nSYS*j]->Clone((hname_down+"_fit").c_str());
	     TH1D *h_up = (TH1D*)histUP[i+nSYS*j]->Clone((hname_up+"_fit").c_str());
	
	     h_down->Scale(sfl[j]);
	     h_up->Scale(sfl[j]);
	     
	     v_sys_down.push_back(h_down);
	     v_sys_up.push_back(h_up);
	  }
	
	v_sys_down_fit.push_back(v_sys_down);
	v_sys_up_fit.push_back(v_sys_up);
     }

   if( nSYS == 0 )
     {
	std::vector<TH1D*> v_sys_down;
	std::vector<TH1D*> v_sys_up;
	
	for(int j=0;j<*nT;j++)
	  {
	     std::string hname_down = histDOWN[j]->GetName();
	     std::string hname_up = histUP[j]->GetName();

	     TH1D *h_down = (TH1D*)histDOWN[j]->Clone((hname_down+"_fit").c_str());
	     TH1D *h_up = (TH1D*)histUP[j]->Clone((hname_up+"_fit").c_str());
	
	     h_down->Scale(sfl[j]);
	     h_up->Scale(sfl[j]);
	     
	     v_sys_down.push_back(h_down);
	     v_sys_up.push_back(h_up);
	  }
	
	v_sys_down_fit.push_back(v_sys_down);
	v_sys_up_fit.push_back(v_sys_up);
     }

   if( !h_comb ) delete h_comb;
   h_comb = (TH1D*)histNOM[0]->Clone("h_comb");
   if( !h_sys_down_combSum ) delete h_sys_down_combSum;
   h_sys_down_combSum = (TH1D*)histNOM[0]->Clone("h_sys_down_combSum");
   if( !h_sys_up_combSum ) delete h_sys_up_combSum;
   h_sys_up_combSum = (TH1D*)histNOM[0]->Clone("h_sys_up_combSum");
   
   for(int j=0;j<nSYS;j++)
     {
	std::string hnameDOWN = "h_sys_down_comb_"+nameSYS[j];
	if( !h_sys_down_comb[j] ) delete h_sys_down_comb[j];
	h_sys_down_comb[j] = (TH1D*)h_comb->Clone(hnameDOWN.c_str());

	std::string hnameUP = "h_sys_up_comb_"+nameSYS[j];
	if( !h_sys_up_comb[j] ) delete h_sys_up_comb[j];
	h_sys_up_comb[j] = (TH1D*)h_comb->Clone(hnameUP.c_str());
     }   

   if( nSYS == 0 )
     {
	std::string hnameDOWN = "h_sys_down_comb";
	if( !h_sys_down_comb[0] ) delete h_sys_down_comb[0];
	h_sys_down_comb[0] = (TH1D*)h_comb->Clone(hnameDOWN.c_str());

	std::string hnameUP = "h_sys_up_comb";
	if( !h_sys_up_comb[0] ) delete h_sys_up_comb[0];
	h_sys_up_comb[0] = (TH1D*)h_comb->Clone(hnameUP.c_str());
     }   

   for(int j=0;j<nSYS;j++)
     {
	combSysLinear(v_fit,
		      v_sys_up_fit[j],
		      v_sys_down_fit[j],
		      h_comb,
		      h_sys_down_comb[j],
		      h_sys_up_comb[j]);
	
	std::string nat_down = "h_sys_down_combNATURAL_"+nameSYSDOWN[j];
	if( !h_sys_down_combNATURAL[j] ) delete h_sys_down_combNATURAL[j];
	h_sys_down_combNATURAL[j] = (TH1D*)h_sys_down_comb[j]->Clone(nat_down.c_str());
	std::string nat_up = "h_sys_up_combNATURAL_"+nameSYSUP[j];
	if( !h_sys_up_combNATURAL[j] ) delete h_sys_up_combNATURAL[j];
	h_sys_up_combNATURAL[j] = (TH1D*)h_sys_up_comb[j]->Clone(nat_up.c_str());

	if( !h_combNATURAL ) delete h_combNATURAL;
	h_combNATURAL = (TH1D*)h_comb->Clone("h_combNATURAL");

	doFracSys(h_comb,h_sys_down_combNATURAL[j],h_sys_up_combNATURAL[j],j,option,"");
	
	totSys(h_comb,
	       h_sys_down_comb[j],
	       h_sys_up_comb[j]
	      );
	
	removeMCSys(h_comb,
		    h_sys_up_comb[j],
		    h_sys_down_comb[j]);
     }

   if( nSYS == 0 )
     {
	combSysLinear(v_fit,
		      v_sys_up_fit[0],
		      v_sys_down_fit[0],
		      h_comb,
		      h_sys_down_comb[0],
		      h_sys_up_comb[0]);

	std::string nat_down = "h_sys_down_combNATURAL";
	if( !h_sys_down_combNATURAL[0] ) delete h_sys_down_combNATURAL[0];
	h_sys_down_combNATURAL[0] = (TH1D*)h_sys_down_comb[0]->Clone(nat_down.c_str());
	std::string nat_up = "h_sys_up_combNATURAL";
	if( !h_sys_up_combNATURAL[0] ) delete h_sys_up_combNATURAL[0];
	h_sys_up_combNATURAL[0] = (TH1D*)h_sys_up_comb[0]->Clone(nat_up.c_str());

	if( !h_combNATURAL ) delete h_combNATURAL;
	h_combNATURAL = (TH1D*)h_comb->Clone("h_combNATURAL");

	doFracSys(h_comb,h_sys_down_combNATURAL[0],h_sys_up_combNATURAL[0],0,option,"");
	
	totSys(h_comb,
	       h_sys_down_comb[0],
	       h_sys_up_comb[0]
	      );
	
	removeMCSys(h_comb,
		    h_sys_up_comb[0],
		    h_sys_down_comb[0]);
     }

   combSys(h_comb,
	   h_sys_down_comb,
	   h_sys_up_comb,
	   h_sys_down_combSum,
	   h_sys_up_combSum);
   
   totSys(h_comb,
	  h_sys_down_combSum,
	  h_sys_up_combSum);

   removeMCSys(h_comb,
	       h_sys_up_combSum,
	       h_sys_down_combSum);
   
   vecM->clear();
   
   for(int i=1;i<=*nBINS;i++)
     {
	double cont = h_comb->GetBinContent(i);
	double contUp = h_sys_up_combSum->GetBinContent(i);
	double contDown = h_sys_down_combSum->GetBinContent(i);
	
	vecM->push_back(i);
     } 
   
   const int vecM_n = vecM->size();

   if( ! rcov )
     {	
	double norm[NBINMAX][NBINMAX];
	
	for(int i1=0;i1<NBINMAX;i1++)
	  {
	     for(int i2=0;i2<NBINMAX;i2++)
	       {
		  norm[i1][i2] = 0.;
	       }	
	  }   

	for(int i1=0;i1<NBINMAX;i1++)
	  {
	     for(int i2=0;i2<NBINMAX;i2++)
	       {
		  cov[i1][i2] = 0.;
	       }
	  }

	double bvMax = 10E-10;
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     int idx1 = vecM->at(i1);
	     
	     double v_nom_i1 = h_combNATURAL->GetBinContent(idx1);
	     if( v_nom_i1 > bvMax ) bvMax = v_nom_i1;
	  }	
	
	// adjust systematics
	TH1D *h_SYM = (TH1D*)h_combNATURAL->Clone("h_SYM");
	TH1D *h_SYM_UP[nSYS];
	TH1D *h_SYM_DOWN[nSYS];
	for(int is=0;is<nSYS;is++)
	  {
	     std::string h_sym_up_name = "h_SYM_UP_"+std::string(Form("%d",is));
	     std::string h_sym_down_name = "h_SYM_DOWN_"+std::string(Form("%d",is));
	     h_SYM_UP[is] = (TH1D*)h_sys_up_combNATURAL[is]->Clone(h_sym_up_name.c_str());
	     h_SYM_DOWN[is] = (TH1D*)h_sys_down_combNATURAL[is]->Clone(h_sym_down_name.c_str());	     
	     
	     adjustSys(h_SYM,h_SYM_DOWN[is],h_SYM_UP[is],is);
	     
	     doFracSys(h_SYM,h_SYM_DOWN[is],h_SYM_UP[is],is,option,"Optimised");
	  }

	// compute covariance matrix
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     int idx1 = vecM->at(i1);
	     
	     double siglin = 0.;
	     
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  int idx2 = vecM->at(i2);
		  
		  double cov_v = 0.;
		  double norm_i1_v = 0.;
		  double norm_i2_v = 0.;

		  for(int is=-1;is<nSYS;is++)
		    {
		       TH1D *nom = h_SYM;
		       TH1D *down = h_SYM;
		       TH1D *up = h_SYM;
		       TH1D *downSum = h_sys_down_combSum;
		       TH1D *upSum = h_sys_up_combSum;

		       if( is != -1 )
			 {
			    down = h_SYM_DOWN[is];
			    up = h_SYM_UP[is];
			 }
		       
		       double v_nom_i1 = nom->GetBinContent(idx1);
		       double v_sysUp_i1 = up->GetBinContent(idx1);
		       double v_sysDown_i1 = down->GetBinContent(idx1);
		       
		       double v_nom_i2 = nom->GetBinContent(idx2);
		       double v_sysUp_i2 = up->GetBinContent(idx2);
		       double v_sysDown_i2 = down->GetBinContent(idx2);
		       
		       double sysUp_i1 = v_sysUp_i1 - v_nom_i1;
		       double sysDown_i1 = v_nom_i1 - v_sysDown_i1;
		       
		       double sysUp_i2 = v_sysUp_i2 - v_nom_i2;
		       double sysDown_i2 = v_nom_i2 - v_sysDown_i2;
		       
		       if( *covMode == COV_BARLOW )
			 {			    
			    if( sysUp_i1*sysDown_i1 < 0 )
			      {
				 if( fabs(sysUp_i1) > fabs(sysDown_i1) )
				   sysDown_i1 = 0.;
				 else
				   sysUp_i1 = 0.;
			      }		  
			    if( sysUp_i2*sysDown_i2 < 0 )
			      {
				 if( fabs(sysUp_i2) > fabs(sysDown_i2) )
				   sysDown_i2 = 0.;
				 else
				   sysUp_i2 = 0.;
			      }		  
			 }		       
		       
		       if( i2 == 0 )
			 siglin += sysUp_i1 - sysDown_i1;
		       
		       double stat_i1 = nom->GetBinError(idx1);
		       
		       double stat_i2 = nom->GetBinError(idx2);
		       
		       double data_stat = h_data->GetBinError(idx1);
		       if( i1 != i2 || is != -1 ) data_stat = 0.;
		       
		       double sigma_i1 = (fabs(sysUp_i1) > fabs(sysDown_i1)) ? sysUp_i1 : sysDown_i1;
		       double sigma_i2 = (fabs(sysUp_i2) > fabs(sysDown_i2)) ? sysUp_i2 : sysDown_i2;

		       if( *covMode == COV_BARLOW )
			 {
			    sigma_i1 = sysUp_i1+sysDown_i1;
			    sigma_i2 = sysUp_i2+sysDown_i2;
			 }		       
		       
		       double alpha_i1 = sysUp_i1-sysDown_i1;
		       double alpha_i2 = sysUp_i2-sysDown_i2;
		       
		       sysUp_i1 = upSum->GetBinContent(idx1) - 
			 h_comb->GetBinContent(idx1);
		       sysDown_i1 = -downSum->GetBinContent(idx1) +
			 h_comb->GetBinContent(idx1);
		       
		       sysUp_i2 = upSum->GetBinContent(idx2) - 
			 h_comb->GetBinContent(idx2);
		       sysDown_i2 = -downSum->GetBinContent(idx2) +
			 h_comb->GetBinContent(idx2);
		       
		       double sigmaAv_i1 = (sysUp_i1+sysDown_i1)/2.;
		       double sigmaAv_i2 = (sysUp_i2+sysDown_i2)/2.;
		       
		       double covAS = 0.;
		       if( is != -1 ) covAS = covFactor[is]*
			 binFactor[is][idx1]*binFactor[is][idx2]*binFactor2D[is][idx1][idx2];
		       if( optMode == OPT_NOCORR ) covAS = 0.;

		       if( covAS > 1. ) covAS = 1.;
		       
		       if( i1 == i2 ) covAS = 1.;
		       
		       if( *covMode == COV_BARLOW )
			 {		   
			    cov_v += covAS*sigma_i1/2.*sigma_i2/2.+
			      covAS*alpha_i1/2.*alpha_i2/2.*(1.-2./PI)+
			      pow(data_stat,2);
			 }
		       else if( *covMode == COV_MAX )
			 {		       
			    cov_v += covAS*sigma_i1*sigma_i2+
			      pow(data_stat,2);
			 }		  
		       
		       if( i1 == i2 && is == -1 )
			 cov_v += pow(stat_i1,2);
		       
		       if( is != -1 )
			 {		       
			    stat_i1 = 0;
			    stat_i2 = 0;
			 }
		       
		       if( *covMode == COV_BARLOW )
			 {		       
			    norm_i1_v += covAS*sigma_i1/2.*sigma_i1/2.+
			      covAS*alpha_i1/2.*alpha_i1/2.*(1.-2./PI) +
			      pow(stat_i1,2)+pow(data_stat,2);
			    
			    norm_i2_v += covAS*sigma_i2/2.*sigma_i2/2.+
			      covAS*alpha_i2/2.*alpha_i2/2.*(1.-2./PI) +
			      pow(stat_i2,2)+pow(data_stat,2);
			 }
		       else if( *covMode == COV_MAX )
			 {		       
			    norm_i1_v += covAS*sigma_i1*sigma_i1+
			      pow(stat_i1,2)+pow(data_stat,2);
			    norm_i2_v += covAS*sigma_i2*sigma_i2+
			      pow(stat_i2,2)+pow(data_stat,2);
			 }
		    }	     
		  
		  norm[i1][i2] = sqrt(norm_i1_v)*sqrt(norm_i2_v);
		  
		  if( i1 == i2 && norm[i1][i2] == 0. )
		    {
		       std::cout << "Input histograms have empty bin #" << i1 << " in both data and MC templates: please adjust the binning" << std::endl;
		       reset();
		       exit(1);		       
		    }		  
		  
		  if( i1 == i2 ) 
		    {
		       norm1D[i1] = sqrt(norm_i1_v);
		    }
		  
		  cov[i1][i2] = cov_v;
	       }
	  }	
	
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  cov[i1][i2] = cov[i1][i2] / norm1D[i1] / norm1D[i2];
	       }
	  }   		
	
	covMIp->ResizeTo(vecM_n,vecM_n);
	covMp->ResizeTo(vecM_n,vecM_n);
	
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  (*covMp)(i1,i2) = cov[i1][i2];
	       }	
	  }
	
	TMatrixDSym covMSym;
	TMatrixDSym covMISym;
	covMSym.ResizeTo(vecM_n,vecM_n);
	covMISym.ResizeTo(vecM_n,vecM_n);
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  covMSym[i1][i2] = cov[i1][i2];
		  covMISym[i1][i2] = cov[i1][i2];
	       }	
	  }
	
//	covMSym.Print();
	
	// matrix inversion
	TDecompBK lu(covMSym);
	lu.Invert(covMISym);
	
//	covMISym.Print();

	TMatrixD unityM = covMSym*covMISym;

	if( fabs(unityM.Determinant()-1) > 10E-6 )
	  {
	     std::cout << "Inversion of covariance matrix failed" << std::endl;
	     std::cout << unityM.Determinant()-1 << std::endl;
	     reset();
	     exit(1);
	  }      
	
	for(int i1=0;i1<vecM_n;i1++)
	  {
	     for(int i2=0;i2<vecM_n;i2++)
	       {
		  (*covMIp)(i1,i2) = covMISym[i1][i2];
	       }	
	  }
	
	if( covMIp->Determinant() < 0 )
	  {
	     std::cout << "Determinant is negative" << std::endl;
	     reset();
	     exit(1);
	  }	

	for(int ie=0;ie<vecM_n;ie++)
	  {
	     (*norm1Dp)[ie] = norm1D[ie];
	  }

	norm1Dp->ResizeTo(vecM_n);
	
	if( covOption == "WRITE" && write )
	  {
	     if( option == "" )
	       {	
		  std::string mname = covName+".root";
		  if( fcov ) {fcov->Close(); delete fcov; fcov = NULL;}
		  fcov = new TFile(mname.c_str(),"RECREATE");
		  if( *verb ) std::cout << "Write "<< mname << std::endl;
	       }
	     else if( option == "tag" )
	       {
		  std::string mname = covName+"_tag.root";
		  if( fcov ) {fcov->Close(); delete fcov; fcov = NULL;}
		  fcov = new TFile(mname.c_str(),"RECREATE");
		  if( *verb ) std::cout << "Write "<< mname << std::endl;
	       }	     
	     else 
	       {
		  std::cout << "Option " << option << " is now known" << std::endl;
		  reset();
		  exit(1);
	       }	     

	     std::string matrixMIName = "covMI";
	     covMIp->Write(matrixMIName.c_str());
	     std::string matrixMName = "covM";
	     covMp->Write(matrixMName.c_str());
	     std::string normName = "norm";
	     norm1Dp->Write(normName.c_str());
	  }
	
	if( h_SYM ) delete h_SYM;
	for(int is=0;is<nSYS;is++)
	  {
	     if( h_SYM_UP[is] ) delete h_SYM_UP[is];
	     if( h_SYM_DOWN[is] ) delete h_SYM_DOWN[is];
	  }   	
     }   
   else
     {
	std::string matrixMIName = "covMI";
	covMIp = std::unique_ptr<TMatrixD>(static_cast<TMatrixD*>(fcov->Get(matrixMIName.c_str())));
	std::string matrixMName = "covM";
	covMp = std::unique_ptr<TMatrixD>(static_cast<TMatrixD*>(fcov->Get(matrixMName.c_str())));
	std::string normName = "norm";
	norm1Dp = std::unique_ptr<TVectorD>(static_cast<TVectorD*>(fcov->Get(normName.c_str())));
     }

   if( producePlots )
     {	
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->SetGrid();
	c1->SetTicks();
//	c1->SetLeftMargin(0.13);
//	c1->SetBottomMargin(0.13);
	c1->SetRightMargin(0.15);
	c1->SetTopMargin(0.15);
	
	gStyle->SetPalette(1,0);

	TMatrixD covM_plot(*covMp);
	
	gStyle->SetPaintTextFormat("3g");	
	
	covM_plot.Draw("COLZ");
//	covM_plot.Draw("TEXTSAME");
	TObject *obj;
	TIter next(c1->GetListOfPrimitives());
	while( (obj=next()) )
	  {
	     TH3D *h_cov = (TH3D*)obj;
	     h_cov->GetXaxis()->SetTitle("Bin number");
	     h_cov->GetYaxis()->SetTitle("Bin number");
//	     h_cov->LabelsOption("d");
	     h_cov->SetLabelOffset(0.011);
	  }	
	if( option == "" ) 
	  {
	     std::string mname = "pics/"+covName+".pdf";
	     c1->Print(mname.c_str());
	  }	
	else 
	  {
	     std::string mname = "pics/"+covName+"_tag.pdf";
	     c1->Print(mname.c_str());
	  }	
	c1->Clear();
	delete c1;
     }

   for(int i=0;i<*nT;i++)
     {
	delete v_fit[i];
     }
   v_fit.clear();
   
   for(int i=0;i<nSYS;i++)
     {
	for(int j=0;j<*nT;j++)
	  {
	     delete v_sys_down_fit[i][j];
	     delete v_sys_up_fit[i][j];
	  }
	v_sys_down_fit[i].clear();
	v_sys_up_fit[i].clear();
     }
   v_sys_down_fit.clear();
   v_sys_up_fit.clear();
}

// minuit fit
void CFIT::cfit::doFit(double *chis,std::string option)
{
   double pval[*nT], perr[*nT];
   
   TMinuit *gMinuit = new TMinuit(100);
   gMinuit->SetFCN(fcnSysM);
   gMinuit->SetPrintLevel(-1);
   
   Double_t arglist[10];
   Int_t ierflg = 0;
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
   
   arglist[0] = 5000;
   arglist[1] = 0.1;
   
   for(int i=0;i<*nTFIT;i++)
     {
	std::string pname = "p"+std::string(Form("%d",i));
	gMinuit->mnparm(i,pname.c_str(),1.,pow(10.,-4),0.00001,1000,ierflg);
     }

   // do the first fit
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

   sf.clear();
   sferr.clear();   
   for(int i=0;i<*nT;i++)
     {   
	Double_t parg, perg;
	gMinuit->GetParameter(PARIDXFIT->at(i),parg,perg);
	sf.push_back(parg);
	sferr.push_back(perg);
//	std::cout << parg[i] << "+" << perg[i] << std::endl;
     }   

   // recalculate correlation matrix
   applySys(sf,sferr,1,option);

   // refit with updated correlation matrix
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("MINOS", arglist ,2,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("MINOS", arglist ,2,ierflg);

   double paramErrP[*nT];
   double paramErrM[*nT];
   double paramErrS[*nT];
   double paramErrG[*nT];
  
   for(int i=0;i<*nT;i++)
     {	
	gMinuit->mnerrs(PARIDXFIT->at(i),
			paramErrP[PARIDXFIT->at(i)],
			paramErrM[PARIDXFIT->at(i)],
			paramErrS[PARIDXFIT->at(i)],
			paramErrG[PARIDXFIT->at(i)]);
     }   
   
   arglist[0] = 3;
   gMinuit->mnexcm("CALL FCN",arglist,1,ierflg);
   
   double amin,edm,errdef;
   int nvpar,nparx,icstat;   
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   chis[0] = amin;
   chis[1] = *nBINS - gMinuit->GetNumFreePars();

   for(int i=0;i<*nT;i++)
     {	   
	gMinuit->GetParameter(PARIDXFIT->at(i),pval[i],perr[i]);
	if( *verb ) std::cout << histNOM.at(i)->GetName() << " " << pval[i] << " +- " << perr[i] << std::endl;
     }   

   sf.clear();
   sferr.clear();
   for(int i=0;i<*nT;i++)
     {
	sf.push_back(pval[i]);
	sferr.push_back(perr[i]);
     }
   
   delete gMinuit;
}

// user main function to perform a fit
void CFIT::cfit::Run(std::string option)
{
   if( dfile && !dfile->IsOpen() )
     {
	std::cout << "Please specify input file" << std::endl;
	reset();
	exit(1);
     }
   
   processInput(option);
   
   processGlueTemplates(option);

   sf.clear();
   sferr.clear();
   
   for(int i=0;i<*nT;i++)
     {
	sf.push_back(1.);
	sferr.push_back(0.);
     }
   
   applySys(sf,sferr,0,option);
   
   doFit(chis,option);
   
   applySys(sf,sferr,0,option);
   
   applySF(option);
   
   if( drawPPF )
     {	
	drawPrePostFit(option,0);
	drawPrePostFit(option,1);
     }   
   
   // restore flucltuations
   nameSYSVAR = "";
   sysIdx = 0;
}

// analysis of fractional systematic uncertainties
void CFIT::cfit::doFracSys(TH1D *hnom,TH1D *hsysDown,TH1D *hsysUp,int isys,std::string option,std::string postfix)
{
   TH1D *hsysDownf;
   TH1D *hsysUpf;

   TH1D *hsysDownf_cov;
   TH1D *hsysUpf_cov;
   
   hsysDownf = (TH1D*)hnom->Clone("hsysDownf");
   hsysUpf = (TH1D*)hnom->Clone("hsysUpf");
   
   hsysDownf_cov = (TH1D*)hnom->Clone("hsysDownf_cov");
   hsysUpf_cov = (TH1D*)hnom->Clone("hsysUpf_cov");
   
   for(int ib0=1;ib0<=*nBINS;ib0++)
     {
	float b0 = hnom->GetBinContent(ib0);
	float b0err = hnom->GetBinError(ib0);
	float bDown = hsysDown->GetBinContent(ib0);
	float bUp = hsysUp->GetBinContent(ib0);
	float bDown_cov = bDown;
	float bUp_cov = bUp;
	
	if( b0 > 0. )
	  {	     
	     bDown = (b0-bDown)/b0;
	     bUp = (b0-bUp)/b0;
	     
	     bDown_cov = (b0-bDown_cov)/b0;
	     bUp_cov = (bUp_cov-b0)/b0;
	  }
	else
	  {
	     bDown = 0.;
	     bUp = 0.;
	     
	     bDown_cov = 0.;
	     bUp_cov = 0.;
	  }
	
	hsysDownf->SetBinContent(ib0,1.-bDown);
	hsysDownf->SetBinError(ib0,0.);
	hsysUpf->SetBinContent(ib0,1.-bUp);
	hsysUpf->SetBinError(ib0,0.);
	
	hsysDownf_cov->SetBinContent(ib0,bDown_cov);
	hsysDownf_cov->SetBinError(ib0,0.);
	hsysUpf_cov->SetBinContent(ib0,bUp_cov);
	hsysUpf_cov->SetBinError(ib0,0.);
     }	

   if( postfix != "Optimised" && nSYS > 0 )
     {		
	double corr = bbCorr(hsysDownf_cov,hsysUpf_cov);
	covFactor[isys] = corr;
     }
//   else
//     {
//	hsysDownf = (TH1D*)hsysDown->Clone("hsysDownf");
//	hsysUpf = (TH1D*)hsysUp->Clone("hsysUpf");
//     }   

   if( nSYS > 0 )
     {	
	TCanvas *c1 = new TCanvas();
	gStyle->SetOptStat(0);
	c1->Clear();
	c1->SetLogy(0);
	hsysDownf->SetFillStyle(0);
	hsysUpf->SetFillStyle(0);
	hsysDownf->SetLineColor(kBlue);
	hsysDownf->SetLineStyle(9);
	hsysDownf->SetLineWidth(2);
	hsysUpf->SetLineColor(kRed);
	hsysUpf->SetLineStyle(1);
	hsysUpf->SetLineWidth(2);
	hsysDownf->Draw("hist");
	hsysUpf->Draw("hist same");
	double xmin = hsysUpf->GetXaxis()->GetBinUpEdge(*nBINS);
	double xmax = hsysUpf->GetXaxis()->GetBinLowEdge(1);
	TLine *l1 = new TLine(xmin,1.0,xmax,1.0);
	l1->SetLineStyle(2);
	l1->SetLineWidth(1);
	l1->Draw("same");

	hsysDownf->GetXaxis()->SetTitle(runName.c_str());
	hsysDownf->GetYaxis()->SetTitle("#sigma(A)/A");
	
//	gStyle->SetOptTitle(1);
//	gStyle->SetTitleColor(253);
//	gStyle->SetTitleFillColor(253);
//	gStyle->SetTitleX(0.5);
//	gStyle->SetTitleW(0.6);
//	gStyle->SetTitleAlign(23);
	
//	hsysDownf->SetTitle(nameSYS[isys].c_str());      
	
	TLegend *legf = new TLegend(0.73,0.85,0.88,0.68);
	legf->SetFillColor(253);
	legf->SetBorderSize(0);
	legf->SetHeader(nameSYS[isys].c_str());
	legf->AddEntry(hsysDownf,"- #sigma","l");
	legf->AddEntry(hsysUpf,"+ #sigma","l");
	legf->Draw();
	
	double ymax = 0.;
	double ymin = 10.;
	for(int i=1;i<=*nBINS;i++)
	  {
	     double y1 = hsysDownf->GetBinContent(i);
	     double y2 = hsysUpf->GetBinContent(i);
	     
	     if( y1 < ymin ) ymin = y1;	
	     if( y2 < ymin ) ymin = y2;
	     
	     if( y1 > ymax ) ymax = y1;
	     if( y2 > ymax ) ymax = y2;
	  }   
	
	hsysDownf->GetYaxis()->SetRangeUser(1.-2.0*(1.-ymin),2.0*(ymax-1.)+1.);
	
	if( producePlots )
	  {
	     std::string picname = "pics/sys"+postfix+".pdf(";
	     if( option == "tag" ) picname = "pics/sys"+postfix+"_tag.pdf(";
	     if( option != "tag" )
	       {	
		  if( isys > 0 && isys != nSYS-1 ) picname = "pics/sys"+postfix+".pdf";
		  else if( isys == nSYS-1 ) picname = "pics/sys"+postfix+".pdf)";
	       }
	     else
	       {
		  if( isys > 0 && isys != nSYS-1 ) picname = "pics/sys"+postfix+"_tag.pdf";
		  else if( isys == nSYS-1 ) picname = "pics/sys"+postfix+"_tag.pdf)";
	       }
	     c1->Print(picname.c_str());
	  }   
	gStyle->SetOptTitle(0);
	c1->Clear();
	delete hsysDownf;
	delete hsysUpf;
	delete hsysDownf_cov;
	delete hsysUpf_cov;
	delete l1;
	delete legf;
	delete c1;
     }   
}

// general correlation factor for each systematics
double CFIT::cfit::bbCorr(TH1D *hsysDown,TH1D *hsysUp)
{
   std::vector<int> bb;
   bb.clear();
   
   for(int ib=1;ib<=*nBINS;ib++)
     {
	double v1 = hsysDown->GetBinContent(ib);
	double v2 = hsysUp->GetBinContent(ib);
	
	if( v1 != 0 || v2 != 0 ) bb.push_back(ib);
     }   
   
   int sz = bb.size();
         
   double rho = 0.;

   if( sz < 2 ) return 1.;
   
   for(int ib1=1;ib1<=sz-1;ib1++)
     {
	for(int ib2=ib1+1;ib2<=sz;ib2++)
	  {	     
	     int id1 = bb[ib1-1];
	     int id2 = bb[ib2-1];
	     
	     double si_Down = hsysDown->GetBinContent(id1);
	     double sj_Down = hsysDown->GetBinContent(id2);
	     
	     double si_Up = hsysUp->GetBinContent(id1);
	     double sj_Up = hsysUp->GetBinContent(id2);
	     
	     double rhoij = (si_Down*sj_Down + si_Up*sj_Up);
	     double norm1 = sqrt(si_Down*si_Down+si_Up*si_Up);
	     double norm2 = sqrt(sj_Down*sj_Down+sj_Up*sj_Up);
	     if( norm1 > 0 && norm2 > 0 ) rhoij /= (norm1*norm2);
	     else rhoij = 0.;	     
	     
	     rho += fabs(rhoij);
	  }
     }   
   
   rho /= (sz*(sz-1)/2);
   
   return rho;
}

// apply fit parameters to the templates
void CFIT::cfit::applySF(std::string option)
{
   TCanvas *c1 = new TCanvas();
   c1->Clear();
   
   THStack *h_draw_st_fit = new THStack();

   TLegend *legf = new TLegend(0.60,0.90,0.88,0.60);
   legf->SetFillColor(253);
   legf->SetBorderSize(0);   

   bool isTag = 0;
   if( strcmp(option.c_str(),"tag") == 0 ) isTag = 1;
         
   TH1D *h_clone[*nT];
   int hasLine[*nT];
   for(int it=0;it<*nT;it++) hasLine[it] = 0;
   
   for(int i=*nT-1;i>=0;i--)
     {
	std::string hname = "hist"+nameT_CURRENT[i]+"_clone";
	h_clone[i] = (TH1D*)(histNOM[i])->Clone(hname.c_str());
	h_clone[i]->Scale(sf[i]);
	PAR[i] = sf[i];
	ERR[i] = sferr[i];
	h_clone[i]->SetMarkerSize(0);

	int colour = 2+i;
	if( isTag ) colour = colourTAG[i];
	else colour = colourT[i];
	
	if( ISMERGEDFIT[i] && isTag ) colour = colourGLUETAG[ISMERGEDFIT[i]-1];
	else if( ISMERGEDFIT[i] && !isTag ) colour = colourGLUE[ISMERGEDFIT[i]-1];
	
	h_clone[i]->SetLineColor(colour);
	h_clone[i]->SetFillColor(colour);
	
	h_clone[i]->Scale(_NDATA);
	h_draw_st_fit->Add(h_clone[i]);
     }   

   h_data->Scale(_NDATA);
   
   h_data->SetMarkerStyle(20);
   h_data->SetMarkerSize(0.8);
//   h_draw_st_fit->Draw("hist e1");
   h_draw_st_fit->Draw("hist");
   h_data->Draw("pe1x0 same");
//   h_data->Draw("pe1 same");
   
   legf->SetHeader(legendName.c_str());
   
   legf->AddEntry(h_data.get(),"Data","ep");

   int mergeCur = 0;
   for(int i=0;i<*nT;i++)
     {
	if( hasLine[i] )
	  {
	     h_clone[i]->SetLineColor(1);
	     h_clone[i]->SetLineWidth(2);	     
	  }	
	else
	  {
	     h_clone[i]->SetLineWidth(0);
	  }	
	
	std::string legName = titleT[i];
	int idx = ISMERGEDFIT[i];
	if( idx > 0 )
	  {
	     if( isTag ) legName = labelGLUETAG[idx-1];
	     else legName = labelGLUE[idx-1];
	     if( mergeCur != idx ) legf->AddEntry(h_clone[i],legName.c_str(),"f");
	     mergeCur = idx;
	  }	
	else legf->AddEntry(h_clone[i],legName.c_str(),"f");
     }   

   h_comb->Scale(_NDATA);
   h_sys_up_combSum->Scale(_NDATA);
   h_sys_down_combSum->Scale(_NDATA);
   
   TGraphAsymmErrors *gr_mc_merged = makeErrorBand(h_comb,
						   h_sys_up_combSum,
						   h_sys_down_combSum);

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   
   h_comb->SetLineColor(kBlack);
   h_comb->Draw("hist same");

//   gr_mc_merged->SetFillStyle(3004);
//   gStyle->SetHatchesLineWidth(2);
//   gStyle->SetHatchesSpacing(0.5);
   //gr_mc_merged->SetFillStyle(3017);
   gr_mc_merged->SetFillStyle(3005);
   gr_mc_merged->SetFillColor(kGray+3);
   gr_mc_merged->Draw("2SAME");

   h_draw_st_fit->GetXaxis()->SetTitle(runName.c_str());
   h_draw_st_fit->GetYaxis()->SetTitle("");
   
   double max = h_draw_st_fit->GetMaximum();
   h_draw_st_fit->SetMaximum(max*1.1);
   
   legf->Draw();

   float maxMC = h_draw_st_fit->GetMaximum();
   float maxData = h_data->GetMaximum();
   h_draw_st_fit->SetMaximum(1.3*std::max(maxMC,maxData));
   
   if( producePlots )
     {	
	std::string fsave = "pics/result.pdf";
	if( option == "tag" ) fsave = "pics/result_tag.pdf";
	c1->Print(fsave.c_str());
	c1->Clear();
	
	// template shape comparison
	float max = 0.;
	
	std::vector<int> isMerged;
	
	TH1D* hCopy[*nT];
	
	for(int i=0;i<*nT;i++)
	  {
	     float integ = h_clone[i]->Integral();
	     
	     hCopy[i] = (TH1D*)h_clone[i]->Clone("hCopy");
	     
	     bool foundMerged = 0;
	     for(unsigned int im=0;im<isMerged.size();im++)
	       {
		  if( ISMERGEDFIT[i] == isMerged[im] )
		    {
		       foundMerged = 1;
		    }		  
	       }	     
	     if( ISMERGEDFIT[i] && !foundMerged ) isMerged.push_back(ISMERGEDFIT[i]);
	     if( foundMerged ) continue;
	     
	     if( ISMERGEDFIT[i] )
	       {
		  for(int it=0;it<*nT;it++)
		    {
		       if( it == i ) continue;
		       if( ISMERGEDFIT[it] == ISMERGEDFIT[i] )
			 {
			    hCopy[i]->Add(h_clone[it]);
			    integ += h_clone[it]->Integral();
			 }		       
		    }		  
	       }	     

	     if( integ > 0. ) hCopy[i]->Scale(1./integ);
	     if( i == 0 ) hCopy[i]->Draw("hist e1");
	     else hCopy[i]->Draw("hist e1 same");
	     hCopy[i]->SetFillColor(0);
	     float hmax = hCopy[i]->GetMaximum();
	     if( hmax > max ) max = hmax;
	     
	     int colour = 2+i;
	     if( isTag ) colour = colourTAG[i];
	     else colour = colourT[i];
	     
	     if( ISMERGEDFIT[i] && isTag ) colour = colourGLUETAG[ISMERGEDFIT[i]-1];
	     else if( ISMERGEDFIT[i] && !isTag ) colour = colourGLUE[ISMERGEDFIT[i]-1];
	     
	     hCopy[i]->SetLineWidth(2);
	     hCopy[i]->SetLineColor(colour);
	  }   	

	TLegend *legShape = new TLegend(0.60,0.90,0.88,0.60);
	legShape->SetFillColor(253);
	legShape->SetBorderSize(0);
	legShape->SetHeader(legendName.c_str());
	
	int mergeCur = 0;	
	for(int i=0;i<*nT;i++)
	  {
	     if( i == 0 ) 
	       {
		  hCopy[i]->SetMaximum(1.1*max);
		  hCopy[i]->GetXaxis()->SetTitle(runName.c_str());
	       }
	     
	     std::string legName = titleT[i];
	     int idx = ISMERGEDFIT[i];
	     
	     if( idx > 0 )
	       {
		  if( isTag ) legName = labelGLUETAG[idx-1];
		  else legName = labelGLUE[idx-1];
		  if( mergeCur != idx ) legShape->AddEntry(hCopy[i],legName.c_str(),"f");
		  mergeCur = idx;
	       }	     
	     else legShape->AddEntry(hCopy[i],legName.c_str(),"f");
	  }   
	legShape->Draw();
	
	std::string fsaveShape = "pics/shape.pdf";
	if( option == "tag" ) fsaveShape = "pics/shape_tag.pdf";
	c1->Print(fsaveShape.c_str());
	c1->Clear();

	delete legShape;
	
	for(int i=0;i<*nT;i++)
	  delete hCopy[i];
     }   
   
   delete legf;
   delete h_draw_st_fit;
   delete gr_mc_merged;
   delete c1;
}

// draw pre and post fit distributions
void CFIT::cfit::drawPrePostFit(std::string option,bool postfit)
{
   TCanvas *c1 = new TCanvas();
   c1->Clear();
   
   THStack *h_draw_st_fit = new THStack();

   TLegend *legf = new TLegend(0.60,0.90,0.88,0.60);
   legf->SetFillColor(253);
   legf->SetBorderSize(0);   

   bool isTag = 0;
   if( strcmp(option.c_str(),"tag") == 0 ) isTag = 1;
         
   TH1D *h_clone[*nT];
   int hasLine[*nT];
   for(int it=0;it<*nT;it++) hasLine[it] = 0;

   float totSum = 0.;
   for(int i=0;i<*nT;i++)
     {
	totSum += _NMC[i];
     }   
   
   for(int i=*nT-1;i>=0;i--)
     {
	std::string hname = "hist"+nameT_CURRENT[i]+"_clone";
	h_clone[i] = (TH1D*)(histNOMTRUE[i])->Clone(hname.c_str());
	h_clone[i]->Scale(totSum);
	if( postfit ) h_clone[i]->Scale(PAR[i]);
	h_clone[i]->SetMarkerSize(0);

	int colour = 2+i;
	if( isTag ) colour = colourTAG[i];
	else colour = colourT[i];
	
	if( ISMERGEDFIT[i] && isTag ) colour = colourGLUETAG[ISMERGEDFIT[i]-1];
	else if( ISMERGEDFIT[i] && !isTag ) colour = colourGLUE[ISMERGEDFIT[i]-1];
	
	h_clone[i]->SetLineColor(colour);
	h_clone[i]->SetFillColor(colour);
	
	h_draw_st_fit->Add(h_clone[i]);
     }   

   h_data->Scale(_NDATA);
   
   h_data->SetMarkerStyle(20);
   h_data->SetMarkerSize(0.8);
   h_draw_st_fit->Draw("hist");
   h_data->Draw("pe1x0 same");
   
   legf->SetHeader(legendName.c_str());
   
   legf->AddEntry(h_data.get(),"Data","ep");

   int mergeCur = 0;
   for(int i=0;i<*nT;i++)
     {
	if( hasLine[i] )
	  {
	     h_clone[i]->SetLineColor(1);
	     h_clone[i]->SetLineWidth(2);	     
	  }	
	else
	  {
	     h_clone[i]->SetLineWidth(0);
	  }	
	
	std::string legName = titleT[i];
	int idx = ISMERGEDFIT[i];
	if( idx > 0 )
	  {
	     if( isTag ) legName = labelGLUETAG[idx-1];
	     else legName = labelGLUE[idx-1];
	     if( mergeCur != idx ) legf->AddEntry(h_clone[i],legName.c_str(),"f");
	     mergeCur = idx;
	  }	
	else legf->AddEntry(h_clone[i],legName.c_str(),"f");
     }   

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   h_draw_st_fit->GetXaxis()->SetTitle(runName.c_str());
   h_draw_st_fit->GetYaxis()->SetTitle("");
   
   double max = h_draw_st_fit->GetMaximum();
   h_draw_st_fit->SetMaximum(max*1.1);
   
   legf->Draw();

   float maxMC = h_draw_st_fit->GetMaximum();
   float maxData = h_data->GetMaximum();
   h_draw_st_fit->SetMaximum(1.3*std::max(maxMC,maxData));
   
   if( producePlots )
     {	
	std::string fsave = "pics/prefitStack.pdf";
	if( postfit ) fsave = "pics/postfitStack.pdf";
	if( option == "tag" ) fsave = "pics/prefitStack_tag.pdf";
	if( option == "tag" && postfit ) fsave = "pics/postfitStack_tag.pdf";
	c1->Print(fsave.c_str());
	c1->Clear();
	
	// template shape comparison
	float max = 0.;
	
	std::vector<int> isMerged;
	
	TH1D* hCopy[*nT];
	
	for(int i=0;i<*nT;i++)
	  {
	     float integ = h_clone[i]->Integral();
	     
	     hCopy[i] = (TH1D*)h_clone[i]->Clone("hCopy");

	     bool foundMerged = 0;
	     for(unsigned int im=0;im<isMerged.size();im++)
	       {
		  if( ISMERGEDFIT[i] == isMerged[im] )
		    {
		       foundMerged = 1;
		    }		  
	       }	     
	     if( ISMERGEDFIT[i] && !foundMerged ) isMerged.push_back(ISMERGEDFIT[i]);
	     if( foundMerged ) continue;
	     
	     if( ISMERGEDFIT[i] )
	       {
		  for(int it=0;it<*nT;it++)
		    {
		       if( it == i ) continue;
		       if( ISMERGEDFIT[it] == ISMERGEDFIT[i] )
			 {
			    hCopy[i]->Add(h_clone[it]);
			    integ += h_clone[it]->Integral();
			 }		       
		    }		  
	       }	     

//	     if( integ > 0. ) hCopy[i]->Scale(1./integ);
	     if( i == 0 ) hCopy[i]->Draw("hist e1");
	     else hCopy[i]->Draw("hist e1 same");
	     hCopy[i]->SetFillColor(0);
	     float hmax = hCopy[i]->GetMaximum();
	     if( hmax > max ) max = hmax;
	     
	     int colour = 2+i;
	     if( isTag ) colour = colourTAG[i];
	     else colour = colourT[i];
	     
	     if( ISMERGEDFIT[i] && isTag ) colour = colourGLUETAG[ISMERGEDFIT[i]-1];
	     else if( ISMERGEDFIT[i] && !isTag ) colour = colourGLUE[ISMERGEDFIT[i]-1];
	     
	     hCopy[i]->SetLineWidth(2);
	     hCopy[i]->SetLineColor(colour);
	  }   	

	TLegend *legShape = new TLegend(0.60,0.90,0.88,0.60);
	legShape->SetFillColor(253);
	legShape->SetBorderSize(0);
	legShape->SetHeader(legendName.c_str());
	
	int mergeCur = 0;	
	for(int i=0;i<*nT;i++)
	  {
	     if( i == 0 ) 
	       {
		  hCopy[i]->SetMaximum(1.1*max);
		  hCopy[i]->GetXaxis()->SetTitle(runName.c_str());
	       }
	     
	     std::string legName = titleT[i];
	     int idx = ISMERGEDFIT[i];
	     
	     if( idx > 0 )
	       {
		  if( isTag ) legName = labelGLUETAG[idx-1];
		  else legName = labelGLUE[idx-1];
		  if( mergeCur != idx ) legShape->AddEntry(hCopy[i],legName.c_str(),"f");
		  mergeCur = idx;
	       }	     
	     else legShape->AddEntry(hCopy[i],legName.c_str(),"f");
	  }   
	legShape->Draw();
	
	std::string fsaveShape = "pics/prefit.pdf";
	if( postfit ) fsaveShape = "pics/postfit.pdf";
	if( option == "tag" ) fsaveShape = "pics/prefit_tag.pdf";
	if( option == "tag" && postfit ) fsaveShape = "pics/postfit_tag.pdf";
	c1->Print(fsaveShape.c_str());
	c1->Clear();

	delete legShape;
	
	for(int i=0;i<*nT;i++)
	  delete hCopy[i];
     }   
   
   delete legf;
   delete h_draw_st_fit;
   delete c1;
}

// uncertainty band creation
TGraphAsymmErrors* CFIT::cfit::makeErrorBand(TH1D* tot, TH1D* plus, TH1D* minus)
{   
   int nbins = tot->GetNbinsX();
   
   const int s = NBINMAX;
   
   double x[s];
   double xerr[s];
   double y[s];
   double yp[s];
   double ym[s];
   
   for(int bin=1;bin<=nbins;++bin)
     {
	int index = bin-1;
	xerr[index] = tot->GetBinWidth(bin)/2.0;
	x[index] = tot->GetBinLowEdge(bin) + xerr[index];
	
	y[index] = tot->GetBinContent(bin);
	yp[index] = plus->GetBinContent(bin)-y[index];
	ym[index] = y[index]-minus->GetBinContent(bin);
	if(y[index] - ym[index] < 0) ym[index] = y[index];
     }
   
   TGraphAsymmErrors* error = new TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp);
   
   return error;
}

void CFIT::cfit::removeMCSys(TH1D *hnom,TH1D *hsysUp,TH1D *hsysDown)
{
   int nb = hnom->GetXaxis()->GetNbins();
   for(int ib=1;ib<=nb;ib++)
     {
        double vb0 = hnom->GetBinContent(ib);
        double ve0 = hnom->GetBinError(ib);
        double vbUp = hsysUp->GetBinContent(ib);

        if( vbUp-vb0 > vb0 )
          {
             hsysUp->SetBinContent(ib,2*vb0);
          }
     }
}

void CFIT::cfit::SetVerbose(int v) { *verb = v; }

double CFIT::cfit::GetChisq() { return *CHISQ; }

int CFIT::cfit::GetNDOF() { return *NDOF; }

double CFIT::cfit::GetPar(int idx) 
{ 
   if( idx >= *nT )
     {
	std::cout << "The number of templates used is " << *nT << std::endl;
	reset();
	exit(1);
     }   
   
   float nMC = 0.;
   for(int i=0;i<*nT;i++)
     {
	nMC += _NMC[i];
     }
   float k = _NDATA/nMC;
   
   if( *VALID ) return PAR[idx]*k;
   else return -1.;
}

double CFIT::cfit::GetParErr(int idx)
{ 
   if( idx >= *nT )
     {
	std::cout << "The number of templates used is " << *nT << std::endl;
	reset();
	exit(1);
     }   

   float nMC = 0.;
   for(int i=0;i<*nT;i++)
     {
	nMC += _NMC[i];
     }
   float k = _NDATA/nMC;
   
   if( *VALID ) return ERR[idx]*k;
   else return -1.;
}

int CFIT::cfit::GetNPar()
{ 
   return *nT;
}

float CFIT::cfit::GetNData()
{ 
   return _NDATA;
}

float CFIT::cfit::GetNTemplate(std::string templateName)
{ 
   int idx = getTemplateId(templateName);
   return _NMC[idx];
}

float CFIT::cfit::GetErrTemplate(std::string templateName)
{ 
   int idx = getTemplateId(templateName);
   return _NMCERR[idx];
}

int CFIT::cfit::getTemplateId(std::string templateName)
{
   int idx = -1;
   
   for(int i=0;i<*nT;i++)
     {
	if( titleT[i] == templateName )
	  {
	     idx = i;
	     break;
	  }	
     }  
   
   return idx;
}

// introduce bin correlation factors and apply general optimisations to systematics
void CFIT::cfit::adjustSys(TH1D *hnom,TH1D *hsys_down,TH1D *hsys_up,int isys)
{
   bool OPT_doMORPH = (optMode == OPT_MORPH || optMode == OPT_MORPH_SGN || 
		       optMode == OPT_MORPH_SIGMA || optMode == OPT_MORPH_SGN_SIGMA);
   bool OPT_doSGN = (optMode == OPT_SGN || optMode == OPT_SGN_SIGMA || 
		     optMode == OPT_MORPH_SGN || optMode == OPT_MORPH_SGN_SIGMA);
   bool OPT_doSIGMA = (optMode == OPT_SIGMA || optMode == OPT_SGN_SIGMA || 
		       optMode == OPT_MORPH_SIGMA || optMode == OPT_MORPH_SGN_SIGMA);
   bool OPT_doSYM = (*covMode != COV_BARLOW);

   // default correlation factors (100%)
   for(int i=0;i<NBINMAX;i++)
     {
	binFactor[isys][i] = 1.;
	
	for(int j=0;j<NBINMAX;j++)
	  {
	     binFactor2D[isys][i][j] = 1.;
	  }	
     }	
   
   Int_t x_nbins = hnom->GetXaxis()->GetNbins();
   
   double val_max = -666.;
   double sigma_av = 0.;
   double sum = 0.;
   
   for(int i=1;i<=x_nbins;i++)
     {
	double v_nom = hnom->GetBinContent(i);
	double v_sys_down = hsys_down->GetBinContent(i);
	double v_sys_up = hsys_up->GetBinContent(i);
	
	double sysUp = (v_nom > 0.) ? (v_sys_up-v_nom)/v_nom : 0.;
	double sysDown = (v_nom > 0.) ? (v_nom-v_sys_down)/v_nom : 0.;
	
	double sign = (sysUp > 0) ? 1. : -1.;
	double sigma = (fabs(sysUp) > fabs(sysDown)) ? fabs(sysUp) : fabs(sysDown);
	sum += v_nom;
	sigma_av += sigma*sigma*v_nom;
	sigma *= sign*v_nom;

	// set systematic variation to the max variation with the sign of sysUp
	if( OPT_doSYM )
	  {	     
	     hsys_up->SetBinContent(i,v_nom+sigma);
	     hsys_down->SetBinContent(i,v_nom-sigma);
	  }	
	
	if( v_nom > val_max ) val_max = v_nom;

	if( OPT_doSGN )
	  {
	     // same-sign up and down variations
	     if( sysUp*sysDown < 0 )
	       {
		  binFactor[isys][i] = 0.;
	       }
	  }
     }
   
   sigma_av = sqrt(sigma_av/sum);
   
   binFactor[isys][x_nbins] = 0.;
   
   for(int i1=1;i1<=x_nbins;i1++)
     {
	double v_nom_i1 = hnom->GetBinContent(i1);
	double v_sys_down_i1 = hsys_down->GetBinContent(i1);
	double v_sys_up_i1 = hsys_up->GetBinContent(i1);
	
	double sysUp_i1 = v_sys_up_i1-v_nom_i1;
	double sysDown_i1 = v_nom_i1-v_sys_down_i1;
	
	double sysUp_rel_i1 = (v_nom_i1 > 0.) ? sysUp_i1/v_nom_i1 : 0.;
	double sysDown_rel_i1 = (v_nom_i1 > 0.) ? sysDown_i1/v_nom_i1 : 0.;
	
	double signUp_i1 = (sysUp_i1 > 0) ? 1. : -1.;
	double signDown_i1 = (sysDown_i1 > 0) ? 1. : -1.;
	
	if( OPT_doMORPH )
	  {
	     if( optMorphMode == OPTMORPH_CUTOFF )
	       {
		  double val_rel = (val_max > 0.) ? v_nom_i1/val_max : 666.;
		  if( val_rel < optMorphFrac ) binFactor[isys][i1] = 0.;
	       }	     
	     else if( optMorphMode == OPTMORPH_GEOMETRIC )
	       {
		  for(int i2=1;i2<=x_nbins;i2++)
		    {		  
		       double v_nom_i2 = hnom->GetBinContent(i2);
		       if( val_max > 0. ) binFactor2D[isys][i1][i2] = sqrt(v_nom_i1*v_nom_i2)/val_max;
		    }		  
	       }	     
	     else if( optMorphMode == OPTMORPH_ARITHMETIC )
	       {
		  for(int i2=1;i2<=x_nbins;i2++)
		    {		  
		       double v_nom_i2 = hnom->GetBinContent(i2);
		       if( val_max > 0. ) binFactor2D[isys][i1][i2] = (v_nom_i1+v_nom_i2)/val_max/2.;
		    }		  
	       }	     
	     else if( optMorphMode == OPTMORPH_FUNC )
	       {
		  for(int i2=1;i2<=x_nbins;i2++)
		    {
		       double v_nom_i2 = hnom->GetBinContent(i2);
		       if( val_max > 0. ) binFactor2D[isys][i1][i2] = v_nom_i1*v_nom_i2/val_max/val_max;
		    }		  
	       }	     
	  }
	
	if( OPT_doSIGMA )
	  {
	     double thres = 0.1*sigma_av;
	     
	     if( fabs(sysUp_rel_i1) < thres )
	       {
		  sysUp_i1 = v_nom_i1*signUp_i1*thres;
		  binFactor[isys][i1] = 0.;
	       }
	     if( fabs(sysDown_rel_i1) < thres )
	       {
		  sysDown_i1 = v_nom_i1*signDown_i1*thres;
		  binFactor[isys][i1] = 0.;
	       }
	     
	     hsys_up->SetBinContent(i1,v_nom_i1+sysUp_i1);
	     hsys_down->SetBinContent(i1,v_nom_i1-sysDown_i1);
	  }
	
	if( i1 < x_nbins )
	  {
	     int i1n = i1+1;
	     
	     double v_nom_i1n = hnom->GetBinContent(i1n);
	     double v_sys_down_i1n = hsys_down->GetBinContent(i1n);
	     double v_sys_up_i1n = hsys_up->GetBinContent(i1n);
	     
	     double sysUp_i1n = v_sys_up_i1n-v_nom_i1n;
	     double sysDown_i1n = v_nom_i1n-v_sys_down_i1n;
	  }
     }
}

void CFIT::cfit::addbin(TH1D *h)
{
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

void CFIT::cfit::adjust(TH1D *h)
{
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   for(int ib=1;ib<=x_nbins;ib++)
     {
	double v = h->GetBinContent(ib);
	if( v == 0 )
	  {
	     h->SetBinContent(ib,1.);
	     h->SetBinError(ib,1.);
	  }	
     }   
}

double* CFIT::cfit::mergeBins(TH1D *h,int ibmin,int ibmax,double xnew[])
{
   bool removeBin = 1;
   
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   
   if( ibmax < x_nbins || ibmin > 1 )
     {
	double sumMax = 0.;
	double sumErrMax = 0.;
	
	double sumMin = 0.;
	double sumErrMin = 0.;
	
	int idx=0;
	for(int ib=1;ib<=ibmax;ib++)
	  {
	     if( ibmin >= ib && ib != 1 ) continue;
	     xnew[idx] = h->GetXaxis()->GetBinLowEdge(ib);
	     idx++;	     
	  }	
	
	for(int ib=1;ib<=ibmin;ib++)
	  {
	     double v = h->GetBinContent(ib);
	     double err = h->GetBinError(ib);
	     sumMin += v;
	     sumErrMin += err*err;
	     
	     if( removeBin )
	       {
		  h->SetBinContent(ib,0.);
		  h->SetBinError(ib,0.);
	       }	     
	  }	

	for(int ib=ibmax;ib<=x_nbins+1;ib++)
	  {
	     double v = h->GetBinContent(ib);
	     double err = h->GetBinError(ib);
	     sumMax += v;
	     sumErrMax += err*err;
	     
	     if( removeBin )
	       {
		  h->SetBinContent(ib,0.);
		  h->SetBinError(ib,0.);
	       }	     
	  }	
	
	xnew[idx] = h->GetXaxis()->GetBinUpEdge(ibmax);
//	xnew[idx] = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetNbins());
	
	sumErrMin = sqrt(sumErrMin);
	sumErrMax = sqrt(sumErrMax);
	
	h->SetBinContent(1,sumMin);
	h->SetBinError(1,sumErrMin);
	
	h->SetBinContent(ibmax,sumMax);
	h->SetBinError(ibmax,sumErrMax);	
     }   
   
   return xnew;
}

void CFIT::cfit::removeMCNeg(TH1D *h1)
{
   int nb = h1->GetXaxis()->GetNbins();
   
   for(int ib=1;ib<=nb;ib++)
     {
	double vb = h1->GetBinContent(ib);
	double ve = h1->GetBinError(ib);
	
	if( vb < 0 || vb <= ve )
	  {
	     h1->SetBinContent(ib,ve);
	     h1->SetBinError(ib,ve);
	  }	
     }   
}

void CFIT::cfit::normalise(TH1D *h1)
{
   double integ = h1->Integral(0,h1->GetXaxis()->GetNbins()+1);
   if( integ > 0. ) h1->Scale(1./integ);
}

void CFIT::cfit::rescale(TH1D *hnom,TH1D *hsysUp,TH1D *hsysDown)
{
   double integ_nom = hnom->Integral(0,hnom->GetXaxis()->GetNbins()+1);
   double integ_up = hsysUp->Integral(0,hsysUp->GetXaxis()->GetNbins()+1);
   double integ_down = hsysDown->Integral(0,hsysDown->GetXaxis()->GetNbins()+1);
   
   if( integ_up > 0. ) hsysUp->Scale(integ_nom/integ_up);
   if( integ_down > 0. ) hsysDown->Scale(integ_nom/integ_down);
}

void CFIT::cfit::SetPlotStyle()
{
   TStyle *styleCurrent = gROOT->GetStyle("PLOT");
   if( !styleCurrent )
     {	
	static TStyle* plotStyle = 0;
	if( plotStyle == 0 ) plotStyle = PlotStyle();
	gROOT->SetStyle("PLOT");
	gROOT->ForceStyle();
     }   
}

TStyle* CFIT::cfit::PlotStyle()
{
   TStyle *plotStyle = new TStyle("PLOT","Plot style");
   
   // use plain black on white colors
   Int_t icol=0; // WHITE
   plotStyle->SetFrameBorderMode(icol);
   plotStyle->SetFrameFillColor(icol);
   plotStyle->SetCanvasBorderMode(icol);
   plotStyle->SetCanvasColor(icol);
   plotStyle->SetPadBorderMode(icol);
   plotStyle->SetPadColor(icol);
   plotStyle->SetStatColor(icol);
   
   // set the paper & margin sizes
   plotStyle->SetPaperSize(20,26);
   
   // set margin sizes
   plotStyle->SetPadTopMargin(0.05);
   plotStyle->SetPadRightMargin(0.05);
   plotStyle->SetPadBottomMargin(0.16);
   plotStyle->SetPadLeftMargin(0.16);
   
   // set title offsets (for axis label)
   plotStyle->SetTitleXOffset(1.4);
   plotStyle->SetTitleYOffset(1.4);
   
   // use large fonts
   Int_t font=42;
   Double_t tsize=0.05;
   plotStyle->SetTextFont(font);
   
   plotStyle->SetTextSize(tsize);
   plotStyle->SetLabelFont(font,"x");
   plotStyle->SetTitleFont(font,"x");
   plotStyle->SetLabelFont(font,"y");
   plotStyle->SetTitleFont(font,"y");
   plotStyle->SetLabelFont(font,"z");
   plotStyle->SetTitleFont(font,"z");
   
   plotStyle->SetLabelSize(tsize,"x");
   plotStyle->SetTitleSize(tsize,"x");
   plotStyle->SetLabelSize(tsize,"y");
   plotStyle->SetTitleSize(tsize,"y");
   plotStyle->SetLabelSize(tsize,"z");
   plotStyle->SetTitleSize(tsize,"z");
   
   // use bold lines and markers
   plotStyle->SetMarkerStyle(20);
   plotStyle->SetMarkerSize(1.2);
   plotStyle->SetHistLineWidth(2.);
   plotStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
   
   // get rid of X error bars
   // plotStyle->SetErrorX(0.001);
   // get rid of error bar caps
   plotStyle->SetEndErrorSize(0.);
   
   // do not display any of the standard histogram decorations
   plotStyle->SetOptTitle(0);
   //plotStyle->SetOptStat(1111);
   plotStyle->SetOptStat(0);
   //plotStyle->SetOptFit(1111);
   plotStyle->SetOptFit(0);
   
   // put tick marks on top and RHS of plots
   plotStyle->SetPadTickX(1);
   plotStyle->SetPadTickY(1);
   
   return plotStyle;
}
