import os
from CFIT import *
from ROOT import *
from tabulate import tabulate
#print tabulate(chi2Table,"firstrow")

def FitAndGetEff(cf,b_templates,denominator_templates,all_templates):
    print "Running fit on untagged histograms"
    cf.Run()
    printResults(cf,all_templates)
    chi2 = cf.GetChisq()
    print "chi2=%s "% chi2

    ndata = cf.GetNData()
    nmc_btmps         = nmc_all        = 0
    nmc_btmps_tagged  = nmc_all_tagged = 0
    for tmp in b_templates:
        nmc_btmps += cf.GetNTemplate(tmp)
    for tmp in denominator_templates:
        nmc_all   += cf.GetNTemplate(tmp)
    
    print "Running fit on tagged histograms"
    cf.Run("tag")   
    chi2 = cf.GetChisq()
    print "chi2=%s "% chi2
    printResults(cf,all_templates)
    
    ndata_tag = cf.GetNData()
    for tmp in b_templates:
        nmc_btmps_tagged += cf.GetNTemplate(tmp)
    for tmp in denominator_templates:
        nmc_all_tagged   += cf.GetNTemplate(tmp)
    
    # Fraction of btagged templates in all QCD templates 
    #fr     = nmc_btmps/nmc_all
    fr     = nmc_btmps/nmc_all
    # Fraction of btagged templates in all QCD templates among the tagged jets
    fr_tag = nmc_btmps_tagged/nmc_all_tagged
    
    effMC   = nmc_btmps_tagged/nmc_btmps
    coeff   = fr_tag/fr
    effDATA = coeff* (ndata_tag/ndata)
    print "Post fit"
    print "nMC tagged   = %.3f                    nMC_total   = %.3f                                 effMC=nMC_tagged/nMC_total   = %.4f"           %(nmc_btmps_tagged,nmc_btmps, effMC)
    print "fr_tag(b)    = %.3f                    fr(b)       = %.3f" %(fr_tag,fr)
    print "nData tagged = %.3f                    nData_total = %.3f  coeff=fr_tag(b)/fr(b)= %.3f   effData=coeff*(nData tagged/nData_total) = %.4f"%(ndata_tag,ndata, coeff, effDATA)
    print "sf      = %.4f"%(effDATA/effMC)

def printResults(cfit,templates):
    postFit_table=[]
    postFit_header=["Template","PostFit yield","Fitted Weight"]
    postFit_table.append(postFit_header)
    for i,template in enumerate(templates):
        fit_result = "%.3f +- %.3f"%(cfit.GetPar(i),cfit.GetParErr(i))
        postFit_table.append([template['label'],cf.GetNTemplate(template['label']),fit_result])
    print tabulate(postFit_table,"firstrow")
    return
    #for i in range(0,cfit.GetNPar()):
    #    print "par %i = %.3f +- %.3f"%(i,cfit.GetPar(i), cfit.GetParErr(i))
    #return

def ListToVector(pyList,Type):
    if Type=="string":
        v = vector("string")()
    if Type=="double":
        v = vector("double")()
    for element in pyList:
        v.push_back(element)
    return v

def SetDataAndTemplates(cf,pTbin,WP,templates,glueTemplates=None,glueTemplatesTag=None):
    tFile = TFile.Open(CFITinput)
    tagList = ["all_"+pTbin,WP+"pass_"+pTbin,WP+"fail_"+pTbin]
    histoNames =[]
    #dataHeader = "UNWEIGHTED__DATA__FatJet_JP_"
    #QCDheader  = "UNWEIGHTED__QCD__FatJet_JP_"
    dataHeader = "DATA__FatJet_JP_"
    QCDheader  = "QCD__FatJet_JP_"
    for tag in tagList :
        if "all" in tag:
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            #print "Setting data preTag histo = ", dataHistoName
            cf.SetData(dataHistoName)
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                #print "Setting template preTag histo = ", QCDhistoname
                cf.AddTemplate(template['label'], QCDhistoname,  template['color'])
        if "pass" in tag:
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            cf.SetDataTag(dataHistoName)
            #print "Setting data Tagged histo = ", dataHistoName
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                #print "Setting template Tagged histo = ", QCDhistoname
                cf.AddTemplateTag(template['label'], QCDhistoname,  template['color'])
        if "fail" in tag:
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            cf.SetDataUntag(dataHistoName)
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                cf.AddTemplateUntag(template['label'], QCDhistoname,  template['color'])
    if glueTemplates is not None:
        for g in glueTemplates:
            cf.GlueTemplates(ListToVector(g['glueList'],"string"),g['label'],g['color'])
    if glueTemplatesTag is not None:
        for g in glueTemplatesTag:
            cf.GlueTemplatesTag(ListToVector(g['glueList'],"string"),g['label'],g['color'])

    #Make table
    preFit_yield = [["Name","Prefit yield(untagged)","Prefit yield(pass)","Prefit yield(fail)"]]
    data_row    =["data_"+pTbin+"_"+WP]
    for tag in tagList:
        dataHistoName = dataHeader+tag+"_data"
        nEvent = tFile.Get(dataHistoName).GetEntries()
        data_row.append(nEvent)
    preFit_yield.append(data_row) 
    for template in templates:
        template_row=["QCD_"+pTbin+"_"+WP+"_"+template['suffix']]
        for tag in tagList:
            QCDhistoname = QCDheader+tag+"_"+template['suffix']
            nEvent = tFile.Get(QCDhistoname).GetEntries()
            template_row.append(nEvent)
        preFit_yield.append(template_row) 
    print "All histo names:"
    for hname in histoNames:
        print hname
    print tabulate(preFit_yield,"firstrow")

gInterpreter.Declare("#include \"RecoBTag/CFIT/interface/cfit.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoBTagCFIT.so"))
cf = CFIT.cfit("JP discriminator")
cf.SetLegendHeader("pT=[350,2000] DoubleBL=0.75")
cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
cf.SetMorphing(OPTMORPH_GEOMETRIC)
cf.ProducePlots(1)

#CFITinput = "/afs/cern.ch/work/k/kakwok/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/CFIT_btagval_histograms_fixQCDnorm.root"
#CFITinput = "/afs/cern.ch/user/k/kakwok/work/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/CFIT_btagval_histograms.root"
#CFITinput = "/afs/cern.ch/user/k/kakwok/work/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/CFitInput/CFIT_mcJPcalib.root"
#CFITinput = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/Final_histograms_sysMerged.root"
CFITinput = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/Final_histograms_btagval.root"

cf.SetInputFile(CFITinput)
#cf.AddSys("SYS1","_sys1_down","_sys1_up")
#cf.AddSys("SYS2","_sys2_down","_sys2_up")

cf.SetMatrixOption("WRITE")
templates =[
{'label':'b'     , 'suffix':'b'     ,'color':2},
{'label':'bfromg', 'suffix':'bfromg','color':3},
{'label':'c'     , 'suffix':'c'     ,'color':4},
{'label':'cfromg', 'suffix':'cfromg','color':5},
{'label':'l'     , 'suffix':'l'     ,'color':6}
]
glueTemplates=[
{'label':'b+cfromg'  , 'glueList':['b','cfromg']  ,'color':221},
{'label':'c+light'   , 'glueList':['c','l']       ,'color':93}
]
glueTemplatesTag=[
{'label':'others'  , 'glueList':['b','c','l','cfromg']  ,'color':42},
]
pTbin = "pt350to2000"
#pTbin = "pt200to350"
WP    = "DoubleBM2"   # Benedikt uses M2 for 0.75
print "Using input file: ",CFITinput

# Set the data and QCD templates 
SetDataAndTemplates(cf, pTbin, WP, templates,glueTemplates,glueTemplatesTag)

#labels of templates counted to be efficient
b_templates =['bfromg']
#labels of templates to be counted as total
denominator_templates=['b','bfromg','c','cfromg']
# Fit and get efficiencies
FitAndGetEff(cf, b_templates, denominator_templates,templates)

# perform statistical variation
cf.SetMatrixOption("READ")
cf.SetStatVariation(667)
cf.ProducePlots(1)   
cf.Run()
#printResults(cf,templates)
cf.SetStatVariation(667)
cf.Run("tag")
#printResults(cf,templates)
   
#do the calculation of SF here again ....

#perform systematic variation
#cf.SetMatrixOption("READ")
#cf.SetSysVariation("_sys1_up")
#cf.Run()
#printResults(cf)
#cf.SetSysVariation("_sys1_up")
#cf.Run("tag")
#printResults(cf)
   

chi2 = cf.GetChisq()
print "chi2=%s "% chi2
#printResults(cf)

del cf
gApplication.Terminate()
