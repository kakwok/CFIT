import os
from ROOT import *
from tabulate import tabulate
#print tabulate(chi2Table,"firstrow")
def errAB(product,a,errA,b,errB):
    ferr    = ((errA/a)**2+(errB/b)**2)**(0.5)
    return product * ferr

def FitAndGetEff(cf,b_template,all_templates,printOut=None):
    print "Running fit on untagged histograms"
    cf.Run()
    #printResults(cf,all_templates)
    chi2pDOF = cf.GetChisq()
    Ndof = cf.GetNDOF()
    chi2pDOF_str = "%.3f/%s = %.3f"%(chi2pDOF*Ndof,Ndof,chi2pDOF)

    #Table for results of yields
    resultTable =[]
    resultTable_cols = ['TemplateName','preFitN(Inclusive)','weight(inclusive)','postFitN(Inclusive)','preFitN(Tagged)','weight(Tagged)','postFitN(Tagged)']
    resultTable.append(resultTable_cols)
    #Table for results of fits
    fitTable    =[]
    fitTable_cols   =['Fit','chi2/dof','data','MC(sum)']
    fitTable.append(fitTable_cols)

    ndata = cf.GetNData()
    mc_sum=0
    #List of dic. to contain the results
    template_results=[]
    #Put fit parameters and 
    for i,tmp in enumerate(all_templates):
        tmp_result ={}
        tmp_result['label']                  = tmp['label']
        b   =tmp_result['par_noTag']         = cf.GetPar(i)
        errb=tmp_result['parErr_noTag']      = cf.GetParErr(i)
        a   =tmp_result['preFitN_noTag']     = cf.GetNTemplate(tmp['label'])
        erra=tmp_result['preFitN_err_noTag'] = cf.GetErrTemplate(tmp['label'])
        tmp_result['postFitN_noTag']         = cf.GetNTemplate(tmp['label']) * cf.GetPar(i)
        tmp_result['postFitN_err_noTag']     = errAB( a*b, a, erra, b,errb)
        mc_sum += tmp_result['postFitN_noTag']
        template_results.append(tmp_result)
    fitTable_row =['Inclusive',chi2pDOF_str,ndata,mc_sum] 
    fitTable.append(fitTable_row)
    print "Running fit on tagged histograms"
    cf.Run("tag")   
    chi2pDOF = cf.GetChisq()
    Ndof = cf.GetNDOF()
    chi2pDOF_str = "%.3f/%s = %.3f"%(chi2pDOF*Ndof,Ndof,chi2pDOF)


    #printResults(cf,all_templates)
    ndata_tag = cf.GetNData()
    mc_sum = 0
    #Save Tagged results
    for i,tmp in enumerate(all_templates):
        for result in template_results:
            if(result['label']==tmp['label']):
                b   =result['par_Tag']          = cf.GetPar(i)
                errb=result['parErr_Tag']       = cf.GetParErr(i)
                a   =result['preFitN_Tag']      = cf.GetNTemplate(tmp['label'])
                erra=result['preFitN_err_Tag']  = cf.GetErrTemplate(tmp['label'])
                result['postFitN_Tag']          = cf.GetNTemplate(tmp['label']) * cf.GetPar(i)
                result['postFitN_err_Tag']      = errAB( a*b,a,erra,b,errb) 
                mc_sum += result['postFitN_Tag']
    fitTable_row =['Tagged',chi2pDOF_str,ndata_tag,mc_sum] 
    fitTable.append(fitTable_row)

    #Calculate eff and SF from results
    for result in template_results:
        resultTable_row = []
        resultTable_row.append(result['label'])
        resultTable_row.append("%.1f"%result['preFitN_noTag'])
        resultTable_row.append("%.3f"%result['par_noTag'])
        resultTable_row.append("%.1f"%result['postFitN_noTag'])
        resultTable_row.append("%.1f"%result['preFitN_Tag'])
        resultTable_row.append("%.3f"%result['par_Tag'])
        resultTable_row.append("%.2f"%result['postFitN_Tag'])
        resultTable.append(resultTable_row)
        if result['label'] == b_template:
            effMC       = result['preFitN_Tag']  / result['preFitN_noTag']
            effData     = result['postFitN_Tag'] / result['postFitN_noTag']
            sf          = effData/effMC
            effMC_err   = errAB(effMC  , result['preFitN_Tag'], result['preFitN_err_Tag'], result['preFitN_noTag'], result['preFitN_err_noTag'])
            effData_err = errAB(effData, result['postFitN_Tag'], result['postFitN_err_Tag'], result['postFitN_noTag'], result['postFitN_err_noTag'])
            sf_err      = errAB(sf     , effMC, effMC_err, effData, effData_err) 
            if(printOut is not "noPrint"):
                print "effMC  = %.3f +- %.3f"%(effMC,effMC_err)
                print "effData= %.3f +- %.3f"%(effData,effData_err)
                print "sf     = %.3f +- %.3f"%(sf,sf_err)

    if(printOut is not "noPrint"):
        print tabulate(resultTable,"firstrow") 
        print tabulate(fitTable,"firstrow")
    return {"effMC":effMC,"errData":effData,"sf":sf,"effMC_err":effMC_err,"effData_err":effData_err,"sf_err":sf_err,"chi2_noTag":fitTable[0][1],"chi2_Tag":fitTable[0][1]} 

def printResults(cfit,templates):
    postFit_table=[]
    postFit_header=["Template","PostFit GetNTemplate","Fit Par.","Post-fit yield"]
    postFit_table.append(postFit_header)
    for i,template in enumerate(templates):
        fit_result = "%.3f +- %.3f"%(cfit.GetPar(i),cfit.GetParErr(i))
        postFit_yield = cf.GetNTemplate(template['label'])*cfit.GetPar(i)
        postFit_table.append([template['label'],cf.GetNTemplate(template['label']),fit_result,postFit_yield])
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
            cf.SetData(dataHistoName)
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                cf.AddTemplate(template['label'], QCDhistoname,  template['color'])
        if "pass" in tag:
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            cf.SetDataTag(dataHistoName)
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
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
        nEvent = tFile.Get(dataHistoName).Integral()
        data_row.append(nEvent)
    preFit_yield.append(data_row) 
    for template in templates:
        template_row=["QCD_"+pTbin+"_"+WP+"_"+template['suffix']]
        for tag in tagList:
            QCDhistoname = QCDheader+tag+"_"+template['suffix']
            nEvent = tFile.Get(QCDhistoname).Integral()
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
#CFITinput = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/Final_histograms_btagval.root"
CFITinput = "/Final_histograms_btagval.root"

sysList = [
    'BFRAG',
    'CD',
    'CFRAG',
    'K0L',
    'NTRACKS',
    'PU'
]
templates =[
{'label':'bfromg', 'suffix':'bfromg','color':2},
{'label':'b'     , 'suffix':'b'     ,'color':3},
{'label':'cfromg', 'suffix':'cfromg','color':5},
{'label':'c'     , 'suffix':'c'     ,'color':4},
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
#WP    = "DoubleBH"   # Benedikt uses H for 0.9
print "Using input file: ",CFITinput


cf.SetInputFile(CFITinput)
cf.AddSys("Bfrag","_BFRAG_down","_BFRAG_up")

cf.SetMatrixOption("WRITE")
# Set the data and QCD templates 
SetDataAndTemplates(cf, pTbin, WP, templates,glueTemplates,glueTemplatesTag)
#SetDataAndTemplates(cf, pTbin, WP, templates)

#labels of templates counted to be efficient
b_templates =['bfromg']
# Fit and get efficiencies
normSFdict=FitAndGetEff(cf, 'bfromg',templates)

# perform statistical variation
print '-----------------------------'
print 'Perform statistical variation'
print '-----------------------------'
cf.SetMatrixOption("READ")
avgSF = 0
nStat = 10
for i in range(667,667+nStat):
    cf.SetStatVariation(678)
    if i==667+nStat:
        cf.ProducePlots(1)   
    sf_dict = FitAndGetEff(cf, 'bfromg',templates,"noPrint")
    avgSF += float(sf_dict['sf'])
print "Average SF = %.3f "% (avgSF/nStat)
#cf.Run()
#cf.SetStatVariation(667)
#cf.Run("tag")
   
#do the calculation of SF here again ....

#perform systematic variation
#cf.SetMatrixOption("READ")
#cf.SetSysVariation("_sys1_up")
#cf.Run()
#printResults(cf)
#cf.SetSysVariation("_sys1_up")
#cf.Run("tag")
#printResults(cf)

del cf
gApplication.Terminate()
