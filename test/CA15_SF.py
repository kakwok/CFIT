import os
import numpy as np
from ROOT import *
from tabulate import tabulate
import argparse
#print tabulate(chi2Table,"firstrow")
def errAB(product,a,errA,b,errB):
    ferr    = ((errA/a)**2+(errB/b)**2)**(0.5)
    return product * ferr
def GetSFtableRow( SFname,  SFdict, normSFdict ):
    row=[SFname]
    row.append('%s'        %(SFdict['chi2_noTag']))
    row.append('%s'        %(SFdict['chi2_Tag']))
    row.append('%.2f +- %.3f'%(SFdict['effMC'],SFdict['effMC_err']))
    row.append('%.2f +- %.3f'%(SFdict['effData'],SFdict['effData_err']))
    row.append('%.5f +- %.5f'%(SFdict['sf'],SFdict['sf_err']))
    if(SFname is not "Nominal"):
        deltaSF = SFdict['sf'] - normSFdict['sf']
        row.append('%.4f'%(deltaSF))
    else:
        row.append(" - ")
    return row

    

def FitAndGetEff(cf,b_template,all_templates,nStat,SysName,printOut=None):
    if nStat is not -1 :        cf.SetStatVariation(nStat)
    #if SysName in all_templates.keys:
    #    
    #elif SysName not "NA" :
    #    cf.SetSysVariation(SysName)
    if SysName is not "NA":
        cf.SetSysVariation(SysName)
        
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
        erra=tmp_result['preFitN_err_noTag'] = 0
        #erra=tmp_result['preFitN_err_noTag'] = cf.GetErrTemplate(tmp['label'])
        tmp_result['postFitN_noTag']         = cf.GetNTemplate(tmp['label']) * cf.GetPar(i)
        tmp_result['postFitN_err_noTag']     = errAB( a*b, a, erra, b,errb)
        mc_sum += tmp_result['postFitN_noTag']
        template_results.append(tmp_result)
    fitTable_row =['Inclusive',chi2pDOF_str,ndata,mc_sum] 
    fitTable.append(fitTable_row)

    print "Running fit on tagged histograms"
    if nStat is not -1 :        cf.SetStatVariation(nStat)
    if SysName is not "NA" :        cf.SetSysVariation(SysName)
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
                erra=result['preFitN_err_Tag']  = 0
                #erra=result['preFitN_err_Tag']  = cf.GetErrTemplate(tmp['label'])
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
        resultTable_row.append("%.3f+-%.3f"%(result['par_noTag'],result['parErr_noTag']))
        resultTable_row.append("%.1f"%result['postFitN_noTag'])
        resultTable_row.append("%.1f"%result['preFitN_Tag'])
        resultTable_row.append("%.3f+-%.3f"%(result['par_Tag'],result['parErr_Tag']))
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
                print "effMC  = %.2f +- %.3f"%(effMC,effMC_err)
                print "effData= %.2f +- %.3f"%(effData,effData_err)
                print "sf     = %.2f +- %.3f"%(sf,sf_err)

    if(printOut is not "noPrint"):
        print tabulate(resultTable,"firstrow") 
        print tabulate(fitTable,"firstrow")
    return {"effMC":effMC,"effData":effData,"sf":sf,"effMC_err":effMC_err,"effData_err":effData_err,"sf_err":sf_err,"chi2_noTag":fitTable[1][1],"chi2_Tag":fitTable[2][1]} 

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

def SetDataAndTemplates(cf,histograms,templates,glueTemplates=None,glueTemplatesTag=None):
# histograms['data']            = {'all':h_nameAll,'tag':h_nameTag,'untag':h_nameUntag}
# histograms[template['label']] = {'all':h_nameAll,'tag':h_nameTag,'untag':h_nameUntag}
    dataHists     = histograms['data']
    cf.SetData(     dataHists['all'])
    cf.SetDataTag(  dataHists['tag'])
    cf.SetDataUntag(dataHists['untag'])
    for template in templates:
        templateHists = histograms[template['label']]
        cf.AddTemplate(template['label']     , templateHists['all'  ],  template['color'])
        cf.AddTemplateTag(template['label']  , templateHists['tag'  ],  template['color'])
        cf.AddTemplateUntag(template['label'], templateHists['untag'],  template['color'])
    if glueTemplates is not None:
        for g in glueTemplates:
            cf.GlueTemplates(ListToVector(g['glueList'],"string"),g['label'],g['color'])
    if glueTemplatesTag is not None:
        for g in glueTemplatesTag:
            cf.GlueTemplatesTag(ListToVector(g['glueList'],"string"),g['label'],g['color'])

def ModifyHistograms(cf,CFITinput,templates,histograms):
    tFile =  TFile.Open(CFITinput)
    CFITupdatedROOTname = CFITinput.replace(".root","_updated.root")
    if os.path.exists(CFITupdatedROOTname):
        print "removing previous scaled file"
        os.system('rm %s'%CFITupdatedROOTname)
    CFITupdatedROOT     = TFile(CFITinput.replace(".root","_updated.root"),"RECREATE")
    for template in templates:
        for hName in template['histoNames']:
            #print "modifying ",hName
            h = tFile.Get(hName)
            hClone = h.Clone()
            hClone.Scale(template['scale'])
            hClone.Write()
    for key in histograms['data'].keys():
        h = tFile.Get(histograms['data'][key])
        hClone = h.Clone()
        hClone.Write()
    tFile.Close()
    CFITupdatedROOT.Close()
    cf.SetInputFile(CFITupdatedROOTname)
    return CFITupdatedROOTname

    
def getHistograms(cf,CFITinput,pTbin,WP,templates,systematics,printTable=None):
    print "reading histograms from ",CFITinput
    tFile = TFile.Open(CFITinput)
    tagList = ["all_"+pTbin,WP+"pass_"+pTbin,WP+"fail_"+pTbin]
    histoNames =[]
    histograms = {}
    histograms['data'] = {}
    for template in templates:
        histograms[template['label']] = {}
        #reset template histoNames
        template['histoNames'] =[]
    

    if doSYS:
        dataHeader = "UNWEIGHTED__DATA__FatJet_JP_"
        QCDheader  = "UNWEIGHTED__QCD__FatJet_JP_"
    else:
        dataHeader = "UNWEIGHTED__DATA__FatJet_JP_"
        QCDheader  = "UNWEIGHTED__QCD__FatJet_JP_"
        #dataHeader = "DATA__FatJet_JP_"
        #QCDheader  = "QCD__FatJet_JP_"
    for tag in tagList :
        if "all" in tag:
            #dataHistoName = dataHeader+tag+"_data_opt"
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            histograms['data']['all'] = dataHistoName
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                template['histoNames'].append(QCDhistoname)
                histograms[template['label']]['all'] = QCDhistoname
        if "pass" in tag:
            #dataHistoName = dataHeader+tag+"_data_opt"
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            histograms['data']['tag'] = dataHistoName
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                template['histoNames'].append(QCDhistoname)
                histograms[template['label']]['tag'] = QCDhistoname
        if "fail" in tag:
            #dataHistoName = dataHeader+tag+"_data_opt"
            dataHistoName = dataHeader+tag+"_data"
            histoNames.append(dataHistoName)
            histograms['data']['untag'] = dataHistoName
            for template in templates:
                QCDhistoname = QCDheader+tag+"_"+template['suffix']
                histoNames.append(QCDhistoname)
                template['histoNames'].append(QCDhistoname)
                histograms[template['label']]['untag'] = QCDhistoname
    #Make table
    preFit_yield = [["Name","Prefit yield(all)","Prefit yield(pass)","Prefit yield(fail)"]]
    data_row    =["data_"+pTbin+"_"+WP]
    for tag in tagList:
        #print dataHistoName
        #dataHistoName = dataHeader+tag+"_data_opt"
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
    if printTable=="printTable":
        print "All histo names:"
        for hname in histoNames:
            print hname
        print tabulate(preFit_yield,"firstrow")
    for template in templates:
        hNames_withoutSys = []
        for hName in template['histoNames']:
            hNames_withoutSys.append(hName)
        for sys in systematics:
            for hName in hNames_withoutSys:
                if not (sys['suffix_up'] in hName or sys['suffix_down']in hName):
                    #print "adding histnames for systematic = ",hName+sys['suffix_up']
                    if not hName+sys['suffix_up'] in template['histoNames']:
                        template['histoNames'].append(hName+sys['suffix_up'])
                    if not sys['suffix_down']==sys['suffix_up']:
                        if not hName+sys['suffix_down'] in template['histoNames']:
                            template['histoNames'].append(hName+sys['suffix_down'])
    return histograms 

#parser = argparse.ArgumentParser()
#parser.add_argument("--inputFile", help="path to the histogram file")
#parser.add_argument("--doSYS"    , help="Run systematic variation")
#parser.add_argument("--doStat"   ,  help="Run statistic variation")
#args = parser.parse_args()
doSYS  =True
doStat =False  # Need to set doSYS true

gInterpreter.Declare("#include \"RecoBTag/CFIT/interface/cfit.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoBTagCFIT.so"))
cf = CFIT.cfit("JP discriminator")
cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
cf.SetMorphing(OPTMORPH_GEOMETRIC)
cf.ProducePlots(1)

#CFITinput = "/afs/cern.ch/work/k/kakwok/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/CFIT_btagval_histograms_fixQCDnorm.root"
#CFITinput = "/afs/cern.ch/user/k/kakwok/work/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/CFIT_btagval_histograms.root"
#CFITinput = "/afs/cern.ch/user/k/kakwok/work/public/Hbb_ISR/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/CFitInput/CFIT_mcJPcalib.root"
#CFITinput = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/Final_histograms_sysMerged.root"
#CFITinput = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/Mu_350_merged/Final_histograms_btagval.root"
#CFITinput = args.inputFile 
#CFITinput = "Final_histograms_sysMerged.root"
#CFITinput  = "/afs/cern.ch/work/b/bmaier/public/MonoHiggs/BTV/CMSSW_8_0_23/src/RecoBTag/BTagValidation/test/SFtemplates_dataJPcalib/Final_histograms_btagval_optimized_doublemu_BTagMu_QCDMuEnriched_HLTAK8Jet300_mcJPcalib_pt350_DoubleBM2.root"

# Use this file to get nominal values+systematics and statistics, doubleB>0.9
#CFITinput = "Final_histograms_sysMerged_rebinned.root"
CFITinput = "Final_histograms_sysMerged__rebinned.root"
#CFITinput = "Final_histograms_sysMerged__rebinned_updated.root"
# Use this file to get DataJPcalib systematics
#CFITinput = "CFIT_btagval_histograms_fixQCDnorm_DataJPcalib_rebin__rebinned.root"
#CFITinput = "Both_ADDBINNING.root"

sysList = [
    'BFRAG',
    'CD',
    'CFRAG',
    'K0L',
    'NTRACKS',
    'PU'
]
templates =[
{'label':'g #rightarrow b#bar{b}', 'suffix':'bfromg','color':2,'scale':1, 'histoNames':[]},
{'label':'b'                     , 'suffix':'b'     ,'color':3,'scale':1, 'histoNames':[]},
{'label':'g #rightarrow c#bar{c}', 'suffix':'cfromg','color':5,'scale':1, 'histoNames':[]},
{'label':'c'                     , 'suffix':'c'     ,'color':4,'scale':1, 'histoNames':[]},
{'label':'l'                     , 'suffix':'l'     ,'color':6,'scale':1, 'histoNames':[]}
]
#templates =[
#{'label':'g #rightarrow b#bar{b}', 'suffix':'bfromg_opt','color':2,'scale':1, 'histoNames':[]},
#{'label':'b'                     , 'suffix':'b_opt'     ,'color':3,'scale':1, 'histoNames':[]},
#{'label':'g #rightarrow c#bar{c}', 'suffix':'cfromg_opt','color':5,'scale':1, 'histoNames':[]},
#{'label':'c'                     , 'suffix':'c_opt'     ,'color':4,'scale':1, 'histoNames':[]},
#{'label':'l'                     , 'suffix':'l_opt'     ,'color':6,'scale':1, 'histoNames':[]}
#]
glueTemplates=[
{'label':'b+g #rightarrow c#bar{c}'  , 'glueList':['b','g #rightarrow c#bar{c}']  ,'color':221},
{'label':'c+light'                   , 'glueList':['c','l']                       ,'color':93}
]
glueTemplatesTag=[
{'label':'other flavours'  , 'glueList':['b','c','l','g #rightarrow c#bar{c}']  ,'color':42},
]
systematics = [
{'label':'bFrag', 'suffix_up':'_BFRAGUP', 'suffix_down':'_BFRAGDOWN' },
{'label':'PU'   , 'suffix_up':'_PUUP', 'suffix_down':'_PUDOWN' },
{'label':'CD'   , 'suffix_up':'_CD', 'suffix_down':'_CD' },
{'label':'CFRAG', 'suffix_up':'_CFRAG', 'suffix_down':'_CFRAG' },
{'label':'K0L'  , 'suffix_up':'_K0L', 'suffix_down':'_K0L' },
{'label':'Ntrks', 'suffix_up':'_NTRACKS', 'suffix_down':'_NTRACKS' },
]
#glueTemplates=None
#glueTemplatesTag=None
pTbin = "pt350to2000"
#pTbin = "pt200to350"
WP    = "DoubleBM2"   # Benedikt uses M2 for 0.75
#WP    = "DoubleBH"   # Benedikt uses H for 0.9
cf.SetLegendHeader("pT=[350,2000] DoubleB=0.75")
#cf.SetLegendHeader("pT=[350,2000] DoubleBL=0.9")
print "Using input file: ",CFITinput

sfTable =[]
sfTable_cols=['Name','chi2(inclusive)','chi2(passed)','effMC','effData','sf','delta_SF']
sfTable.append(sfTable_cols)

cf.SetInputFile(CFITinput)
if doSYS:
    for sys in systematics:
        cf.AddSys( sys['label'] , sys['suffix_up'], sys['suffix_down'])
histograms =  getHistograms(cf,CFITinput,pTbin,WP,templates,systematics,"printTable")

#Collect the histograms
# Set the data and QCD templates 
CFITinput_updated  = ModifyHistograms (cf,CFITinput, templates,histograms)
histograms         = getHistograms(cf,CFITinput_updated,pTbin,WP,templates,systematics,"printTable")

SetDataAndTemplates(cf, histograms, templates,glueTemplates,glueTemplatesTag)
# Do not glue Templates
#SetDataAndTemplates(cf, pTbin, WP, templates)

cf.SetMatrixOption("WRITE")

#labels of templates counted to be efficient
b_template ="g #rightarrow b#bar{b}"
# Fit and get efficiencies
normSFdict  = FitAndGetEff(cf, b_template,templates,-1,"NA")
sfTable.append( GetSFtableRow("Nominal",normSFdict,normSFdict) )

#for template in templates:
#    print template['label'],template['label']=='b'
#    if not (template['label']=='b' or  template['label']=='g #rightarrow c#bar{c}') :
#        continue
#    scaleUp   = 1.5
#    scaleDown = 0.5
#    for scale in [scaleDown,scaleUp]:
#        print "scaling template = %s, by %s"%(scale,template['label'])
#        #Manipulate histograms
#        template['scale'] = scale
#        CFITinput_updated  = ModifyHistograms (cf,CFITinput, templates,histograms)
#        histograms         = getHistograms(cf,CFITinput_updated,pTbin,WP,templates,systematics,"printTable")
#        SetDataAndTemplates(cf, histograms, templates,glueTemplates,glueTemplatesTag)
#
#        SFdict  = FitAndGetEff(cf, b_template,templates,-1,"NA")
#        if scale == scaleUp:
#            tableHeader = template['suffix']+"_up"
#        if scale == scaleDown:
#            tableHeader = template['suffix']+"_down"
#        sfTable.append( GetSFtableRow(tableHeader,SFdict,normSFdict) )
#    #reset scales of all templates to 1
#    for template in templates:
#        template['scale'] = 1
        

# perform statistical variation
if doStat:
    print '-----------------------------'
    print 'Perform statistical variation'
    print '-----------------------------'
    cf.SetMatrixOption("READ")
    avgSF = 0
    nStat = 100
    SFarray =[]
    #The random seed starts from 667
    for i in range(667,667+nStat):
        if i==667+nStat:
            cf.ProducePlots(1)   
        sf_dict = FitAndGetEff(cf, b_template,templates,i,"NA","noPrint")
        SFarray.append(float(sf_dict['sf']))
    npSF  = np.array(SFarray)
    avgSF = np.mean(npSF)
    sdSF = np.std(npSF)
    
    print "Average SF = %.3f , S.D. = %.3f"% (avgSF,sdSF)
    sfTable.append(['Statistics',' - ', '-', '-','-','%.3f'%avgSF, '%.3f'%(sdSF)])
   
#do the calculation of SF here again ....
#perform systematic variation
if doSYS:
    print '-----------------------------'
    print 'Perform systematic variation'
    print '-----------------------------'
    cf.SetMatrixOption("READ")
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_BFRAGUP")
    sfTable.append( GetSFtableRow("BFragUP",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_BFRAGDOWN")
    sfTable.append( GetSFtableRow("BFragDOWN",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_PUUP")
    sfTable.append( GetSFtableRow("PU up",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_PUDOWN")
    sfTable.append( GetSFtableRow("PU down",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_CD")
    sfTable.append( GetSFtableRow("CD",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_K0L")
    sfTable.append( GetSFtableRow("K0L",SFdict,normSFdict) )
    SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_NTRACKS")
    sfTable.append( GetSFtableRow("NTrks",SFdict,normSFdict) )
    #SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_JESup")
    #sfTable.append( GetSFtableRow("JESup",SFdict,normSFdict) )
    #SFdict  = FitAndGetEff(cf, b_template,templates,-1,"_JESdown")
    #sfTable.append( GetSFtableRow("JESdown",SFdict,normSFdict) )

print tabulate(sfTable,"firstrow")

if doSYS:
    for row in sfTable:
        sumw2 = 0.0
        if not sfTable.index(row)==0 and not row[0]=='Nominal':
            sumw2 += float(row[-1])**2
    print "total sys err = ",sumw2**0.5
del cf
gApplication.Terminate()
