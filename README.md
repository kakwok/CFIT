CFIT is a tool to perform a template fit with including
systematic uncertainties in a correlation matrix. Systematic
correlations between templates are taken into account. The tool also
allows one to estimate statistical and systematic uncertainties
on the fit parameters and to perform the measurement of the scale factors.

Documentation is [available](./doc/doc.pdf)

To use CFIT in CMSSW, 
```
cd CMSSW_8_X_X/src
mkdir RecoBTag
git clone https://github.com/kakwok/CFIT.git
cd CFIT
git checkout pyCFIT
scram b -j8
```
then `CFIT` can be used with python.

The `test/CA15_SF.py` script takes a BTV NTuple to measure the data/MC S.F. with LT method. \n
To run the script, simply do:
`python CA15_SF.py -b`

## Input Constants
```
#Input histogram file
CFITinput = "Final_histograms_sysMerged__rebinned.root"

#Strings to build the histogram names
pTbin = "pt350to2000"
WP    = "DoubleBH"   
LegendHeader = "pT=[350,2000] DoubleB=0.9"

#Options for different run-modes, for a start, turn-off all options
doSYS                = True
doStat               = True   # Need to set doSYS true
doIndividualTemplate = True
do5template          = True
doFittingOpts        = True
```

## Configuration of templates:
The script takes a list of dictionaries to configure the templates.
For each template,
```
label  : label to be used on drawing plots
suffix : suffix used to built histogram names
color  : color of the template
scale  : To be used in systematic unc. calculattion, it scales the histogram normalization before gluing
hisotNames: List of histoNames (including systematics of the template) to be scaled together
```

## SF calculation

MC templates are weighted by fit parameters to match to data.
The post-fit MC events is obtained by 
`result = cf.GetNTemplate(tmp['label']) * cf.GetPar(i)`
Then SF is calculated based on the ratio of pre-fit eff, over post-fit eff
```
effMC       = result['preFitN_Tag']  / result['preFitN_noTag']
effData     = result['postFitN_Tag'] / result['postFitN_noTag']
sf          = effData/effMC
```
