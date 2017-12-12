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
