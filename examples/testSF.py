import ROOT as r

r.gSystem.Load("../libCFIT.so")
print "loaded cfit lib"
cf = cfit("test")
print "created cf obj"
