from ROOT import *

def ListToVector(pyList,Type):
    if Type=="string":
        v = vector("string")()
    if Type=="double":
        v = vector("double")()
    for element in pyList:
        v.push_back(element)
    return v

l = ['a','n','b']
v = ListToVector(l,"string")
for i in range(0,v.size()):
    print v[i]
