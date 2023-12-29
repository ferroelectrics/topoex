'''
Example of how free energy functional and its first variation should be implemented
to run simulation and estimate energy of the system
'''


def FreeEnergy(alphas, gs, p, phi):
    a1   = alphas['1']
    px,py,pz = [pi for pi in p]
    return a1 * (px**2 + py**2 + pz**2)


def FreeEnergyF(alphas, gs, p, v, phi):
    a1   = alphas['1']
    px,py,pz = [pi for pi in p]
    vx,vy,vz = [vi for vi in v]
    return 2 * a1*(px*vx + py*vy + pz*vz)




