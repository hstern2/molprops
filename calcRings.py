#  
# Calculation of molecular descriptors and complexity index  
#  
# Implements RDKit  
#  
# Bryon Drown, May 2015  
# Updated Oct. 9, 2015  
# University of Illinois, Urbana-Champaign  
#  
# https://github.com/HergenrotherLab/ctdTools 

import sys  
from rdkit import Chem  
from rdkit.Chem import Descriptors  
from rdkit.ML.Descriptors import MoleculeDescriptors  
from collections import defaultdict  
from collections import OrderedDict  

def calcRingDescriptors(m):
      
  nBonds = m.GetNumBonds()  
  nAtoms = m.GetNumAtoms()  
  cyclomatic = nBonds - nAtoms + 1  
  if(cyclomatic < 1):
       return 0,0

  ri = m.GetRingInfo()  
  if(ri.NumRings() < 1):
       return 0,0
  # get total ring path and nBondRings  
  totalRing = 0  
  Bonds = []  
  Bridges = []  
  for ring in ri.BondRings():
        
    for id in ring:
          
      if (ri.NumBondRings(id) > 1):
            
        Bridges.append(id)  
      totalRing += 1  
      Bonds.append(id)  

  # remove duplicates, then get length    
  nBondRings = len(OrderedDict.fromkeys(Bonds).keys())  
  nBridgeEdges = len(OrderedDict.fromkeys(Bridges).keys())  

  # get nAtomRings  
  Atoms = []  
  for ring in ri.AtomRings():
        
    for id in ring:
          
      Atoms.append(id)      
  nAtomRings = len(OrderedDict.fromkeys(Atoms).keys())  

  # descriptors  
  ringFusionDensity = 2 * float(nBridgeEdges) / float(nAtomRings)  
  ringComplexityIndex = float(totalRing) / float(nAtomRings)  
  molecularCyclizedDegree = float(nAtomRings) / float(nAtoms)  
  nRingSystems = (nBonds - nBondRings) - (nAtoms - nAtomRings) + 1  
  if(nRingSystems < 1):
         
    ringFusionDegree = 0  
  else:
        
    ringFusionDegree = float(cyclomatic) / float(nRingSystems)  

  # set props  
  #m.SetProp('TotalRing', str(totalRing))  
  #m.SetProp('NumBridges', str(nBridgeEdges))  
  #m.SetProp('nBondRings', str(nBondRings))  
  #m.SetProp('nAtomRings', str(nAtomRings))  
  #m.SetProp('ringFusionDensity', str(ringFusionDensity))
  #m.SetProp('ringFusionDegree', str(ringFusionDegree))  
  #m.SetProp('ringComplexityIndex', str(ringComplexityIndex))
  #m.SetProp('molecularCyclizedDegree', str(molecularCyclizedDegree))  
  #m.SetProp('NumRingSystems', str(nRingSystems))  

  return ringFusionDensity, ringComplexityIndex

if __name__ == '__main__':
      
    file_in  = sys.argv[1]  
    file_out = file_in + ".descr.sdf"  
    ms = [x for x in Chem.SDMolSupplier(file_in) if x is not None]  
    ms_wr = Chem.SDWriter(file_out)  

    nms = ('BalabanJ', 'BertzCT', 'FractionCSP3', 'MolWt', 'RingCount')  
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)  

    for i in range(len(ms)):
          
      descrs = calc.CalcDescriptors(ms[i])  
      calcRingDescriptors(ms[i])  
      for x in range(len(descrs)):
            
        ms[i].SetProp(str(nms[x]), str(descrs[x]))  
      ms_wr.write(ms[i])  
