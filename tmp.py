from CHOMP import PDB, NewSystem, RotamerLibrary, RotamerGenerator, FasterPacker, SimulatedAnnealingPacker, System, ResidueGenerator
from chomp_setup import dunbrack_lib, energy_function, NewNeighborMap
from chomptools.design import DesignSystem, DesignEnergyGraph, Palette 
from time import time
import os

#################################################################################
###############################   LOAD INPUT FILE (DERIVED FROM 1zaa.pdb) #######
#################################################################################
### Note: I've renumbered the provided pdb to match the numbering referred to in the paper
inputpdb = '2l7iCH1.pdb'
assert os.path.isfile(inputpdb)  ### What does this line do?

f=open('log.txt','w')

### Load a protein normally
P = PDB(inputpdb)
P.RemoveHydrogens()
S = System.Join(NewSystem(P))

from CHOMP import AminoAcidIndex, ResidueBurialState, ExtraChiType

def invariant_huge_rotamer_generator(neighbormap, dunbrack_lib):
  rotgen = RotamerGenerator(neighbormap, dunbrack_lib)
  rgp = rotgen.GetParameters()
  rgp.baseBuriedCutoff = 0
  rgp.buriedCutoff = 0
  rgp.maxPercentSurface = rgp.maxPercentBuried
  rgp.maxRotamersSurface = rgp.maxRotamersBuried
  rotgen.SetParameters(rgp)
  for aa in AminoAcidIndex.values.values():
    for burial in ResidueBurialState.values.values():
      for chi in [0,1,2]:
        # enable 1xSD chis 0,1 for ALL amino acids and ALL burial states
        rotgen.SetExtraChiFlag(aa, burial, chi, ExtraChiType.EXTRA_CHI_1xSD)
        rotgen.SetExtraChiFlag(aa, burial, chi, ExtraChiType.EXTRA_CHI_2xSD)
  return rotgen

def invariant_large_rotamer_generator(neighbormap, dunbrack_lib):
  rotgen = RotamerGenerator(neighbormap, dunbrack_lib)
  rgp = rotgen.GetParameters()
  rgp.baseBuriedCutoff = 0
  rgp.buriedCutoff = 0
  rgp.maxPercentSurface = rgp.maxPercentBuried
  rgp.maxRotamersSurface = rgp.maxRotamersBuried
  rotgen.SetParameters(rgp)
  for aa in AminoAcidIndex.values.values():
    for burial in ResidueBurialState.values.values():
      for chi in [0,1]:
        # enable 1xSD chis 0,1 for ALL amino acids and ALL burial states
        rotgen.SetExtraChiFlag(aa, burial, chi, ExtraChiType.EXTRA_CHI_1xSD)
        rotgen.SetExtraChiFlag(aa, burial, chi, ExtraChiType.EXTRA_CHI_2xSD)
  return rotgen

n=NewNeighborMap(S)

### Calculate the energy function score for the System before making any changes
### Use the ScoreSystem method of energy_function
wtscore = energy_function.ScoreSystem(S, n) 
out= 'WT E '+str(wtscore)
print out
f.write(out+'\n')

#################################################################################
###############################   SETUP DESIGN PALETTE HERE  ####################
#################################################################################

### The goal is to mimic the Dahiyat & Mayo, (1997) paper
### First, create an empty design palette
### Note that the Palette object wants a System for setup
dp = Palette(S) 
### Now add polar amino acids to the surface positions
### Try using the dp.AddAA method
core='AVLIFWM'
#surf='ASTHDNQEKR'

corepos = [280,284,287,291,294,312,315,319,322,326,329]
Bcorepos = [x + 58 for x in corepos]
active = corepos + Bcorepos

print dp.View()

dp.FreezeAll()
for rid in corepos: dp.AddAA(rid, core)
for rid in Bcorepos: dp.AddAA(rid, core)

### Use the Palette SearchSize method to quantify the number of 
### candidate sequences in the combinatorial space
seqspacesize = dp.SearchSize() 
out= str(seqspacesize)+' possible sequences'
print out
f.write(out+'\n')

### Here we diverge from a simple repacking calculation.
### DesignSystem converts the provided System into a SystemVector: composite
### composite is essentially a  list of Systems (some of which are 1-residue Systems) 
### composite has an attribute, 'inactive_pairs' which specifies with portions of composite do not coexist
### (i.e. different amino acids for the same residue)
composite = DesignSystem(S, dp)
 
### Proceed with rotamer generation normally
n = NewNeighborMap(composite)
rl = RotamerLibrary(composite)
#rg = RotamerGenerator(n, dunbrack_lib)
rg = invariant_large_rotamer_generator(n, dunbrack_lib)
for r in ResidueGenerator(composite):
    if r.ID in dp.frozen: continue
    rl.GenerateRotamers(rg, r)

### Note: the RotamerLibrary SearchSize is not accurate
### since the different amino acid alternatives for 
### design positions cannot coexist 
 
### Here we diverge from a simple repacking calculation.
### DesignEnergyGraph constructs and fills an EnergyGraph that corresponds to the design problem.
### It uses the 'inactive_pairs' attribute of composite to only calculate meaningful energy terms.
### The resulting EnergyGraph has additional meta-data compared to a standard EnergyGraph.
### Most important, it has a PoseLibrary: pl which takes the place of a RotamerLibrary
### for any downstream combinatorial optimization.
eg = DesignEnergyGraph(composite, rl, n, energy_function)
 
### One meta-data item is the combinatorial search size
out='Combinatorial optimization. Search space = %.3e' % eg.searchsize
#print out
#f.write(out+'\n')

### Let's time how long the Packer takes for optimization
before = time() 
### For small design problems you can use the FasterPacker
#fp = FasterPacker()
### For really small design problems you might use the CplexPacker
#cpp = CplexPacker()
sa = SimulatedAnnealingPacker()
combo, E = sa.Pack(eg.pl, eg)  ## Pack with eg.pl (PoseLibrary) not rl (the rotamer library)

out=str(time()-before)+' seconds for combinatorial optimization'
print out
f.write(out+'\n') 
### Here is the final difference for a protein design calculation relative to simple repacking.
designedS = rl.CreateSystemFromRotamers(combo) ## Instead of rl.CreateSystem 
 
### Let's rescore the redesigned output. It will vary a bit from the E that came out of 
### combinatorial optimization. First, there is the usual source of deviation (neglect of
### small interactions in the EnergyGraph. Another source of small deviations is that when
### we create a system using CreateSystemFromRotamers coordinate precision is reduced to 
### the precision of a PDB file.

### Use the ScoreSystem function to score the new System: designedS
### You will need to make a NeighborMap for the new System
nn = NewNeighborMap(designedS)
score = energy_function.ScoreSystem(designedS,nn) 
out='E = '+str(score)
print out
f.write(out+'\n')

f.close()

open('designed3.pdb','w').write(str(designedS))
open('designed3.seq','w').write(designedS.GetSequence())

print ' or resi '.join(map(str,active))
