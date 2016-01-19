import os
from colorama import Fore
from CHOMP import PDB, NewSystem, RotamerLibrary, RotamerGenerator, FasterPacker, SimulatedAnnealingPacker, ResidueGenerator
from chomp_setup import dunbrack_lib, energy_function, NewNeighborMap
from chomptools.design import DesignSystem, DesignEnergyGraph 
#from chomptools.cplexpacker import CplexPacker
from time import time
from CHOMP import System
inputpdb = '2zta.pdb'
assert os.path.isfile(inputpdb)
from chomptools.design import Palette
 
### Load a protein normally
P = PDB(inputpdb)
P.RemoveHydrogens()
S = NewSystem(P)
print 'WT E', energy_function(S)
dimer = NewSystem(P)
S = System.Join(dimer)
print 'WT E', energy_function(S)
S.WritePDB('2zta_renum.pdb')


### Create a design palette 
### TODO: make a version of the chomptools.design tools that
### identifies residues by ID and chain, rather than just ID
### this will facilitate the design of oligomers
dp = Palette(S)
dp.AddWT()
## Manually identified interface residues in PyMOL
## Use this line to freeze all non-interfacial resides
#dp.FreezeAll()
#interfacial = set([5, 8, 9, 12, 15, 16, 19, 22, 23, 26, 29, 30, 36, 39, 40, 43, 46, 47, 50, 53, 54, 57, 61])
interfacial = set([5, 8, 9, 12, 15, 36, 39, 40, 43, 46])
for resid in interfacial:
    ## Use this line to just repack interfacial residues
    #dp.AddAA(resid,[S.GetResidue(resid).ResName])
    ## Use this line to do a minimal redesign
    dp.AddAA(resid,'ADFI')

### Here we diverge from a simple repacking calculation.
### DesignSystem converts the provided System into a SystemVector: composite
### composite is essentially a  list of Systems (some of which are 1-residue Systems) 
### composite has an attribute, 'inactive_pairs' which specifies with portions of composite do not coexist
### (i.e. different amino acids for the same residue)
composite = DesignSystem(S, dp)
 
### Proceed with rotamer generation normally
n = NewNeighborMap(composite)
rl = RotamerLibrary(composite)
rg = RotamerGenerator(n, dunbrack_lib)
for r in ResidueGenerator(composite):
    if r.ID in dp.frozen: continue
    rl.GenerateRotamers(rg,r)
 
### Here we diverge from a simple repacking calculation.
### DesignEnergyGraph constructs and fills an EnergyGraph that corresponds to the design problem.
### It uses the 'inactive_pairs' attribute of composite to only calculate meaningful energy terms.
### The resulting EnergyGraph has additional meta-data compared to a standard EnergyGraph.
### Most important, it has a PoseLibrary: pl which takes the place of a RotamerLibrary
### for any downstream combinatorial optimization.
eg = DesignEnergyGraph(composite, rl, n, energy_function)
 
rotdetail = rl.View()

bias = -10
for rotidI,(resI,aaI,chiI) in rotdetail.items():
    for rotidJ,(resJ,aaJ,chiJ) in rotdetail.items():
        if rotidJ <= rotidI: continue
        if resJ != resI + 31: continue  ## only apply sym benefit if residue matches
        if aaJ != aaI: continue  ## only apply sym benefit if amino acid matches
        print rotidI, rotidJ, resI, resJ, aaI, aaJ
        try:
            edgeE = eg.GetEdgeEnergy(rotidI, rotidJ)
            eg.SetEdgeEnergy(rotidI, rotidJ, edgeE + bias)
            print 'Updating edge energy:',edgeE,bias
        except:
            eg.AddEdge(rotidI, rotidJ, bias)
            print 'Adding edge:',bias


### One meta-data item is the combinatorial search size
print 'Combinatorial optimization. Search space = %.3e' % eg.searchsize
before = time()

### For small design problems you can use the FasterPacker
fp = FasterPacker()
#sa = SimulatedAnnealingPacker()
#cpp = CplexPacker()
#combo, E = sa.Pack(eg.pl, eg)  ## Pack with eg.pl (PoseLibrary) not rl (the rotamer library)
combo, E = fp.Pack(eg.pl, eg)  ## Pack with eg.pl (PoseLibrary) not rl (the rotamer library)
print 'Packer E',E
print 'Retest E',eg[combo]

### When you are considering a large design problem with many rotamers the SimulatedAnnealingPacker is faster
#sa = SimulatedAnnealingPacker()
#combo, E = sa.Pack(eg.pl, eg)  ## Pack with eg.pl (PoseLibrary) not rl (the rotamer library)

print time()-before,'seconds for combinatorial optimization'
 
### Here is the final difference for a protein design calculation relative to simple repacking.
#designedS = reconstitute_design(rl, combo)  ## Instead of rl.CreateSystem 
designedS = rl.CreateSystemFromRotamers(combo) ## Instead of rl.CreateSystem 
 
designedP = PDB(designedS)
seqA = designedP['chain A'].seq
seqB = designedP['chain B'].seq

print Fore.GREEN
print 'DESIGN RESULTS'
print 'Chain A:',seqA
print 'Chain B:',seqB
if seqA == seqB:
    print 'Homo-dimer!'
else:
    print 'Hetero-dimer'
print Fore.RESET

### Let's rescore the redesigned output. It will vary a bit from the E that came out of 
### combinatorial optimization. First, there is the usual source of deviation (neglect of
### small interactions in the EnergyGraph. Anoter source of small deviations is that when
### we create a system using CreateSystemFromRotamers coordinate precision is reduced to 
### the precision of a PDB file.
score = energy_function(designedS)
print 'E =', score
open('designed.pdb','w').write(str(designedS))
open('designed.seq','w').write(designedS.GetSequence())

S.Label = 'Wild-type'
designedS.Label = 'Redesigned'
energy_function.CompareSystems(S, designedS)

