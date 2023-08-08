from sys import stderr
from re import sub
from rdkit.Chem import Mol, MolFromSmiles, MolToSmiles, RDKFingerprint
from rdkit.Chem.QED import qed
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem.Lipinski import FractionCSP3
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from rdkit.Chem.Descriptors3D import NPR1, NPR2
import npscorer, sascorer
from openbabel.openbabel import OBConversion, OBMol, PerceiveStereo
from bottchscore.bottchscore3 import BottchScore
from calcRings import calcRingDescriptors
from sammon.sammon import sammon
import matplotlib.pyplot as plt

npmodel = npscorer.readNPModel()
bottchScore = BottchScore()

def MolProps(smiles):

    tab = '''
<table>
<tr><th></th><th>QED</th><th>NP</th><th>SA</th><th>B&ouml;ttcher</th><th>RFD</th><th>RCI</th><th>fsp3</th><th>SMILES</th></tr>
'''

    bad_smiles = ''

    for s in smiles.splitlines():
        s = sub(r'\s', '', s)
        if not s:
            continue
        m = MolFromSmiles(s)
        if not m:
            bad_smiles += f'Bad smiles: {s}<br>\n'
            continue
        mc = Mol(m.ToBinary())
        rdDepictor.Compute2DCoords(mc)
        drawer = rdMolDraw2D.MolDraw2DSVG(200,200)
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        obmol = OBMol()
        PerceiveStereo(obmol)
        obconv = OBConversion()
        obconv.SetInFormat('smi')
        obconv.ReadString(obmol, s)
        ringFusionDensity, ringComplexityIndex = calcRingDescriptors(m)
        tab +='<tr>\n'
        tab += f'<td>{svg}</td>\n'
        tab += f'<td>{qed(m):.3f}</td>\n'
        tab += f'<td>{npscorer.scoreMol(m,npmodel):.3f}</td>\n'
        tab += f'<td>{sascorer.calculateScore(m):.3f}</td>\n'
        tab += f'<td>{bottchScore.score(obmol):.3f}</td>\n'
        tab += f'<td>{ringFusionDensity:.3f}</td>\n'
        tab += f'<td>{ringComplexityIndex:.3f}</td>\n'
        tab += f'<td>{FractionCSP3(m):.3f}</td>\n'
        tab += f'<td>{MolToSmiles(m)}</td>\n'
        tab += '</tr>\n'

    tab += '</table>\n'
    return '<p></p>' + bad_smiles + '<p></p>' + tab

######## molplots ########

'''

def markerType(s):
    return "o" if search(r"\d", s) else "v"

for s in smiles.splitlines():
    if s.startswith('#'):
        avgQED = 0
        avgNP = 0
        avgSA = 0
        avgBottcher = 0
        avgRFD = 0
        avgRCI = 0
        avgFSP3 = 0
        n = 0
        key = sub(r'\.[^\.]*$','',fname)
        keys.append(key)
        istart[key] = len(fp)
        continue
            ln = ln.rstrip()
            m = MolFromInchi(ln) if ln.startswith('InChI') else MolFromSmiles(ln)
            fp.append(RDKFingerprint(m))
            avgQED += qed(m)
            avgNP += npscorer.scoreMol(m, npmodel)
            avgSA += sascorer.calculateScore(m)
            obmol = OBMol()
            PerceiveStereo(obmol)
            obconv = OBConversion()
            obconv.SetInFormat('inchi' if ln.startswith('InChI') else 'smi')
            obconv.ReadString(obmol, ln)
            avgBottcher += bottchScore.score(obmol)
            RFD, RCI = calcRingDescriptors(m)
            fsp3 = FractionCSP3(m)
            avgRFD += RFD
            avgRCI += RCI
            avgFSP3 += fsp3
            m2 = AddHs(m)
            EmbedMolecule(m2)
            MMFFOptimizeMolecule(m2)
            npr1.append(NPR1(m2))
            npr2.append(NPR2(m2))
            n += 1
    iend[key] = len(fp)
    if n == 0:
        n = 1
    avg[key] = (avgQED/n, avgNP/n, avgSA/n, avgBottcher/n, avgRFD/n, avgRCI/n, avgFSP3/n)

print(f'writing csv file...', file=stderr)
with (open(a.o+'.csv', 'w')) as f:
    print('file,avgQED,avgNP,avgSA,avgBottcher,avgRFD,avgRCI,avgFSP3', file=f)
    for key in keys:
        avgQED, avgNP, avgSA, avgBottcher, avgRFD, avgRCI, avgFSP3 = avg[key]
        print(f'{key},{avgQED:.3f},{avgNP:.3f},{avgSA:.3f},{avgBottcher:.3f},{avgRFD:.3f},{avgRCI:.3f},{avgFSP3:.3f}', file=f)

colors = ['red','green','blue','orange','violet','goldenrod','cyan','magenta',
            'gray','sandybrown',
            'darkred','darkgreen','darkblue','darkorange','darkviolet',
            'darkgoldenrod','darkcyan','darkmagenta','darkgray','saddlebrown','lightcoral']

print('making QED plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][0] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average QED score')
plt.tight_layout()
plt.savefig(a.o+'.qed.png')

print('making NP plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][1] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average NP-likeness score')
plt.tight_layout()
plt.savefig(a.o+'.NP.png')

print('making SA plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][2] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average synthetic accessibility score')
plt.tight_layout()
plt.savefig(a.o+'.SA.png')

print('making Bottcher plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][3] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average Bottcher score')
plt.tight_layout()
plt.savefig(a.o+'.Bottcher.png')

print('making ring fusion density plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][4] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average ring fusion density')
plt.tight_layout()
plt.savefig(a.o+'.RFD.png')

print('making ring complexity index plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][5] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average ring complexity index')
plt.tight_layout()
plt.savefig(a.o+'.RCI.png')

print('making fsp3 plot...', file=stderr)
plt.clf()
plt.bar(keys, [avg[k][6] for k in keys], color=colors)
plt.xticks(rotation=90)
plt.ylabel('average fsp3')
plt.tight_layout()
plt.savefig(a.o+'.fsp3.png')

print('making sammon plot...', file=stderr)
D = numpy.array([[1-FingerprintSimilarity(i,j) for j in fp] for i in fp])
xy,E = sammon(D,2,inputdist='distance',maxiter=10000)
plt.clf()
for i,key in enumerate(keys):
    plt.scatter([p[0] for p in xy[istart[key]:iend[key]]],
                [p[1] for p in xy[istart[key]:iend[key]]],
                c=colors[i],
                marker=markerType(key),
                label=key)
plt.xticks([])
plt.yticks([])
plt.legend(bbox_to_anchor=(0,1))
plt.tight_layout()
plt.savefig(a.o+'.sammon.png')

print('making PMI plot...', file=stderr)
plt.clf()
for i,key in enumerate(keys):
    plt.scatter(npr1[istart[key]:iend[key]],
                npr2[istart[key]:iend[key]],
                c=colors[i],
                marker=markerType(key),
                label=key)
plt.xticks([])
plt.yticks([])
plt.plot((0,0.5),(1,0.5),'black')
plt.plot((0.5,1),(0.5,1),'black')
plt.plot((0,1),(1,1),'black')
plt.legend()
#plt.legend(bbox_to_anchor=(0,1))
plt.tight_layout()
plt.savefig(a.o+'.pmi.png')

'''