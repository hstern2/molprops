from sys import stderr
from re import search, sub
from rdkit.Chem import Mol, MolFromSmiles, MolToSmiles, RDKFingerprint, AddHs
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io, base64

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

### plots ###

colors = ['red','green','blue','orange','violet','goldenrod','cyan','magenta',
            'gray','sandybrown',
            'darkred','darkgreen','darkblue','darkorange','darkviolet',
            'darkgoldenrod','darkcyan','darkmagenta','darkgray','saddlebrown','lightcoral']

def markerType(s):
    return "o" if search(r"\d", s) else "v"

class Props:
    def __init__(self):
        any = False
        self.QED = []
        self.NP = []
        self.SA = []
        self.Bottcher = []
        self.RFD = []
        self.RCI = []
        self.FSP3 = []
        self.fp = []
        self.NPR1 = []
        self.NPR2 = []

def avg(a):
    return sum(a) / len(a)

def imgstr():
    img = io.BytesIO()
    plt.savefig(img, format='png')
    plt.close()
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf-8')
    return f'<img src="data:image/png;base64,{img_base64}" alt="Plot">\n'

def MolPlots(smiles):
    g = None
    groups = []
    props = {}
    bad_smiles = ''
    for s in smiles.splitlines():
        s = s.strip()
        if not s: # blank line
            continue
        m = search(r'^\#\s*(\S.*)', s)
        if m:
            g = m.group(1)
            if g not in groups:
                groups.append(g)
                props[g] = Props()
            continue
        m = MolFromSmiles(s)
        if not m:
            bad_smiles += f'Bad smiles: {s}<br>\n'
            continue
        if g is None:
            continue
        p = props[g]
        p.any = True
        p.fp.append(RDKFingerprint(m))
        p.QED.append(qed(m))
        p.NP.append(npscorer.scoreMol(m, npmodel))
        p.SA.append(sascorer.calculateScore(m))
        obmol = OBMol()
        PerceiveStereo(obmol)
        obconv = OBConversion()
        obconv.SetInFormat('smi')
        obconv.ReadString(obmol, s)
        p.Bottcher.append(bottchScore.score(obmol))
        RFD, RCI = calcRingDescriptors(m)
        fsp3 = FractionCSP3(m)
        p.RFD.append(RFD)
        p.RCI.append(RCI)
        p.FSP3.append(fsp3)
        m2 = AddHs(m)
        EmbedMolecule(m2)
        MMFFOptimizeMolecule(m2)
        p.NPR1.append(NPR1(m2))
        p.NPR2.append(NPR2(m2))


    groups = [g for g in groups if props[g].any]

    out = ''

    # QED    
    #plt.clf()
    plt.bar(groups, [avg(props[g].QED) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average QED score')
    plt.tight_layout()
    out += imgstr()

    # NP
    #plt.clf()
    plt.bar(groups, [avg(props[g].NP) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average NP-likeness score')
    plt.tight_layout()
    out += imgstr()

    # SA
    #plt.clf()
    plt.bar(groups, [avg(props[g].SA) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average synthetic accessibility score')
    plt.tight_layout()
    out += imgstr()

    # Bottcher
    #plt.clf()
    plt.bar(groups, [avg(props[g].SA) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average Bottcher score')
    plt.tight_layout()
    out += imgstr()

    # ring fusion density
    #plt.clf()
    plt.bar(groups, [avg(props[g].RFD) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average ring fusion density')
    plt.tight_layout()
    out += imgstr()

    # ring complexity index
    #plt.clf()
    plt.bar(groups, [avg(props[g].RCI) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average ring complexity index')
    plt.tight_layout()
    out += imgstr()

    # FSP3 
    #plt.clf()
    plt.bar(groups, [avg(props[g].FSP3) for g in groups], color=colors)
    #plt.xticks(rotation=90)
    plt.ylabel('average fsp3')
    plt.tight_layout()
    out += imgstr()

    # PMI
    #plt.clf()
    for i,g in enumerate(groups):
        plt.scatter(props[g].NPR1,
                    props[g].NPR2,
                    c=colors[i],
                    marker=markerType(g),
                    label=g)
    plt.xticks([])
    plt.yticks([])
    plt.plot((0,0.5),(1,0.5),'black')
    plt.plot((0.5,1),(0.5,1),'black')
    plt.plot((0,1),(1,1),'black')
    plt.legend()
    #plt.legend(bbox_to_anchor=(0,1))
    plt.tight_layout()
    out += imgstr()

    return out
