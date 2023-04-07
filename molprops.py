from sys import stderr
from re import sub
from rdkit.Chem import Mol, MolFromSmiles, MolToSmiles, RDKFingerprint, MolFromInchi
from rdkit.Chem.QED import qed
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit.Chem.Lipinski import FractionCSP3
import npscorer, sascorer
from openbabel.openbabel import OBConversion, OBMol, PerceiveStereo
from bottchscore.bottchscore3 import BottchScore
from calcRings import calcRingDescriptors

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
