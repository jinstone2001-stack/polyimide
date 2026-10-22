"""Independent chemical query regression checks; not the production extractor."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / '.deps' / 'imide'))
from rdkit import Chem, rdBase
from rdkit.Chem import Fragments

query = Chem.MolFromSmarts('[#6;X3](=[#8;X1])[#7;X3][#6;X3](=[#8;X1])')
squery = Chem.MolFromSmarts('[#6;X3](=[O,S])[#7;X3][#6;X3](=[O,S])')
cases = [
    ('succinimide', 'O=C1CCC(=O)N1', True),
    ('N-methylphthalimide', 'CN1C(=O)c2ccccc2C1=O', True),
    ('diacetamide', 'CC(=O)NC(C)=O', True),
    ('amide', 'CC(=O)NC', False),
    ('urea', 'NC(=O)N', False),
    ('anhydride', 'CC(=O)OC(C)=O', False),
    ('imidazole', 'c1ncc[nH]1', False),
    ('thioimide', 'CC(=S)NC(C)=S', False),
    ('acylurea', 'CC(=O)NC(=O)N', True),
    ('boundary_raw', '*N(C)C(=O)CC(=O)*', False),
    ('pendant_imide_polymer', '*CC(*)N1C(=O)CCC1=O', True),
    ('aromatic_uracil_representation', 'O=c1cc[nH]c(=O)[nH]1', True),
]
results = []
for name, smi, expected in cases:
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None, name
    actual = mol.HasSubstructMatch(query)
    assert actual == expected, (name, actual, expected)
    results.append(dict(name=name, smiles=smi, canonical=Chem.MolToSmiles(mol),
                        oxygen_imide=actual, sulfur_allowed_imide=mol.HasSubstructMatch(squery),
                        rdkit_builtin_fr_imide=Fragments.fr_imide(mol)))

mol = Chem.MolFromSmiles('*N(C)C(=O)CC(=O)*')
size = mol.GetNumAtoms()
dummy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
head, tail = dummy
hn = mol.GetAtomWithIdx(head).GetNeighbors()[0].GetIdx()
tn = mol.GetAtomWithIdx(tail).GetNeighbors()[0].GetIdx()
joined = Chem.RWMol(Chem.CombineMols(Chem.CombineMols(mol, mol), mol))
for i in (0, 1):
    joined.AddBond(tn + i * size, hn + (i + 1) * size, Chem.BondType.SINGLE)
for idx in sorted((tail, head+size, tail+size, head+2*size), reverse=True):
    joined.RemoveAtom(idx)
trimer = joined.GetMol()
Chem.SanitizeMol(trimer)
assert trimer.HasSubstructMatch(query)
results.append(dict(name='boundary_linear_trimer', canonical=Chem.MolToSmiles(trimer),
                    oxygen_imide=True, unique_imide_centres=len({m[2] for m in trimer.GetSubstructMatches(query)})))
payload = dict(rdkit_version=rdBase.rdkitVersion, passed=len(results), results=results)
Path(__file__).with_suffix('.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
print(json.dumps(payload, indent=2))
