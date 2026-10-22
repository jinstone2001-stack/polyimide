"""Independent inspection of source applicability and chemical edge cases."""
import sys, json, collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'.deps/imide'))
import pandas as pd
from rdkit import Chem,RDLogger
RDLogger.DisableLog('rdApp.*')
cols=['UUID','smiles_list','smiles_search','DP','class_PIMD','polymer_class']
cols += [f'{p}_{i}' for p in ('smiles','copoly_ratio') for i in range(1,5)]
df=pd.read_csv(ROOT/'huggingface/all_data.csv',usecols=cols,dtype=str,keep_default_na=False)
poly=pd.to_numeric(df.DP,errors='coerce').gt(1)&df.smiles_list.str.contains('*',regex=False)
stats={'total':len(df),'polymer_rows':int(poly.sum()),'multiple_component_polymer_rows':int((poly&df.smiles_2.ne('')).sum()),
 'dp_gt1_without_dummy':int((pd.to_numeric(df.DP,errors='coerce').gt(1)&~df.smiles_list.str.contains('*',regex=False)).sum()),
 'dp_missing_with_dummy':int((df.DP.eq('')&df.smiles_list.str.contains('*',regex=False)).sum())}
for i in range(1,5):
    stats[f'component{i}_present_missing_ratio']=int((poly&df[f'smiles_{i}'].ne('')&df[f'copoly_ratio_{i}'].eq('')).sum())
all_smis=set()
for i in range(1,5): all_smis.update(df.loc[poly,f'smiles_{i}'].unique())
all_smis.discard('')
core=Chem.MolFromSmarts('[#6X3](=[#8X1])[#7X3][#6X3](=[#8X1])')
relaxed=Chem.MolFromSmarts('[#6](=[#8])[#7][#6](=[#8])')
edge=[]; linker_counts=collections.Counter(); invalid=[]
for smi in all_smis:
    mol=Chem.MolFromSmiles(smi)
    if mol is None:
        invalid.append(smi);continue
    if any(a.GetAtomicNum()==7 for a in mol.GetAtoms()):
        ds=[a for a in mol.GetAtoms() if a.GetAtomicNum()==0]
        linker_counts[len(ds)]+=1
        if mol.HasSubstructMatch(relaxed) and not mol.HasSubstructMatch(core):
            edge.append({'smiles':smi,'canonical':Chem.MolToSmiles(mol)})
stats.update(unique_polymer_components=len(all_smis),n_bearing_dummy_counts=dict(linker_counts),
 invalid_unique_component_smiles=invalid,relaxed_only_motifs=edge)
(ROOT/'output/imide_dielectric/notes/source_edgecase_review.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
print(json.dumps(stats,indent=2))
