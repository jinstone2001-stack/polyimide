from pathlib import Path
import pandas as pd
import json, hashlib, collections

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'output/imide_dielectric'
OUT.mkdir(parents=True, exist_ok=True)
path=ROOT/'huggingface/all_data.csv'
header=pd.read_csv(path,nrows=0).columns.tolist()
counts=collections.Counter(); values={}; samples=[]; rows=0
keys=['class_PIMD','polymer_class','DP','temp','press','check_eq','copoly_type','RadonPy_ver','preset_eq_ver','efdp_permittivity_real','efdp_permittivity_imaginary','efdp_dielectric_loss_tan']
for chunk in pd.read_csv(path,dtype=str,keep_default_na=False,chunksize=10000):
    rows+=len(chunk)
    for c in header: counts[c]+=int((chunk[c]!='').sum())
    for c in keys:
        values.setdefault(c,collections.Counter()).update(chunk[c])
    if len(samples)<6: samples.extend(chunk[['UUID','smiles_list','smiles_1','smiles_search','class_PIMD','DP']].head(6).to_dict('records'))
    sel=chunk[chunk['class_PIMD'].isin(['1','1.0','True'])]
    if len(sel): sel.to_csv(OUT/'notes'/'source_pimd_candidates.csv',mode='a',header=not (OUT/'notes'/'source_pimd_candidates.csv').exists(),index=False)
summary={'rows':rows,'columns':len(header),'nonmissing':dict(counts),'value_counts':{k:dict(v.most_common(30)) for k,v in values.items()},'examples':samples[:6]}
(OUT/'notes'/'source_profile.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k!='nonmissing'},indent=2))
print('Related coverage:', {k:v for k,v in counts.items() if any(x in k for x in ['dielectric','efdp','polarizability','dipole','refractive','volume'])})
