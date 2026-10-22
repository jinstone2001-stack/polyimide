"""Independent graph-neighbor check; does not reuse classifier SMARTS/oligomer builder."""
from pathlib import Path
from functools import lru_cache
import sys,json
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'output/imide_dielectric'
sys.path.insert(0,str(ROOT/'.deps/imide'))
import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

@lru_cache(maxsize=180000)
def independent_match(s):
    m=Chem.MolFromSmiles(s)
    if m is None: return None
    dummy=[a for a in m.GetAtoms() if a.GetAtomicNum()==0]
    pair={}
    if len(dummy)==2 and all(a.GetDegree()==1 and a.GetBonds()[0].GetBondType()==Chem.BondType.SINGLE for a in dummy):
        pair={dummy[0].GetIdx():dummy[1].GetNeighbors()[0],dummy[1].GetIdx():dummy[0].GetNeighbors()[0]}
    for n in m.GetAtoms():
        if n.GetAtomicNum()!=7 or n.GetTotalDegree()!=3: continue
        count=0
        for bond in n.GetBonds():
            if bond.GetBondType() not in (Chem.BondType.SINGLE,Chem.BondType.AROMATIC): continue
            c=bond.GetOtherAtom(n)
            if c.GetAtomicNum()==0 and c.GetIdx() in pair: c=pair[c.GetIdx()]
            if c.GetAtomicNum()!=6 or c.GetTotalDegree()!=3: continue
            if any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(c).GetAtomicNum()==8 and b.GetOtherAtom(c).GetTotalDegree()==1 for b in c.GetBonds()): count+=1
        if count>=2: return True
    return False

def run():
    f=pd.read_csv(OUT/'imide_all_dielectric_related.csv',dtype=str,keep_default_na=False)
    n=pd.read_csv(OUT/'imide_10GHz_dielectric.csv',dtype=str,keep_default_na=False)
    got=set(f.UUID); expected=set(); total=0
    cols=['UUID','DP','smiles_list']+[f'smiles_{i}' for i in range(1,5)]+[f'copoly_ratio_{i}' for i in range(1,5)]
    for ch in pd.read_csv(ROOT/'huggingface/all_data.csv',dtype=str,keep_default_na=False,usecols=cols,chunksize=10000):
        for row in ch.to_dict('records'):
            total+=1
            if not row['DP'] or float(row['DP'])<=1 or '*' not in row['smiles_list']: continue
            if any(independent_match(row[f'smiles_{i}']) for i in range(1,5) if row[f'smiles_{i}'] and (not row[f'copoly_ratio_{i}'] or float(row[f'copoly_ratio_{i}'])>0)): expected.add(row['UUID'])
    missing=sorted(expected-got); extra=sorted(got-expected)
    numeric=lambda c:pd.to_numeric(f[c],errors='coerce')
    reconstructed=numeric('static_dielectric_const')+numeric('refractive_index')**2-1
    static=numeric('dielectric_const_dc'); available=static.notna()&reconstructed.notna()
    residual=static-reconstructed
    totaldip=(numeric('qm_dipole_x_monomer1')**2+numeric('qm_dipole_y_monomer1')**2+numeric('qm_dipole_z_monomer1')**2)**.5
    dip=numeric('qm_dipole_monomer1'); da=dip.notna()&totaldip.notna()
    result={'source_rows_reviewed':total,'independent_graph_membership_rows':len(expected),'membership_matches_exactly':not missing and not extra,'missing_uuids':missing,'extra_uuids':extra,
      'static_identity_available':int(available.sum()),'static_identity_close_rtol_1e_4':int((available & np.isclose(static,reconstructed,rtol=1e-4,atol=1e-6)).sum()),'static_identity_max_abs_residual':float(residual.abs().max()),
      'dipole_identity_available':int(da.sum()),'dipole_identity_close_rtol_1e_4':int((da & np.isclose(dip,totaldip,rtol=1e-4,atol=1e-6)).sum()),'dipole_identity_max_abs_residual':float((dip-totaldip).abs().max()),
      'dynamic_unique_monomer_ids':int(n.monomer_ID.nunique()),'dynamic_unique_repeat_smiles':int(n.canonical_smiles_1.nunique()),
      'dynamic_conventional_unique_repeat_smiles':int(n.loc[n.has_conventional_carboximide.eq('True'),'canonical_smiles_1'].nunique()),
      'full_complete_component_signature_count':int(f[['canonical_smiles_1','smiles_2','smiles_3','smiles_4','copoly_ratio_1','copoly_ratio_2','copoly_ratio_3','copoly_ratio_4']].drop_duplicates().shape[0]),
      'conventional_imide_rows_without_dynamic_data':int((f.has_conventional_carboximide.eq('True') & f.efdp_permittivity_real.eq('')).sum()),
      'dynamic_records_all_equilibrated':bool(n.check_eq.eq('True').all())}
    (OUT/'independent_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    assert not missing and not extra

if __name__=='__main__': run()
