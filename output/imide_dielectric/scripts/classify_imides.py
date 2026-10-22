"""RDKit motif screening, keeping a full row-level audit and untouched source strings."""
from pathlib import Path
import sys, json, collections
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'.deps/imide'))
import pandas as pd
from rdkit import Chem, rdBase, RDLogger
RDLogger.DisableLog('rdApp.*')
OUT=ROOT/'output/imide_dielectric'
SMARTS='[#6X3](=[#8X1])[#7X3][#6X3](=[#8X1])'
PATTERN=Chem.MolFromSmarts(SMARTS)
THIO=Chem.MolFromSmarts('[#6X3](=[O,S])[#7X3][#6X3](=[O,S])')
OFFICIAL=Chem.MolFromSmarts('N(-C(=O))-C=O')

def hits(mol, query=PATTERN):
    return mol.GetSubstructMatches(query)

def conventional(mol, h):
    # Other substituent on both carbonyl carbons is C/dummy or H, not N/O/S.
    for ci,oi in [(h[0],h[1]),(h[3],h[4])]:
        if any(a.GetAtomicNum() not in (1,6) for a in mol.GetAtomWithIdx(ci).GetNeighbors() if a.GetIdx() not in (oi,h[2])):
            return False
    return True

def trimer(mol):
    ds=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==0]
    if len(ds)!=2 or any(mol.GetAtomWithIdx(i).GetDegree()!=1 for i in ds): return None
    bonds=[mol.GetAtomWithIdx(i).GetBonds()[0].GetBondType() for i in ds]
    if bonds[0]!=bonds[1] or bonds[0]!=Chem.BondType.SINGLE: return None
    n=mol.GetNumAtoms(); h,t=ds
    hn=mol.GetAtomWithIdx(h).GetNeighbors()[0].GetIdx()
    tn=mol.GetAtomWithIdx(t).GetNeighbors()[0].GetIdx()
    rw=Chem.RWMol(Chem.CombineMols(Chem.CombineMols(mol,mol),mol))
    rw.AddBond(tn,n+hn,bonds[0]); rw.AddBond(n+tn,2*n+hn,bonds[0])
    for i in sorted([t,n+h,n+t,2*n+h],reverse=True): rw.RemoveAtom(i)
    out=rw.GetMol(); Chem.SanitizeMol(out)
    return out

@lru_cache(maxsize=180000)
def check(smi):
    mol=Chem.MolFromSmiles(smi)
    if mol is None: return {'parse':'invalid','imide':False}
    direct=hits(mol); expanded=None; eh=()
    expansion='not_needed'
    # Expand all N-bearing repeat units where a new boundary motif is plausible.
    if any(a.GetAtomicNum()==7 for a in mol.GetAtoms()) and any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):
        try:
            expanded=trimer(mol)
            expansion='linear_trimer' if expanded is not None else 'unsupported_linkers'
            if expanded is not None: eh=hits(expanded)
        except Exception: expansion='sanitization_failed'
    use=mol if direct else expanded
    selected=direct if direct else eh
    cyclic=False; strict=False
    if selected:
        rings=use.GetRingInfo().AtomRings()
        cyclic=any(any({h[0],h[2],h[3]}.issubset(r) and len(r)<=8 for r in rings) for h in selected)
        strict=any(conventional(expanded,h) for h in eh) if expanded is not None else any(conventional(use,h) for h in selected)
    thio=bool(hits(mol,THIO)) or (expanded is not None and bool(hits(expanded,THIO)))
    return {'parse':'ok','imide':bool(selected),'direct':bool(direct),'boundary_only':bool(eh) and not bool(direct),
            'direct_n_count':len({h[2] for h in direct}),'small_cyclic':cyclic,'conventional':strict,
            'os_motif':thio,'official_fr_imide':bool(mol.HasSubstructMatch(OFFICIAL)),
            'canonical':Chem.MolToSmiles(mol),'expansion':expansion,
            'direct_atoms':[list(h) for h in direct]}

def test_rules():
    cases={'succinimide':('O=C1CCC(=O)N1',True), 'phthalimide':('O=C1NC(=O)c2ccccc12',True),
      'acyclic_imide':('CC(=O)NC(C)=O',True),'amide':('CC(=O)NCC',False), 'urea':('NC(=O)N',False),
      'anhydride':('CC(=O)OC(C)=O',False),'imidazole':('c1ncc[nH]1',False),
      'thio_only':('CC(=S)NC(C)=S',False), 'pendant':('*CC(*)N1C(=O)CCC1=O',True),
      'repeat_boundary':('*NC(=O)CC(=O)*',True), 'polyamide_boundary':('*NCCC(=O)*',False),
      'acylurea_broad':('CC(=O)NC(=O)N',True)}
    for name,(s,expected) in cases.items(): assert check(s)['imide']==expected,(name,check(s))
    assert check('*NC(=O)CC(=O)*')['boundary_only']
    assert not check('CC(=O)NC(=O)N')['conventional']
    return len(cases)

def run():
    import faulthandler
    faulthandler.dump_traceback_later(45, repeat=True)
    print('Starting classifier',flush=True)
    tested=test_rules(); selected=[]; audit=[]; excluded=[]; offset=0
    for chunk in pd.read_csv(ROOT/'huggingface/all_data.csv',dtype=str,keep_default_na=False,chunksize=5000):
        for row in chunk.to_dict('records'):
            offset+=1
            if offset % 1000==0: print(f'Row {offset}',flush=True)
            dp=float(row['DP']) if row['DP'] else None
            polymer=dp is not None and dp>1 and '*' in row['smiles_list']
            comps=[]
            for i in range(1,5):
                if not row[f'smiles_{i}']: continue
                ratio=row[f'copoly_ratio_{i}']
                if ratio and float(ratio)<=0: continue
                comps.append((i,check(row[f'smiles_{i}'])))
            confirmed=any(r['imide'] for _,r in comps)
            smi_search_hit=False
            # Complex/ladder linking is represented by supplied connected search structure.
            unsupported=any(r.get('expansion') in ('unsupported_linkers','sanitization_failed') for _,r in comps)
            if polymer and not confirmed and unsupported and row['smiles_search'] and ('N' in row['smiles_search'] or 'n' in row['smiles_search']):
                # Supplied cyclic search graphs can be prohibitively expensive to kekulize.
                # Their explicit bonds/aromatic flags suffice for this local graph query.
                ms=Chem.MolFromSmiles(row['smiles_search'],sanitize=False)
                if ms is not None:
                    ms.UpdatePropertyCache(strict=False)
                    smi_search_hit=bool(hits(ms))
            included=polymer and (confirmed or smi_search_hit)
            matches=[(i,r) for i,r in comps if r['imide']]
            meta={'source_row_number':offset,'is_polymer_record':polymer,'rdkit_imide':confirmed or smi_search_hit,
              'imide_detection_basis':'repeat_unit' if any(r.get('direct') for _,r in matches) else ('linear_trimer_boundary' if matches else ('supplied_smiles_search' if smi_search_hit else 'none')),
              'imide_component_indices':';'.join(str(i) for i,_ in matches),
              'imide_direct_n_count':sum(r.get('direct_n_count',0) for _,r in comps),
              'has_small_cyclic_imide':'' if smi_search_hit else any(r.get('small_cyclic',False) for _,r in matches),
              'has_conventional_carboximide':any(conventional(ms,h) for h in hits(ms)) if smi_search_hit else any(r.get('conventional',False) for _,r in matches),
              'has_os_imide_motif':smi_search_hit or any(r.get('os_motif',False) for _,r in comps),
              'rdkit_official_fr_imide_direct':any(r.get('official_fr_imide',False) for _,r in comps),
              'smiles_parse_status':'ok' if all(r['parse']=='ok' for _,r in comps) else 'invalid_component',
              'class_PIMD_agrees': (row['class_PIMD']=='True') == (confirmed or smi_search_hit) if row['class_PIMD'] else '',
              'canonical_smiles_1':comps[0][1].get('canonical',''),
              'imide_match_atom_indices_json':json.dumps({str(i):r.get('direct_atoms',[]) for i,r in matches},separators=(',',':')),
              'imide_rule_version':'oxygen_imide_v1','rdkit_screen_version':rdBase.rdkitVersion}
            audit.append({'UUID':row['UUID'],**meta,'source_class_PIMD':row['class_PIMD'],'DP':row['DP']})
            if included: selected.append({**meta,**row})
            elif row['class_PIMD']=='True' or confirmed or smi_search_hit: excluded.append({**meta,**row})
        print(f'Processed {offset}; selected {len(selected)}',flush=True)
    pd.DataFrame(selected).to_pickle(OUT/'notes/selected_source.pkl')
    pd.DataFrame(audit).to_pickle(OUT/'notes/classification_audit.pkl')
    pd.DataFrame(excluded).to_pickle(OUT/'notes/excluded_candidates.pkl')
    stats={'source_rows':offset,'imide_polymer_rows':len(selected),'excluded_candidates':len(excluded),'rule_tests_passed':tested,
      'invalid_component_rows':sum(a['smiles_parse_status']!='ok' for a in audit),'polymer_rows':sum(a['is_polymer_record'] for a in audit),
      'detection_basis':dict(collections.Counter(x['imide_detection_basis'] for x in selected)),
      'original_pimd_included':sum(x['class_PIMD']=='True' for x in selected),
      'dynamic_rows':sum(x['efdp_permittivity_real']!='' for x in selected),
      'conventional_rows':sum(x['has_conventional_carboximide'] for x in selected),
      'cyclic_rows':sum(x['has_small_cyclic_imide'] for x in selected)}
    (OUT/'notes/classification_summary.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
    print(json.dumps(stats,indent=2))

if __name__=='__main__': run()
