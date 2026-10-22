"""Produce the two requested CSV extracts and machine-readable validation/provenance."""
from pathlib import Path
import pandas as pd
import numpy as np
import json, hashlib, re
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/imide_dielectric'
SOURCE=ROOT/'huggingface/all_data.csv'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def numeric(d,c): return pd.to_numeric(d[c],errors='coerce')

def run():
    d=pd.read_pickle(OUT/'notes/selected_source.pkl')
    audit=pd.read_pickle(OUT/'notes/classification_audit.pkl')
    excluded=pd.read_pickle(OUT/'notes/excluded_candidates.pkl')
    source_cols=pd.read_csv(SOURCE,nrows=0).columns.tolist()
    real=numeric(d,'efdp_permittivity_real'); imag=numeric(d,'efdp_permittivity_imaginary'); tan=numeric(d,'efdp_dielectric_loss_tan')
    dynamic=real.notna()&imag.notna()&tan.notna()
    d['has_dynamic_dielectric_data']=dynamic
    d['frequency_GHz_inferred']=np.where(dynamic,10,np.nan)
    d['frequency_verified_per_row']=False
    d['frequency_assignment']=np.where(dynamic,'inferred_paper_Fig3c_1084_rows_and_RadonPy_default','no_dynamic_value')
    d['data_origin']='PolyOmics_computational_not_experimental'
    d['permittivity_unit']='dimensionless_relative'
    d['epsilon_imag_derived_real_times_tan']=real*tan
    d['derived_imaginary_method']=np.where(dynamic,'product_of_reported_aggregates_not_original_run_mean','')
    d['raw_imaginary_definition_inconsistent']=np.where(dynamic,~np.isclose(imag,real*tan,rtol=1e-6,atol=1e-10),False)
    d['raw_imaginary_matches_real_div_tan']=np.where(dynamic,np.isclose(imag,real/tan,rtol=1e-6,atol=1e-10),False)
    d['dynamic_nonpositive_real']=dynamic & (real<=0)
    d['dynamic_negative_loss_tan']=dynamic & (tan<0)
    d['dynamic_basic_value_checks_pass']=dynamic & (real>0) & (tan>=0)
    d['equilibration_flag_pass']=d['check_eq'].eq('True')
    d['epsilon_electronic_proxy_n2']=numeric(d,'refractive_index')**2
    d['epsilon_static_corrected_derived']=numeric(d,'static_dielectric_const')+d['epsilon_electronic_proxy_n2']-1
    d['epsilon_static_correction_residual']=numeric(d,'dielectric_const_dc')-d['epsilon_static_corrected_derived']
    d['static_identity_within_rtol_1e_4']=np.where(numeric(d,'dielectric_const_dc').notna() & d['epsilon_static_corrected_derived'].notna(),np.isclose(numeric(d,'dielectric_const_dc'),d['epsilon_static_corrected_derived'],rtol=1e-4,atol=1e-6),False)
    d['structure_scope']=np.where(d['has_conventional_carboximide'].eq(True),'conventional_carboximide_motif_present','broad_imide_like_only_or_unresolved')
    d['dynamic_qc_note']=np.where(dynamic,'raw_imaginary_inconsistent;frequency_inferred;original_efdp_fit_quality_unavailable','dynamic_data_missing')
    d.loc[d['dynamic_nonpositive_real']|d['dynamic_negative_loss_tan'],'dynamic_qc_note']+=';nonphysical_sign_review'
    d.loc[~d['equilibration_flag_pass'],'dynamic_qc_note']+=';source_check_eq_false'
    # Source names and original decimal strings remain unchanged.
    meta=[c for c in d if c not in source_cols]
    front=['UUID','monomer_ID','smiles_list','smiles_1','smiles_search','class_PIMD','polymer_class']
    key=['efdp_permittivity_real','efdp_permittivity_imaginary','efdp_dielectric_loss_tan','static_dielectric_const','dielectric_const_dc','refractive_index','temp','press','check_eq','density','DP']
    common=front+meta+key
    full=d[common+[c for c in source_cols if c not in common]].copy()
    narrow_extra=['smiles_2','smiles_3','smiles_4','copoly_ratio_1','copoly_ratio_2','copoly_ratio_3','copoly_ratio_4','qm_method','charge','forcefield','RadonPy_ver','preset_eq_ver','RDKit_ver','mol_weight_monomer1','qm_polarizability_monomer1','qm_dipole_monomer1','qm_dipole_x_monomer1','qm_dipole_y_monomer1','qm_dipole_z_monomer1','cell_volume','fractional_free_volume']
    narrow=d.loc[dynamic,common+narrow_extra].copy()
    paths=[OUT/'imide_10GHz_dielectric.csv',OUT/'imide_all_dielectric_related.csv']
    narrow.to_csv(paths[0],index=False,encoding='utf-8-sig',na_rep='')
    full.to_csv(paths[1],index=False,encoding='utf-8-sig',na_rep='')
    # Re-open exact delivered files, then independently scan source for UUID/value reconciliation.
    n=pd.read_csv(paths[0],dtype=str,keep_default_na=False)
    f=pd.read_csv(paths[1],dtype=str,keep_default_na=False)
    lookup=f.set_index('UUID',drop=False); requested=set(f.UUID); matched=0; mismatches=[]
    global_dynamic=[]
    for chunk in pd.read_csv(SOURCE,dtype=str,keep_default_na=False,chunksize=10000):
        relevant=chunk[chunk.UUID.isin(requested)]
        for r in relevant.to_dict('records'):
            other=lookup.loc[r['UUID']]; matched+=1
            bad=[c for c in source_cols if r[c]!=other[c]]
            if bad: mismatches.append({'UUID':r['UUID'],'columns':bad})
        global_dynamic.append(chunk.loc[chunk.efdp_permittivity_real.ne(''),['UUID','efdp_permittivity_real','efdp_permittivity_imaginary','efdp_dielectric_loss_tan']])
    gd=pd.concat(global_dynamic); gr=numeric(gd,'efdp_permittivity_real'); gi=numeric(gd,'efdp_permittivity_imaginary'); gt=numeric(gd,'efdp_dielectric_loss_tan')
    narrow_expected=set(d.loc[dynamic,'UUID'])
    nr=numeric(n,'efdp_permittivity_real'); nt=numeric(n,'efdp_dielectric_loss_tan'); ni=numeric(n,'epsilon_imag_derived_real_times_tan')
    assert not f.UUID.duplicated().any() and not n.UUID.duplicated().any()
    assert set(n.UUID)==narrow_expected and set(n.UUID).issubset(set(f.UUID))
    assert matched==len(f) and not mismatches
    assert np.allclose(ni,nr*nt,rtol=1e-12,atol=1e-12)
    assert f.smiles_parse_status.eq('ok').all()
    checks={'source_rows':len(audit),'source_columns':len(source_cols),'imide_rows':len(f),'dynamic_imide_rows':len(n),
       'full_export_columns':len(f.columns),'dynamic_export_columns':len(n.columns),
       'source_values_exactly_preserved_rows':matched,'source_value_mismatches':mismatches,
       'uuid_duplicates_full':int(f.UUID.duplicated().sum()),'uuid_duplicates_dynamic':int(n.UUID.duplicated().sum()),
       'unique_canonical_smiles1_full':f.canonical_smiles_1.nunique(),'duplicate_canonical_smiles1_rows':int(f.canonical_smiles_1.duplicated().sum()),
       'source_dynamic_rows':len(gd),'global_imaginary_matches_division':int(np.isclose(gi,gr/gt,rtol=1e-6,atol=1e-10).sum()),
       'global_imaginary_matches_product':int(np.isclose(gi,gr*gt,rtol=1e-6,atol=1e-10).sum()),
       'full_static_available':int(f.static_dielectric_const.ne('').sum()),'full_corrected_static_available':int(f.dielectric_const_dc.ne('').sum()),
       'full_polarizability_monomer1_available':int(f.qm_polarizability_monomer1.ne('').sum()),
       'full_check_eq_false':int(f.check_eq.eq('False').sum()),'dynamic_check_eq_false':int(n.check_eq.eq('False').sum()),
       'dynamic_nonpositive_real':int((nr<=0).sum()),'dynamic_negative_tan':int((nt<0).sum()),
       'dynamic_conventional_carboximide':int(n.has_conventional_carboximide.eq('True').sum()),
       'full_conventional_carboximide':int(f.has_conventional_carboximide.eq('True').sum()),
       'full_source_PIMD_true':int(f.class_PIMD.eq('True').sum()),'full_source_PIMD_false':int(f.class_PIMD.eq('False').sum()),
       'full_source_PIMD_blank':int(f.class_PIMD.eq('').sum()),
       'source_pimd_excluded':int(excluded.class_PIMD.eq('True').sum()),
       'excluded_nonpolymer_oxygen_imide':int((excluded.rdkit_imide.eq(True)&excluded.is_polymer_record.eq(False)).sum()),
       'excluded_polymer_thio_only':int((excluded.has_os_imide_motif.eq(True)&excluded.rdkit_imide.eq(False)&excluded.is_polymer_record.eq(True)).sum()),
       'frequency_direct_verification_available':False,'checks_passed':True,
       'dynamic_real_range':[float(nr.min()),float(nr.max())],'dynamic_tan_range':[float(nt.min()),float(nt.max())],
       'detection_basis':f.imide_detection_basis.value_counts().to_dict()}
    (OUT/'validation_summary.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
    audit.to_json(OUT/'notes/classification_audit.jsonl',orient='records',lines=True,force_ascii=False)
    excluded.to_json(OUT/'notes/excluded_candidates.jsonl',orient='records',lines=True,force_ascii=False)
    manifest={'created_local_date':'2026-10-22','timezone':'Asia/Seoul','source_url':'https://huggingface.co/datasets/yhayashi1986/PolyOmics','source_license':'CC-BY-4.0',
      'source_files':[{ 'path':str(p.relative_to(ROOT)), 'sha256':sha(p)} for p in [SOURCE,ROOT/'2511.11626v1.pdf',ROOT/'huggingface/README.md',ROOT/'Radonpy/radonpy/core/poly.py',ROOT/'Radonpy/radonpy/sim/preset/ef_dp.py',ROOT/'Radonpy/AutoMD_scripts/6_ef_dp.py']],
      'outputs':[{ 'path':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in paths],
      'rdkit_version':'2026.09.1','python_version':__import__('sys').version,'script_files':[{ 'path':p.name,'sha256':sha(p)} for p in Path(__file__).parent.glob('*.py')]}
    (OUT/'provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    dictionary={}
    readme=(ROOT/'Radonpy/AutoMD_scripts/README.md').read_text(encoding='utf-8')
    for line in readme.splitlines():
        if line.startswith('|'):
            cells=[x.strip() for x in line.split('|')]
            if len(cells)>3 and cells[1]: dictionary[cells[1]]=cells[2]
    info=[]
    for c in f.columns:
        meaning=dictionary.get(c,'')
        if not meaning:
            base=re.sub(r'monomer[234]','monomer1',c)
            meaning=dictionary.get(base,'')
        if 'polarizability' in c: unit='angstrom^3'
        elif 'dipole' in c: unit='Debye'
        elif c=='temp' or c=='tg': unit='K'
        elif c=='press': unit='atm'
        elif c=='density': unit='g/cm^3'
        elif 'permittivity' in c or 'dielectric' in c or 'refractive' in c or c.startswith('epsilon_'): unit='dimensionless_or_flag_see_description'
        else: unit='see_original_description'
        info.append({'column':c,'origin':'source_unchanged' if c in source_cols else 'derived_or_audit','description':meaning,'unit':unit,'nonempty_imide_rows':int(f[c].ne('').sum()),'nonempty_dynamic_rows':int(n[c].ne('').sum()) if c in n else None})
    (OUT/'column_dictionary.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(checks,ensure_ascii=False,indent=2))

if __name__=='__main__': run()
