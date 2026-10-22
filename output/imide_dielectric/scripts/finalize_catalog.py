from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'output/imide_dielectric'
descriptions={
 'source_row_number':'헤더 제외 1-based 원본 레코드 순번',
 'is_polymer_record':'DP>1이고 실제 smiles_list에 반복 연결점이 있음',
 'rdkit_imide':'산소 imide 모티프 독립 RDKit 판정',
 'imide_detection_basis':'repeat_unit 직접 검출 / linear_trimer_boundary 반복 경계 검출',
 'imide_component_indices':'모티프가 검출된 1-based 성분 번호; 세미콜론 구분',
 'imide_direct_n_count':'성분별 원래 반복단위에 직접 매칭되는 N 개수의 합; 조성 가중값 아님',
 'has_small_cyclic_imide':'8원 이하 실제 고리에 N과 두 carbonyl C가 함께 있는 imide 모티프',
 'has_conventional_carboximide':'두 carbonyl의 반대쪽 치환기가 C/H인 모티프를 적어도 하나 포함',
 'has_os_imide_motif':'RadonPy처럼 carbonyl O/S를 모두 허용한 넓은 모티프',
 'rdkit_official_fr_imide_direct':'공식 RDKit fr_imide 패턴의 원래 성분 직접 매칭; aromatic/경계 제한 있음',
 'smiles_parse_status':'포함 성분 SMILES의 RDKit sanitize 파싱 상태',
 'class_PIMD_agrees':'원본 class_PIMD와 새 산소 imide 모티프 여부의 일치',
 'canonical_smiles_1':'첫 번째 반복단위의 RDKit canonical SMILES; 다성분 전체 구조 식별자 아님',
 'imide_match_atom_indices_json':'원본 성분의 0-based RDKit 원자 매칭 index. 경계만 검출 시 직접 목록은 빈 배열',
 'imide_rule_version':'프로젝트 구조 분류 규칙 버전',
 'rdkit_screen_version':'실제 분류에 사용한 RDKit 버전',
 'has_dynamic_dielectric_data':'원본 동적 세 열 모두 숫자로 존재',
 'frequency_GHz_inferred':'문헌 Fig3c, 1084행 및 코드 기본값에서 추정한 주파수; 행별 직접 검증 아님',
 'frequency_verified_per_row':'원시 동적 파일에서 각 행의 주파수 직접 검증 여부; 모두 False',
 'frequency_assignment':'주파수 추정 근거 또는 동적 데이터 결측',
 'data_origin':'PolyOmics 계산값; 실험값 아님',
 'permittivity_unit':'유전율이 상대유전율임을 표시; F/m 아님',
 'epsilon_imag_derived_real_times_tan':'보고된 real*tan 집계값으로 유도한 허수부; 원시 반복 평균 복원 아님',
 'derived_imaginary_method':'허수부 유도 방법과 집계 한계',
 'raw_imaginary_definition_inconsistent':'원본 imaginary가 real*tan과 불일치. 동적 값 결측 시 False=미평가',
 'raw_imaginary_matches_real_div_tan':'원본 imaginary가 real/tan과 일치. 동적 값 결측 시 False=미평가',
 'dynamic_nonpositive_real':'동적 값이 존재하며 real<=0',
 'dynamic_negative_loss_tan':'동적 값이 존재하며 tan<0',
 'dynamic_basic_value_checks_pass':'동적 값 존재, real>0, tan>=0. 적합 정확도/주파수 검증은 아님',
 'equilibration_flag_pass':'원본 check_eq가 True; 동적 적합 품질과 구분',
 'epsilon_electronic_proxy_n2':'굴절률 제곱의 전자분극 근사. 전체 10GHz 유전율 아님',
 'epsilon_static_corrected_derived':'static_dielectric_const+refractive_index^2-1로 유도',
 'epsilon_static_correction_residual':'원본 dielectric_const_dc - 유도 보정 정적값',
 'static_identity_within_rtol_1e_4':'정적 보정 관계 rtol=1e-4, atol=1e-6 일치 여부. 입력 결측 False=미평가',
 'structure_scope':'conventional carboximide motif 있음 / broad imide-like만 존재',
 'dynamic_qc_note':'허수부 불일치, 주파수 추정, EFDP 적합 정보 결측 등',
 'efdp_permittivity_real':'원본 보고 동적 상대유전율 실수부. 추가 전자분극 보정 적용하지 않음',
 'efdp_permittivity_imaginary':'원본 보고 허수부 열. 모든 존재값이 real/tan으로 정의에 불일치. 직접 사용 부적절',
 'efdp_dielectric_loss_tan':'원본 보고 손실탄젠트 tan(delta)',
 'static_dielectric_const':'고정 전하 MD 쌍극자 요동 기반 정적 상대유전율. 전자분극 보정 전',
 'dielectric_const_dc':'static_dielectric_const + refractive_index^2 - 1의 전자분극 보정 정적 상대유전율',
 'r2':'평균제곱 말단간 거리 관련 평형 속성; EFDP fit R^2 아님',
}
dictionary=json.loads((OUT/'column_dictionary.json').read_text(encoding='utf-8'))
booleans={'is_polymer_record','rdkit_imide','has_small_cyclic_imide','has_conventional_carboximide','has_os_imide_motif','rdkit_official_fr_imide_direct','class_PIMD_agrees','has_dynamic_dielectric_data','frequency_verified_per_row','raw_imaginary_definition_inconsistent','raw_imaginary_matches_real_div_tan','dynamic_nonpositive_real','dynamic_negative_loss_tan','dynamic_basic_value_checks_pass','equilibration_flag_pass','static_identity_within_rtol_1e_4'}
for item in dictionary:
    c=item['column']
    if c in descriptions: item['description']=descriptions[c]
    if not item['description']: item['description']='원본 PolyOmics 열 그대로 보존. 원본 README 및 논문 Table S3 참조.'
    if c in booleans: item['unit']='boolean'
    elif c=='frequency_GHz_inferred': item['unit']='GHz'
    elif c.startswith('epsilon_') or c in ['static_dielectric_const','dielectric_const_dc','efdp_permittivity_real','efdp_permittivity_imaginary','efdp_dielectric_loss_tan'] or c.startswith('refractive_index'): item['unit']='dimensionless'
    elif item['origin']=='derived_or_audit': item['unit']='metadata_or_count'
    if 'polarizability_sos_' in c and any(f'_{w}_' in c for w in [486,589,656]): item['description']+=' 광학 파장 486/589/656 nm 표기이며 GHz 아님.'
(OUT/'column_dictionary.json').write_text(json.dumps(dictionary,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
manifest=json.loads((OUT/'provenance.json').read_text(encoding='utf-8'))
for item in manifest['outputs']: assert sha(OUT/item['path'])==item['sha256']
manifest['script_files']=[{'path':str(p.relative_to(OUT)),'sha256':sha(p)} for p in (OUT/'scripts').glob('*.py')]
manifest['documentation']=[{'path':p,'sha256':sha(OUT/p)} for p in ['README_ko.md','validation_summary.json','independent_review.json','column_dictionary.json']]
(OUT/'provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('Column definitions finalized; both CSV hashes unchanged; provenance finalized.')
