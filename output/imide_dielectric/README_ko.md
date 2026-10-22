# PolyOmics imide 계열 유전율 데이터 추출 및 검토

작성일: 2026-10-22 (Asia/Seoul). 원본 `huggingface/all_data.csv` 176,096행, 263열을 분석했다. 입력 CSV, 논문, RadonPy 소스는 수정하지 않았다. 결과는 **계산 데이터**이며 실험값이 아니다.

## 결과 파일과 범위

| 파일 | 레코드 수 | 열 수 | 내용 |
|---|---:|---:|---|
| [imide_10GHz_dielectric.csv](imide_10GHz_dielectric.csv) | 59 | 75 | imide 모티프를 가진 고분자 중 동적 유전율 3개 원본 값이 있는 계산 레코드. 10 GHz는 아래 문헌 근거에 따른 추정이다. |
| [imide_all_dielectric_related.csv](imide_all_dielectric_related.csv) | 14,909 | 299 | 모든 추출 imide 고분자의 정적·동적 유전율, 쌍극자, 분극률 텐서, 굴절률, 광학 분산, 밀도, 분자량, 자유부피 및 계산 조건. 누락 방지와 추적을 위해 원본 263개 열 전체와 추가 분류·검증 36개 열을 보존했다. |

UTF-8 BOM 인코딩, 쉼표 구분, 헤더 포함. 빈 칸은 결측이며 0으로 대체하지 않았다. 원본 열 이름과 값의 십진수 문자열을 그대로 보존했다. 종합 파일에는 관련성 검토 및 추적을 위한 열 외에 원본의 열전도도 등 다른 물성 열도 함께 남아 있다.

59행은 **43개 고유 반복단위 SMILES**에 해당한다. 같은 구조의 독립 계산/조건별 UUID가 있으므로 59종 고분자라고 부르면 안 된다. UUID는 두 파일에서 모두 유일하다. 반복 레코드를 임의로 평균하거나 삭제하지 않았다. 전체 파일의 `canonical_smiles_1`만으로는 셀룰로오스 치환체 등 여러 성분의 구조를 구분할 수 없으므로 이것만으로 중복 제거하지 않는다.

## 1. 10 GHz 유전율에 필요한 변수와 방법

목표 상대 복소유전율은 ε*(ω)=ε′(ω)−iε″(ω), 손실탄젠트는 tanδ=ε″/ε′이다. ε′, ε″, tanδ는 모두 무차원이다. f=10¹⁰ Hz, ω=2πf, 주기는 100 ps이다.

RadonPy의 교류 전기장 MD 방법은 다음과 같다.

1. 온도·압력에서 평형화된 고분자 셀, 원자 전하와 힘장을 준비한다.
2. E(t)=E₀ cos(ωt)의 전기장을 적용하고 셀 전체 쌍극자 μ(t), 부피 V(t), 시간축을 기록한다.
3. P(t)=μ(t)/V(t), D(t)=ε₀E(t)+P(t)를 구하고 D(t)=D₀cos(ωt−δ_D)에 적합한다.
4. ε′MD=D₀cos(δ_D)/(ε₀E₀), ε″MD=D₀sin(δ_D)/(ε₀E₀)를 계산한다.
5. 고정 전하 모델에서 빠지는 전자분극을 굴절률 n으로 보정한다: ε′corr=ε′MD+n²−1. 최종 tanδ=ε″MD/ε′corr.

따라서 주파수, 전기장 진폭·방향, μ(t), V(t), 적합 진폭과 위상, n 또는 분극률·밀도·분자량, 온도·압력, 반복 계산 및 적합 품질이 필요하다. CSV에는 결과 세 값과 일부 정적 입력만 있고, 동적 시계열·전기장·위상·행별 주파수·EFDP 적합 품질이 없다. **정적 유전율이나 분극률만으로 결측 10 GHz 값을 유일하게 계산할 수 없다.** 결측을 정적값 또는 n²으로 채우지 않았다.

코드의 입력 단위는 전기장 V/Å, 셀 부피 Å³, 셀 쌍극자 eÅ이며 분석 시 SI로 변환한다. 단량체 `qm_dipole_*`의 Debye 값은 셀 쌍극자 시계열을 대체하지 않는다. `r2`는 일반 평형 데이터의 말단간 거리 관련 열이며 EFDP 적합 R²가 아니다. `check_eq=True`도 동적 응답 적합의 성공이나 불확실성을 보증하지 않는다.

## 2. 원본의 동적 유전율 데이터와 10 GHz의 근거

전체 CSV에서 `efdp_permittivity_real`, `efdp_permittivity_imaginary`, `efdp_dielectric_loss_tan`이 모두 있는 행은 1,084개다. 제공 논문 PDF 6페이지 Fig. 3c 및 9페이지 본문은 1,084개 고분자의 **10 GHz** 결과를 기술한다. 로컬 RadonPy `AutoMD_scripts/6_ef_dp.py:44`의 기본 주파수도 10 GHz다.

그러나 논문 PDF 26페이지 Table S3는 같은 동적 속성에 대해 10–500 GHz 범위를 적고 있고, CSV에는 행별 주파수 열이 없다. 공개 `efdp_*` 열로 집계하는 원본 변환 코드도 확인되지 않았다. 따라서 출력은 다음처럼 표시했다.

- `frequency_GHz_inferred=10`: 문헌·표본 수·코드 기본값에 근거한 할당.
- `frequency_verified_per_row=False`: 행별 원시 계산 파일에서 검증한 주파수는 아님.
- `frequency_assignment`: 추정 근거 또는 동적 값 결측 상태.

이는 실무적으로 사용할 수 있는 **10 GHz 문헌 추정 데이터 집합**이며, 행별 주파수까지 확정하려면 저자의 원시 EFDP 결과/집계 코드가 추가로 필요하다. 보고된 ε′에는 n²−1을 다시 더하지 않았다. 최신 로컬 코드와 데이터 생성 당시 버전은 다를 수 있으므로 원본 `RadonPy_ver` 등 버전 열을 보존했다.

## 3. RDKit imide 판정

프로젝트 안 `.deps/imide`에 설치한 RDKit 2026.09.1로 검사했다. 공식 RDKit `fr_imide`는 `N(-C(=O))-C=O`이다. 방향족 원자 표기의 차이도 처리하기 위해 다음 산소 imide SMARTS를 사용했다.

`[#6X3](=[#8X1])[#7X3][#6X3](=[#8X1])`

한 질소가 두 carbonyl carbon과 연결된 C(=O)–N–C(=O) 모티프를 검사한다. 고리형, 비고리형, 주쇄·측쇄를 포함한다. 일반 amide 하나, urea 하나, anhydride, imidazole은 이 조건을 만족하지 않는다. acylurea/일부 ureide 등은 넓은 모티프 정의에 포함되므로 별도 보조 분류로 구분한다.

- 실제 시뮬레이션 입력 `smiles_list`에 연결점 `*`가 있고 `DP>1`인 고분자 레코드를 대상으로 했다. 소분자에서 `smiles_1`에 형식상 `*`가 있어도 DP=1이면 제외했다.
- `smiles_1`–`smiles_4`의 양의 비율 성분을 검사했다. 여러 성분이 있다고 모두 공중합체로 간주하지 않았다. 셀룰로오스와 치환체 구성도 포함된다.
- 2개의 말단 dummy와 단일결합이 있는 반복단위는 선형 3mer로 연결하여 경계에서 완성되는 모티프를 검사했다. 12,985행은 반복단위 직접 검출, 1,924행은 반복 경계에서 검출됐다.
- 복잡한 연결점용 `smiles_search` 보조 검사도 수행했으나 최종 포함 행 중 그 검사에만 의존한 행은 **0개**였다.
- 전체 입력 성분의 SMILES 파싱 실패는 0행이다. 별도 독립 점검에서 X3 제한을 완화해도 직접 모티프 추가 검출은 없었다. DP가 비어 있으면서 연결점이 있는 3행은 imide가 아님을 별도 확인했다.
- 원본 `class_PIMD=True`인 10,054행 중 9,956행이 포함됐다. 나머지 98행은 O/S를 허용하는 RadonPy 분류에는 해당하지만 산소 imide 모티프는 없는 thio/mixed-thio 계열이다. 산소 imide 추출에서는 제외하고 감사 기록에 남겼다.
- 원본 `class_PIMD=False`에서 추가된 4,953행이 있다. 원본 분류는 주쇄 중심이며 새 검사는 측쇄/구성 성분도 포함한다. `polymer_class==13` 하나로 필터하면 polyurethane/urea 우선 분류와 측쇄 구조 때문에 누락될 수 있다.
- imide 모티프가 있지만 고분자 조건을 만족하지 않는 1,753행은 제외했다.

### 넓은 imide 계열과 일반적인 carboximide 구분

`has_conventional_carboximide=True`는 두 carbonyl의 반대쪽 치환기가 탄소 또는 수소인 carboximide 모티프를 적어도 하나 갖는다는 보조 판정이다. 반복 경계에서는 연결된 구조로 검사했다. 탄소 대신 N/O/S가 연결되는 acylurea 등의 모티프만 있으면 False다. 이 값은 방향족 polyimide/주쇄 polyimide만을 뜻하지는 않는다.

| 범위 | 전체 레코드 | 동적 값 보유 레코드 |
|---|---:|---:|
| 산소 imide 넓은 모티프 전체 | 14,909 | 59 |
| `has_conventional_carboximide=True` | 11,685 | 51 |
| 넓은 모티프에만 해당 | 3,224 | 8 |

동적 파일에서 일반적인 carboximide만 원하면 `has_conventional_carboximide=True`로 필터하면 된다. 이는 35개 고유 반복단위 SMILES, 51개 계산 레코드다. 순수 thioimide까지 imide 계열에 포함하는 다른 연구 범위에서는 별도 재추출이 필요하다.

## 4. 허수부 원본 오류와 유도값

**전체 동적 데이터 1,084행 모두 원본 `efdp_permittivity_imaginary = real / tan`을 만족했다.** 최대 상대잔차는 약 4.94×10⁻¹⁶이다. 물리적 정의인 `imaginary = real × tan`을 만족하는 행은 0개다. 단순한 반올림 차이로 설명되지 않으며, 이 원본 열을 그대로 손실유전율로 사용하면 안 된다. 어느 집계 단계에서 발생했는지는 확인되지 않았다.

출력에는 다음을 분리했다.

| 열 | 해석 |
|---|---|
| `efdp_permittivity_real` | 원본 보고 실수부, 그대로 보존 |
| `efdp_permittivity_imaginary` | 원본 보고 허수부 열, 정의와 불일치하므로 직접 사용 부적절 |
| `efdp_dielectric_loss_tan` | 원본 보고 손실탄젠트, 그대로 보존 |
| `epsilon_imag_derived_real_times_tan` | 보고된 ε′×tanδ로 재구성한 허수부 |
| `raw_imaginary_definition_inconsistent` | 원본 허수부가 정의에 맞지 않음 |
| `raw_imaginary_matches_real_div_tan` | 원본이 잘못된 division 관계와 일치함 |
| `dynamic_basic_value_checks_pass` | 동적 세 값 보유, ε′>0, tanδ≥0의 기본 범위 검사만 통과 |

`epsilon_imag_derived_real_times_tan`은 **보고된 집계값의 곱**이며, 반복별 원시 ε″의 평균을 복원한 값이라고 단정할 수 없다. 실제 동적 결과의 정확도/오차막대가 검증된 것도 아니다. 59행은 모두 양의 ε′와 비음수 tanδ를 가지며 `check_eq=True`다. 그래도 원본 허수부 오류와 주파수 메타데이터 한계는 그대로 남는다.

## 5. 정적 유전율·분극률 등 관련 데이터

정적 고정 전하 MD는 전체 셀 쌍극자 요동으로 다음을 계산한다.

εstatic = 1 + (〈|μ|²〉−|〈μ〉|²)/(3ε₀kBT〈V〉)

전자분극 보정은 εstatic,corr=εstatic+n²−1이다. 따라서:

- `static_dielectric_const`: 전자분극 보정 전 정적 MD 값. 14,909행 모두 존재.
- `dielectric_const_dc`: n²−1을 더한 정적 값. 14,773행 존재.
- `qm_polarizability_monomer1`: 반복단위 분극률 부피, Å³. 14,774행 존재. xx/yy/zz/xy/xz/yz 및 성분 2–4 열도 모두 보존.
- `qm_dipole_*`: Debye 단위의 단량체 쌍극자 및 성분. 셀 전체 시간 응답과 구분.
- `qm_polarizability_sos_486*`, `_589*`, `_656*`: **nm 단위 가시광 파장**의 SOS 분극률. GHz나 10 GHz 분극률이 아니다. 파장이 없는 `_sos` 값도 원본 이름대로 보존.
- `refractive_index*`, `abbe_number_sos`: 광학 관련 값.
- `density` (g/cm³), `mol_weight_monomer*` (g/mol), `temp` (K), `press` (atm), 부피·자유부피 및 구성/버전/힘장 정보도 보존.

분극률 α(Å³), 밀도 ρ(g/cm³), 반복단위 분자량 M(g/mol)으로 Lorentz–Lorenz 값을 계산하려면 A=(4π/3)ρN_A(Σxᵢαᵢ×10⁻²⁴)/(ΣxᵢMᵢ), n=√[(1+2A)/(1−A)]를 사용한다. 이로 얻는 n²은 전자분극의 근사이며 전체 10 GHz 유전율과 동일하지 않다. 이번 출력의 `epsilon_electronic_proxy_n2` 및 `epsilon_static_corrected_derived`는 진단용 유도값이다.

## 6. 두 파일 검토 결과

- 전달 CSV를 다시 읽어 원본 CSV와 UUID별로 대조했다. 종합 파일의 **14,909행 × 원본 263열 전부 값이 정확히 일치**했다.
- UUID 중복 0개. 동적 파일 59개 UUID 모두 종합 파일에 존재하며, 종합 파일 중 동적 값 보유 행과 정확히 일치한다.
- 별도 그래프 이웃 검사로 전체 176,096행을 다시 검토했다. 이 검사는 추출 SMARTS나 3mer 생성 함수를 재사용하지 않고, N의 carbonyl 이웃 및 주기 경계 연결을 직접 검사했다. 포함 14,909개 UUID가 정확히 일치했고 누락/추가 불일치 0개였다.
- 양성·음성·반복 경계 화학 사례 12개 테스트 및 독립 사례 13개 테스트를 통과했다.
- 정적 보정 관계를 검토할 수 있는 14,773행 전부 일치했다(rtol=10⁻⁴, atol=10⁻⁶; 최대 절대잔차 9.74×10⁻⁸).
- 단량체 쌍극자 크기와 벡터 norm을 비교할 수 있는 14,841행 전부 일치했다(동일 허용오차; 최대 절대잔차 9.58×10⁻⁹).
- 유도 허수부를 다시 계산하여 모든 동적 행이 ε′×tanδ와 일치함을 확인했다.
- 추출 14,909행 모두 `check_eq=True`. 동적 59행에는 ε′≤0 또는 tanδ<0이 없다. 전체 원본 1,084행에는 부호 이상이 있지만 이번 imide 부분집합에는 없다.
- 전체 imide 중 동적 데이터 없는 14,850행은 종합 CSV에서 빈 칸으로 유지했다.

## 추가 열과 감사 자료

`source_row_number`는 헤더 제외 1부터 시작하는 원본 **레코드 순번**이다. `imide_match_atom_indices_json`은 원본 반복단위의 RDKit 0-based 원자 index로 직접 매칭한 기록이다. 경계에서만 검출되면 직접 매칭 목록은 빈 배열이며 `imide_detection_basis=linear_trimer_boundary`로 구분한다. `imide_direct_n_count`는 성분들에 대한 직접 모티프의 N 개수 합이며 중합체 전체의 기능기 수나 조성 가중 수가 아니다.

`has_small_cyclic_imide`는 원래 반복단위 또는 선형 3mer의 8원 이하 실제 고리 내 모티프 여부다. 큰 고리를 모두 분류하는 일반적인 cyclic-imide 플래그는 아니다. `class_PIMD_agrees`는 기존 플래그와 넓은 산소 모티프 여부의 일치 비교다. `rdkit_official_fr_imide_direct`는 공식 RDKit 패턴으로 원래 단위만 검사한 결과라 aromatic 표기와 반복 경계를 놓칠 수 있다.

동적 값이 없는 행의 동적 QC boolean False는 정상 판정이 아닌 미평가를 포함하므로 `has_dynamic_dielectric_data`와 함께 읽어야 한다. 정적 보정 입력이 결측인 경우에도 `static_identity_within_rtol_1e_4=False`는 미평가일 수 있다. `epsilon_static_correction_residual`이 비어 있는지 함께 확인한다.

- [validation_summary.json](validation_summary.json): 출력 및 원본 대조 결과.
- [independent_review.json](independent_review.json): 독립 구조 검사와 수식 검토.
- [column_dictionary.json](column_dictionary.json): 원본/추가 열 구분, 설명·단위, 결측 범위.
- [provenance.json](provenance.json): 원본/출력 SHA-256, 버전과 파일 정보.
- `notes/classification_audit.jsonl`: 전체 행의 구조 판정 기록.
- `notes/excluded_candidates.jsonl`: 소분자 및 thio 계열 제외 내역.
- `scripts/classify_imides.py`, `scripts/export_imides.py`, `scripts/review_final.py`: 재현 가능한 분류·추출·검토 코드.

## 출처

1. [PolyOmics 원본 데이터셋](https://huggingface.co/datasets/yhayashi1986/PolyOmics), DOI 10.57967/hf/7475, CC-BY-4.0. 이번 추출의 제어 원본은 사용자가 제공한 로컬 CSV이며 온라인 파일로 교체하지 않았다.
2. R. Yoshida, Y. Hayashi, et al., *Omics-scale polymer computational database transferable to real-world artificial intelligence applications*, [arXiv:2511.11626v1](https://arxiv.org/abs/2511.11626v1). 제공 PDF의 물리적 페이지 6, 9, 14–15, 26–27, 56–57, 61을 참조했다.
3. [RadonPy](https://github.com/RadonPy/RadonPy): 로컬 `radonpy/sim/preset/ef_dp.py:245–334`, `sim/helper.py:739–741`, `core/calc.py:1166–1201`, `core/poly.py:3926–3969`, `AutoMD_scripts/6_ef_dp.py:44`, `AutoMD_scripts/README.md`.
4. [RDKit 공식 fragment 정의](https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv), [공식 부분구조 검색 문서](https://www.rdkit.org/docs/GettingStartedInPython.html#substructure-searching). 분류 정의 자체는 위에 명시한 연구 범위와 aromaticity 확장을 포함한다.

상세 수식·단위·페이지 근거는 [방법론 근거](notes/methodology_evidence.md), 화학 판정의 경계 조건은 [구조 판정 검토](notes/imide_rule_review.md)에 기록했다.
