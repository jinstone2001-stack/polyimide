# PolyOmics 유전율 방법론 및 출처 검증

검토일: 2026-10-22. 입력 파일은 변경하지 않았다. D:/AGENTS.md, D:/Codex/AGENTS.md, 프로젝트 및 하위 경로에서 적용할 AGENTS.md는 발견되지 않았다. 아래 페이지 번호는 PDF 파일의 1-based 물리적 페이지이며, 괄호 안은 SI 인쇄 페이지 번호다.

## 가장 중요한 결론

1. `efdp_permittivity_real`, `efdp_permittivity_imaginary`, `efdp_dielectric_loss_tan`은 동적 상대유전율 관련 열이며 모두 무차원이다. 절대유전율(F/m)이 아니다.
2. 로컬 CSV의 세 동적 열에 데이터가 있는 행은 각각 1,084개다. 논문 Fig. 3c는 1,084종에 대한 10 GHz 결과라고 명시한다. RadonPy 자동 실행 스크립트 기본값도 10 GHz다. 따라서 이 표를 10 GHz 집합으로 해석할 근거는 강하지만, **CSV 자체에 행별 주파수 필드가 없고 Table S3에는 10-500 GHz라고 되어 있어 행별 주파수의 직접 검증은 불가능하다.** 출력에는 `10 GHz (문헌 및 개수 일치에 근거한 추정)`임을 표시해야 한다.
3. **원본 `efdp_permittivity_imaginary` 열은 정의와 모순된다.** 완전한 1,084행 모두 `imaginary = real / tan`을 수치 정밀도 내 만족한다. 물리적 정의는 `imaginary = real * tan`이다. 원본 열은 보존하되 손실유전율로 곧바로 사용하지 않는다.
4. `real * tan`으로 별도 유도한 값을 제공할 수 있으나, 이는 보고된 두 집계값에서 재구성한 값이다. 각 반복 계산의 허수부 평균을 복원한 값 또는 검증된 원래 측정값이라고 할 수 없다.
5. 정적 유전율, 단량체 쌍극자, 분극률, 밀도만으로 10 GHz 유전율을 유일하게 계산할 수 없다. 시간 응답/완화 스펙트럼이 추가로 필요하다.

## 10 GHz 계산에 필요한 데이터와 식

목표 주파수 f=10^10 Hz, 각주파수 ω=2πf, 주기 1/f=100 ps.

RadonPy 구현에서 전기장 E(t)=E0 cos(ωt), 전체 셀 쌍극자 μ(t), 셀 부피 V(t)를 이용한다. 논문은 sin 표현을 사용한다. 같은 시간 기준에서 위상차를 맞추면 sin/cos 기준 변경은 물리적 결과를 바꾸지 않는다.

- P(t)=μ(t)/V(t), 단위 C/m².
- D(t)=ε0 E(t)+P(t), 단위 C/m².
- D(t)를 D0 cos(ωt-δ_D)에 적합한다.
- ε′MD = D0 cos(δ_D)/(ε0 E0).
- ε″MD = D0 sin(δ_D)/(ε0 E0).
- 전자분극 보정: ε′corr = ε′MD + n² - 1.
- 최종 loss tangent = ε″MD / ε′corr.

소스 근거: `Radonpy/radonpy/sim/preset/ef_dp.py` 245-260, 322-334행. 분석은 로그의 후반 절반을 사용한다(245행). 전기장 입력값 V/Å는 ×10^10 하여 V/m로, 부피 Å³는 ×10^-30 하여 m³로, 쌍극자 eÅ는 e×10^-10 하여 C·m로 바꾼다(246-251행). 위상 δ_D는 적합된 MD 변위의 위상이며, 보정 후 loss tangent를 단순히 이 δ_D의 tan으로 대체하면 안 된다.

필수 또는 재현성에 필요한 항목:

| 항목 | 단위/의미 | all_data.csv에서의 상태 |
|---|---|---|
| 목표 주파수 f | Hz | 행별 필드 없음; 문헌 추정 |
| E0, 전기장 방향, 시간축, μ(t), V(t) | V/m, x/y/z, s, C·m, m³ | 동적 재계산에 필요한 시계열/설정 없음 |
| D0, δ_D | C/m², rad | 원시 적합 결과 없음 |
| n 또는 α, ρ, M | 무차원, Å³, g/cm³, g/mol | refractive_index, qm_polarizability_monomer*, density, mol_weight_monomer* |
| T, pressure | K, atm | temp, press |
| 평형 구조/원자 전하/힘장 | MD 입력 | CSV만으로 완전 재현 불가; forcefield/charge 등 일부 메타데이터 |
| 반복별 값, 방향별 값, fit R²/RMSPE | 품질 및 불확실성 평가 | 동적 속성에 대응하는 필드 없음 |

`r2`는 문맥 없이 EFDP fit R²로 사용하면 안 된다. 일반 평형 속성 쪽 `r2`는 고분자 말단간 평균제곱 거리 관련 데이터일 수 있으며 동적 스크립트는 별도 이름 `ef_dp_r2_10GHz_z`를 쓴다. 마찬가지로 `check_eq`는 평형 검사이지 동적 fit의 성공/품질 판정이 아니다.

논문 PDF 14-15페이지와 61페이지(SI 40): NPT, 1 atm, AC 전기장, 방향별 응답, 5개 독립 초기 구조의 평균을 설명한다. `AutoMD_scripts/6_ef_dp.py` 38-50행은 기본 T=300 K, P=1 atm, f=10 GHz, axis=z, 10주기이며 모두 사용자 설정으로 변경될 수 있다. 따라서 이 기본값을 모든 기존 CSV 행의 실제 실행값이라고 단정하지 않는다.

## 주파수 식별의 증거 강도

- PDF 6페이지 Fig. 3c 캡션: 1,084 polymers at 10 GHz. 그림과 캡션을 렌더링하여 직접 확인했다.
- PDF 9페이지 본문: 1,084 isotropic amorphous polymers의 10 GHz 분포를 설명한다.
- PDF 26페이지(SI 5), Table S3: 세 `efdp_*` 열의 수가 각각 1,084이며 범위는 10-500 GHz라고 기재한다. 렌더링하여 직접 확인했다.
- PDF 14페이지 및 61페이지: 워크플로는 10-500 GHz를 계산할 수 있다.
- `Radonpy/AutoMD_scripts/6_ef_dp.py:44`: `RadonPy_EFDP_Freq` 환경변수 기본값 10.0e9 Hz.
- `huggingface/README.md:334-336`: 세 열 이름은 나오지만 그 열 자체의 주파수는 지정하지 않는다.
- CSV 헤더에는 frequency/freq/GHz, EFDP 방향/전기장/fit 파라미터가 없다.

권장 레이블: `frequency_GHz_inferred=10`, `frequency_verified_per_row=False`, `frequency_evidence=paper_Fig3c_1084_matches_CSV_and_RadonPy_default`.

현재 소스 출력 이름은 `ef_dp_dynamic_dielectric_constant_corr_10GHz_z` 등이며, 공개 CSV의 간략한 `efdp_*` 이름으로 집계/변환한 코드는 로컬 Python/Markdown 소스에서 찾지 못했다. 주파수·방향·반복 평균·보정 매핑을 끝까지 입증하려면 원 저자의 집계 코드 또는 원시 EFDP 결과 파일이 필요하다. 논문 정의는 전자분극을 포함하므로, 이미 보고된 `efdp_permittivity_real`에 n²-1을 다시 더하지 않는다.

## 허수부 열의 독립 수치 검증

전체 CSV에서 세 열만 읽고 세 값이 모두 있는 행을 조사했다.

| 검사 | 결과 |
|---|---:|
| 완전한 동적 데이터 행 | 1,084 |
| imaginary ≈ real / tan 만족 | 1,084 / 1,084 |
| imaginary ≈ real * tan 만족 | 0 / 1,084 |
| division 관계 최대 상대잔차 | 4.936502254564192e-16 |
| 원본 real 범위 | -1.701949715 ~ 5.064715193 |
| 원본 imaginary 범위 | -412115.8696169088 ~ 1063555.0242537311 |
| 원본 tan 범위 | -0.063236696 ~ 0.234346329 |

근접 검사 rtol=1e-6, atol=1e-10. 평균/반올림 정도의 불일치가 아니라 모든 값이 division 식과 일치하는 현상이다. 다만 어떤 데이터 처리 단계에서 생긴 문제인지는 이 증거만으로 확정할 수 없다.

예: ε′=3.157339427, 원본 ε″=359.5012416090762, tan=0.008782555이면 정의상 재구성 ε″=0.0277295072 수준이다. 원본은 ε′/tan에 해당한다.

원본 열과 유도 열을 분리하고 `raw_imaginary_definition_inconsistent=True` 등의 QC 필드를 둔다. ε′≤0, tan<0 등은 별도 검토 대상으로 표시한다. 부호 이상은 무조건 삭제하거나 절대값을 취하지 말고 원본을 보존한다. 특히 빈 동적 데이터에 정적값을 대신 넣지 않는다.

## 정적 유전율, 전자분극, 광학 속성

PDF 57페이지(SI 36), 식 (16)-(17):

εstatic = 1 + (〈|μ|²〉 - |〈μ〉|²)/(3 ε0 kB T 〈V〉)

εstatic,corr = εstatic + n² - 1

- `static_dielectric_const`: 고정 전하 MD의 쌍극자 요동 기반 값. 전자분극 보정 전.
- `dielectric_const_dc`: 전자분극 n²-1을 더한 정적 값.
- 코드 `Radonpy/radonpy/sim/helper.py:739-741`, 설명 `Radonpy/AutoMD_scripts/README.md:226-227`이 이 관계를 명시한다.
- `Radonpy/radonpy/sim/lammps.py:2125-2134`는 mdtraj.static_dielectric을 사용하며, 2593행이 결과 필드를 저장한다.
- n²은 전자분극 기여의 근사이며, 전체 10 GHz 상대유전율과 동일하다고 가정할 수 없다.

Lorentz-Lorenz 식의 구현(`Radonpy/radonpy/core/calc.py:1166-1201`):

A = (4π/3) ρ NA [Σ xi αi ×10^-24]/[Σ xi Mi]

n = sqrt[(1+2A)/(1-A)]

여기서 ρ는 g/cm³, αi는 Å³, Mi는 g/mol이고, xi는 정규화된 반복단위 비율이다. SI 변환 없이 α(Å³)를 그대로 대입하면 안 된다. 논문 PDF 56페이지 식 (14)는 해당 관계를 설명하지만, 표기에서 NA 등이 축약되어 있으므로 수치 구현에는 단위가 명시된 소스를 따른다.

분극률은 단량체/반복단위의 polarizability volume(Å³)이며 SI 분극률(C·m²/V)로 쓸 때 αSI=4πε0 αvolume×10^-30 변환이 필요하다. 근거: `Radonpy/radonpy/sim/psi4_wrapper.py:932` 및 PDF 27페이지 Table S3. 단량체 쌍극자는 Debye이며 전체 셀 쌍극자 시계열과 다르다(`psi4_wrapper.py:1222-1237`).

`qm_polarizability_sos_486*`, `_589*`, `_656*`는 가시광 파장 nm에 대한 TD-DFT/SOS 분극률이다(PDF 27페이지 Table S3). 이 숫자는 GHz 주파수가 아니며 10 GHz 값을 의미하지 않는다. 파장 없는 `*_sos` 열은 명시적인 파장 지정 없이 보존한다. 재해석하여 10 GHz 분극률로 이름을 붙이지 않는다.

## 결과 작성 권고

- 좁은 파일: imide 구조로 판정되고 동적 값이 있는 행, 원본 동적 세 열, 정의에 따른 ε″ 유도 열, 주파수 추정 근거, 원본 이상/QC, T/P, 주요 전자분극 입력, 식별자와 SMILES.
- 넓은 파일: imide 구조를 가진 모든 고분자 행. 정적/동적 유전율, dipole 벡터와 크기, polarizability 텐서/광학 분산, refractive index/Abbe, 밀도/분자량/구조와 관련 메타데이터. 동적 값이 없는 행은 결측으로 보존.
- 두 파일 모두 문헌 계산값임을 명시한다. 실험 데이터가 아니다. 논문은 13개 polyimide의 실험값을 0% RH로 외삽하여 비교했다(PDF 15페이지); 이 사실이 전체 CSV의 습도·실험 조건을 의미하지 않는다.
- 원본을 유지하면서 유도 열에 접두사/접미사를 써서 출처값과 구분한다. 독립적으로 평균한 집계값 간 비선형 연산은 개별 시뮬레이션 수준 재계산과 같지 않다.

## 사용한 원본

- `D:/Codex/03.PolyOmic/2511.11626v1.pdf` (arXiv:2511.11626v1)
- `D:/Codex/03.PolyOmic/huggingface/README.md`
- `D:/Codex/03.PolyOmic/huggingface/all_data.csv`
- 위에 명시한 `D:/Codex/03.PolyOmic/Radonpy/` 소스 파일
- 공개 원본: https://huggingface.co/datasets/yhayashi1986/PolyOmics , https://github.com/RadonPy/RadonPy

이번 검토는 해당 로컬 스냅샷을 기준으로 한다. 현재 RadonPy 코드와 원본 데이터의 생성 당시 버전은 동일하다고 보장하지 않는다.
