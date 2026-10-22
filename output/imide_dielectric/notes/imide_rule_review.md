# Imide structure classification — independent review

검토일: 2026-10-22 (Asia/Seoul). 원본 자료는 변경하지 않았으며, 프로젝트와 상위 `D:/Codex`, `D:/`에서 적용되는 AGENTS.md는 발견되지 않았다.

## 권장 범위와 SMARTS

이번 추출의 `imide 계열`은 반복단위 또는 반복단위 연결부에 **하나의 질소가 두 carbonyl carbon에 연결된 `C(=O)-N-C(=O)` 구조를 포함하는 고분자**로 정의한다. 고리형 및 비고리형, 주쇄 및 측쇄 imide를 포함하되 각각의 근거/범위는 별도 열로 남긴다. 이는 IUPAC의 diacyl nitrogen 정의 및 RDKit의 공식 `fr_imide` descriptor에 근거한 구조 기반의 넓은 분류이며, 모든 행이 상용 통상 의미의 aromatic polyimide인 것은 아니다.

- 공식 RDKit `Data/FragmentDescriptors.csv`: `fr_imide` = `N(-C(=O))-C=O`.
- 권장 aromaticity-tolerant oxygen-imide SMARTS: `[#6;X3](=[#8;X1])[#7;X3][#6;X3](=[#8;X1])`.
- 원자 번호 `#6`, `#7`는 C/c 및 N/n 표기에 독립적이다. 생략된 SMARTS bond는 single 또는 aromatic을 허용한다. 불필요하게 모든 결합형을 허용하는 `~`보다 이 표현이 명확하다.
- 원본 RadonPy `PIMD` SMARTS: `[X3&C,X3&c](=[O,S])[X3&N,X3&n][X3&C,X3&c](=[O,S])`.
- 따라서 RadonPy에는 O/S mixed 및 thioimide도 포함된다. `O=C-N-C=O`와 `O/S=C-N-C=O/S`를 각각 검사하고, 산소 imide가 없지만 S 허용 규칙만 통과하는 행은 별도 `thioimide_only`로 보고한다. 산소 imide와 구분 없이 섞지 않는다.

RDKit는 판정 엔진과 SMARTS의 공식 예제를 제공하지만, 이 프로젝트의 고분자 범위와 aromaticity 확장 규칙 자체는 분석상의 명시적 정의이다. `MolFromSmiles`로 sanitize된 분자에 `MolFromSmarts` 및 `GetSubstructMatches`/`HasSubstructMatch`를 적용하며, 문자열 검색을 최종 판정에 사용하지 않는다.

## 중요한 경계 조건

1. **일반 amide/urea/anhydride와의 구별.** 하나의 N에 두 carbonyl carbon이 직접 연결되어야 한다. 단일 amide `CC(=O)NC`, urea `NC(=O)N`, anhydride `CC(=O)OC(C)=O`는 이 규칙을 통과하지 않는다. 두 carbonyl이 분자 어딘가에 존재한다는 사실만으로 통과시키면 안 된다.
2. **Acylurea / cyclic ureide.** 예를 들어 `CC(=O)NC(=O)N`은 중간 N에 두 carbonyl이 결합하므로 RDKit broad motif를 통과한다. 이는 broad imide-like 구조 범위 안에서는 알고리즘 오류가 아니지만, 전통적인 polyimide만 원하는 후속 분석에서는 구분이 필요하다. Hydantoin/barbiturate류도 imide/urea 특성이 함께 나타날 수 있다. 필요하면 carbonyl의 반대쪽 치환 원자가 carbon/H인지, N/O/S인지 기록한다.
3. **반복 경계.** `*N(C)C(=O)CC(=O)*`처럼 반복단위의 N과 다음 반복단위 carbonyl 사이 결합으로 imide가 완성될 수 있다. 원본 한 단위만 검사하면 누락된다. 정확히 두 terminal dummy가 있고 각 dummy의 degree가 1인 구조는 head–tail 선형 3mer를 만들어 경계부까지 검사한다. dummy 연결 bond order를 보존하고 결합 후 sanitize한다. `*`는 실제 탄소가 아니므로 탄소로 치환하여 검사하지 않는다.
4. **고리형 판정.** 원본 반복단위 또는 선형 oligomer에서 N 및 두 carbonyl C가 동일한 실제 ring에 속하는지 검사한다. `smiles_search`나 RadonPy가 만든 cyclic oligomer는 인위적인 polymer macrocycle을 포함하므로 단순 `IsInRing()`으로 cyclic-imide를 단정하면 안 된다.
5. **Copolymer.** `smiles_1`–`smiles_4`와 대응 `copoly_ratio_i`를 확인한다. 양의 비율인 구성 반복단위에서 motif가 있으면 imide-containing copolymer에 포함한다. 비율 0인 단위만의 match는 포함 근거가 될 수 없다. 서로 다른 repeat unit 사이에서만 motif가 생기는 경우는 sequence/head–tail 연결 정보가 있어야 확정할 수 있으며, 단독 monomer match와 분리하여 기록한다. 누락된 비율을 임의로 0이나 1로 확정하지 않는다.
6. **소분자 제외.** imide motif만으로 고분자를 판정하지 않는다. 반복 연결점, DP, 데이터셋 source/category 등 원본 문맥으로 polymer 여부를 별도로 확인한다. 예를 들어 succinimide는 motif positive이나 반복점 없는 소분자이다. 다만 dummy 0개라는 이유만으로 명확한 polymer metadata를 가진 fully specified oligomer를 즉시 버리지 않는다.
7. **비표준 연결점.** dummy 1개, 3개 이상, disconnected fragments, 잘못된 valence/parse failure는 임의 결합하지 않고 review 상태로 남긴다. ladder의 다중 연결점을 2개처럼 처리하면 구조를 왜곡한다.
8. **Charged/tautomeric representations.** neutral carbonyl-imide SMARTS는 `C(=O)-[N-]-C(=O)`처럼 X2인 이온형 또는 lactim tautomer를 놓칠 수 있다. 해당 표현이 있는지 별도 점검하고, 자동 tautomerization/중성화는 원본 chemical state를 바꿀 수 있으므로 무조건 시행하지 않는다.

## 기존 분류 열은 검증 보조 자료

`Radonpy/radonpy/core/poly.py:3926–3969`의 `polyinfo_classifier`는 주쇄를 추출하고 cyclic 4mer를 생성하여 functional-group flags를 구한다. 따라서 sidechain-only imide는 `class_PIMD`와 불일치할 수 있다. 같은 파일의 약 4055행부터 최종 `polymer_class` 선택 우선순위가 `PURT` → `PURA` → `PIMD`이므로 **`polymer_class == 13`만 사용하는 추출은 imide를 가진 혼합 기능성 고분자를 누락할 수 있다**. 기존 `class_PIMD`/`polymer_class`는 원본값을 보존하고 독립 RDKit 판정과 cross-tab으로 비교해야 한다.

권장 최소 감사 열: source row number/UUID, 원본 smiles 열들, canonical repeat SMILES, parse 상태, dummy 수, 직접 motif 여부, 경계 motif 여부, matched N 수, oxygen-imide/thio-only 구분, 기존 class_PIMD, 기존 polymer_class, polymer 근거, 포함 이유/검토 사유. 서로 대칭인 match가 여러 개 나와도 imide 중심 N index를 deduplicate하면 기능기 수를 안정적으로 셀 수 있다.

## 최종 검토용 예제

다음 사례 및 aromatic uracil 표현 `O=c1cc[nH]c(=O)[nH]1`을 독립 스크립트 `verify_imide_examples.py`로 실행하였다. RDKit **2026.09.1**, 총 **13건 통과**. 결과는 같은 폴더 `verify_imide_examples.json`에 저장했다. aromatic uracil 사례는 확장 oxygen SMARTS=True이지만 공식 `Fragments.fr_imide`=0이므로 N/n 및 C/c 처리를 확장하는 실제 이유를 확인했다. 이 결과는 분류법의 사례 검증이며 전체 데이터셋 검증은 별개이다.

| 분류 | SMILES | oxygen motif 예상 | polymer 포함 |
|---|---|---:|---|
| Succinimide | `O=C1CCC(=O)N1` | True | False, 소분자 |
| N-methylphthalimide | `CN1C(=O)c2ccccc2C1=O` | True | False, 소분자 |
| Acyclic diacetamide | `CC(=O)NC(C)=O` | True | False, 소분자 |
| Amide | `CC(=O)NC` | False | False |
| Urea | `NC(=O)N` | False | False |
| Anhydride | `CC(=O)OC(C)=O` | False | False |
| Imidazole (imide 아님) | `c1ncc[nH]1` | False | False |
| Thioimide only | `CC(=S)NC(C)=S` | False; S-rule True | False, 소분자 |
| Acylurea broad motif | `CC(=O)NC(=O)N` | True | False, 소분자 |
| 경계에서 완성 | `*N(C)C(=O)CC(=O)*` | 단위 False, trimer True | True |
| Pendant imide polymer | `*CC(*)N1C(=O)CCC1=O` | True | True, pendant |

## 확인한 1차 근거

- [RDKit 공식 fragment descriptor source](https://github.com/rdkit/rdkit/blob/master/Data/FragmentDescriptors.csv), GitHub plugin fetch 기준 blob SHA `d3de26bcd432d64bce3c306f19dd87f46afb0307`.
- [RDKit 공식 substructure searching 설명](https://www.rdkit.org/docs/GettingStartedInPython.html#substructure-searching).
- [RDKit Book: SMARTS reference 및 aromaticity](https://www.rdkit.org/docs/RDKit_Book.html#smarts-reference).
- [IUPAC Gold Book: imides](https://goldbook.iupac.org/terms/view/I02948/plain).
- [RadonPy polymer classifier source](https://github.com/RadonPy/RadonPy/blob/develop/radonpy/core/poly.py); 실제 라인 검토는 이 프로젝트의 `Radonpy/radonpy/core/poly.py`에 근거한다. GitHub API에서 default branch는 `develop`로 확인했다. local snapshot과의 일치 여부는 별도 provenance 검사에서 확인해야 한다.
