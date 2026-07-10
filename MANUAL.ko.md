# ShereKhan 단계별 매뉴얼 — Individual vs Global fit + Jackknife 검증

> 🌐 **한국어** · [English](MANUAL.en.md)

CPMG 완화 분산 데이터에서 **모든 잔기가 하나의 교환 과정(kex)을 공유하는지(global fit)**,
아니면 **잔기마다 자기 kex를 갖는지(individual fit)**를 판정하고, **잭나이프**로 global fit의
견고성을 검증하는 실전 절차.

- 대상: `sk-prepare` → `sk-run` 파이프라인 사용자
- 사전 지식: `README.md`(설치), `USAGE.md`(config/모델 상세)
- 이 문서는 원시 `.dat`에서 판정·검증까지의 **순서**를 다룬다.

---

## 0. 설치 (한 번만)

```bash
pip install .
```

콘솔 스크립트 3개가 설치된다: `sk-synth`, `sk-prepare`, `sk-run`.
설치 없이 실행하려면 `python sk_prepare.py ...`, `python sk_run.py ...` 형태로도 된다.

요구사항: Python ≥ 3.8, `numpy`, `scipy`, `matplotlib`.

---

## 1. 입력 데이터 준비 (`.dat`)

자기장(field)마다 파일 1개. 형식:

```
 60.12                        # field (MHz, 1H Larmor)
0.040000                      # tcp (s, CPMG 총 echo time)
#nu_cpmg(Hz)  R2(1/s)  Esd(R2)
# R1                          # 잔기 라벨
  50.000  35.80  0.36
 100.000  35.14  0.35
 ...
# R2
 ...
```

- 같은 잔기는 **모든 파일에서 같은 라벨**이어야 field 간에 묶인다.
- kex를 잘 결정하려면 **field 2개 이상** 필요 (잭나이프·비교도 field ≥ 2에서 의미).
- 테스트 데이터가 없으면 `sk-synth`로 합성 데이터 생성.

---

## 2. Step 1 — config 생성 (`sk-prepare`)

```bash
sk-prepare f60.dat f90.dat > run.conf
```

`sk-prepare`가 하는 일:
- 잔기별 alpha 계산 → fast/slow 교환 체제 추정, 모델 제안(fast→Meiboom, slow→London, 애매→Matrix)
- **`compare_aic: true`, `jackknife: true`를 기본으로 넣는다** (개별-전역 비교 + 잭나이프 활성)

출력은 stdout JSON. 파일로 저장 후 다음 단계에서 검토.

---

## 3. Step 2 — config 검토·수정 (`run.conf`)

핵심 필드:

| 필드 | 값 | 의미 |
|------|-----|------|
| `exchange` | `"fast"` / `"slow"` / `"undefined"` | 교환 체제; `"undefined"`(애매하거나 단일 필드 데이터)는 모델 `Matrix`와 짝이 되고 모든 잔기를 자동으로 `flag: "off"`로 둠 — 수동으로 켜야 함 |
| `model` | `"Meiboom"` / `"London"` / `"Matrix"` | 피팅 모델 (Matrix=정확, 느림) |
| `init.mode` | `"guess"` / `"values"` | 초기값: 그리드 탐색 / 직접 지정 |
| `residues[].flag` | `"on"` / `"off"` | 잔기 포함 여부 |
| `compare_aic` | `true` / `false` | 잔기별 개별 vs 전역 비교 |
| `jackknife` | `true` / `false` | leave-one-residue-out 검증 |

검토 포인트:
- 분석할 잔기는 `flag: "on"`인지 확인. 명백한 잡음 잔기는 `"off"`.
- `compare_aic`/`jackknife`가 켜져 있는지 확인 (판정·검증의 핵심 출력).
- 초기값을 알면 `init.mode: "values"`로 그리드 탐색을 건너뛰어 빠르게. 예:
  ```json
  "init": {"mode": "values", "kex": 1000.0, "pB": 0.02, "csd": 2.5}
  ```
  (Meiboom은 `pB`/`csd` 대신 `"phi"` 지정.)

---

## 4. Step 3 — 실행 (`sk-run`)

```bash
sk-run run.conf
```

산출물:
- `<Project Name>.log` — 사람이 읽는 전체 리포트 (전역 fit 통계 + 비교표 + 잭나이프)
- `<Project Name>.pdf` — 잔기별 분산 곡선 플롯 (exp 점 + calc 곡선)
- **stdout** — JSON 블록 (프로그램 연동용, §9)

---

## 5. Step 4 — 전역(global) fit 결과 읽기 (`.log`)

`.log`의 통계 블록:

```
npar=13 nvar=66 ndof=53 chi2=  63.629 chi2/dof=   1.201
kex: 1051.054  +-    94.327
 pB:    0.020  +-     0.002
resId:  csd [ppm]   R2_0 (60.1 MHz)   R2_0 (90.2 MHz)
  R1:   7.669 ...
```

- `chi2/dof ≈ 1` 이면 전역 fit이 데이터·오차와 일관. ≫ 1 이면 모델 부적합 또는 오차 과소평가.
- `kex ± std`, `pB ± std`가 공유 교환 파라미터. 잔기별 `csd`(또는 phi)와 `R2_0`가 개별 파라미터.

이 블록은 **전역 fit** 결과다. 다음 단계에서 개별 fit과 비교한다.

---

## 6. Step 5 — 개별 vs 전역 판정 읽기 (핵심)

`.log` 하단 / stdout에 비교 블록:

```
Model               chi2      k      n          AIC         AICc
global            86.139     16    110      118.139      123.989
individual        83.084     20    110      123.084      132.523
delta AICc (individual - global):        8.534
Preferred model (whole dataset): global
---
Per-residue preference (individual vs global fit):
 resId    chi2_glob   chi2_indiv        dAICc       z(kex)     better
   R1        23.31        21.12          4.23         1.55     global
   R5      1164.47        11.38      -1146.68        26.31 individual
---
Residues preferring global: 4   preferring individual: 1
```

**전체(whole dataset) 판정** — `Preferred model`:
- `delta AICc > 0` → **global** (공유 kex가 파라미터 대비 낫다)
- `delta AICc < 0` → **individual** (잔기마다 kex 필요)

**잔기별 판정** — `better` 열:
- `chi2_glob` = 공유 kex 하에서 그 잔기의 chi² 기여분
- `chi2_indiv` = 그 잔기만 단독 피팅했을 때 chi²
- `dAICc` = AICc(individual) − AICc(global): **> 0 → global**, **< 0 → individual**
  - 규약: individual은 자기 교환 파라미터 비용을 지불, global은 공유 kex를 무료로 빌린다.
    → "이 잔기가 자기만의 kex 값어치가 있나?"
- `z(kex)` = (kex_individual − kex_leave-one-out) / σ(kex_individual)
  - **leave-one-out 기준**(그 잔기 뺀 나머지로 적합한 kex)이라 편향 없음.
  - **|z| > 2** → 그 잔기가 공유 kex와 통계적으로 불일치 (outlier 후보).

읽는 순서: ① 전체 판정으로 global/individual 큰 그림 → ② `better` 열로 어느 잔기가 갈리는지
→ ③ 갈리는 잔기의 `z(kex)`·`chi2_glob` 크기로 심각도 판단.

---

## 7. Step 6 — 잭나이프로 global fit 검증

```
Jackknife validation of the global fit (leave-one-residue-out)
kex (full global fit):            1134.895
kex (jackknife mean):             1133.920
jackknife std error (kex):          67.180
relative std error:                  5.92%
---
 resId      kex(-res)      delta_kex        |z|
   R1        1166.713        -31.818       0.47
   R4        1179.261        -44.366       0.66
---
(* = influential: dropping this residue shifts kex by > 2 jackknife SE)
```

각 잔기를 하나씩 빼고 전역 재적합 → 공유 kex 변화 추적.

- **relative SE** 작음(예: < 몇 %) + `*` 플래그 잔기 없음 → **global fit 견고**. 어느 한 잔기도
  kex를 지배하지 않음.
- **relative SE 큼** (예: > 15–20%) → kex가 불안정. 잔기 간 kex 불일치 신호.
- `delta_kex` = kex_full − kex(그 잔기 제외). `|delta_kex| > 2 SE`면 `*`로 표시 → 영향력 큰 잔기.

> 주의: 강한 outlier가 **하나** 있으면 그 outlier가 SE 자체를 부풀려 `*` 플래그가 가려질 수
> 있다(single-deletion 잭나이프 breakdown point = 1). 이 경우 §6의 `z(kex)`·`chi2_glob`이 더
> 민감하다. relative SE는 전반적 안정성 지표로 여전히 유효.

---

## 8. Step 7 — 이질성 진단: 개별 잔기(N개) & 그룹 구조(K개)

전체 판정이 individual로 나왔다고 끝이 아니다. **왜** 그런지 두 가지로 파고든다.

> **핵심 원리**: 최소자승 global fit은 outlier에 **비robust**하다. 소수의 발산 잔기 하나만
> 있어도 그 잔기가 전체 chi²를 지배해 공유 kex를 끌어당기고, 그러면 정상 잔기들까지 오염된
> global kex에 안 맞아 **전부 `better = individual`로 뒤집힌다**. 그래서 이진 `preferred`만
> 보면 원인을 못 찾는다. **`chi2_glob`의 크기(랭킹)**가 오염에 강한 진짜 신호다.

### 8-A. 개별 individual 잔기 진단 (N개) → 제외 → 재검증

**① 범인 지목 (chi2_glob 랭킹)** — 비교표를 `chi2_glob` 내림차순으로 보면, 참 individual
잔기들이 상위에 몰린다. individual이 여러 개(N개)여도 상위 N개가 그들이다. 오염이 심한
경우(예: 6개 중 4개가 individual)에도 랭킹은 유지된다. `z(kex)`도 보조 신호.

```
ranked chi2_global:  R6=1241*  R3=659*  R4=298*  R5=288*  R2=82  R1=56   (* = 참 individual)
```
→ 상위 4개(R6,R3,R4,R5)가 individual. 정상 잔기(R1,R2)와 크게 벌어짐.

**② 제외 후 재적합** — 지목한 잔기들의 `flag`를 `"off"`로:
```bash
# run.conf에서 individual 잔기들의 "flag": "on" -> "off"
sk-run run.conf
```

**③ 재검증** — 남은 잔기들이 진짜 하나의 공유 kex였다면:
- 전체 판정이 **global**로 복귀, 남은 잔기 전부 `better = global`
- 잭나이프 relative SE 급감 (예: 12.6% → 0.6%, N개 제거 시에도 ≤ ~3%)
- `kex`가 참 합의값으로 수렴

### 8-B. 여러 global 그룹(K개) 구조 판별

제외해도 global이 복원되지 않고, **fit된 kex가 몇 개의 값에 뭉치면** — 데이터가 하나의
global fit이 아니라 **K개의 서로 다른 global fit**으로 이뤄졌을 수 있다(예: 잔기들이 kex≈600
그룹과 kex≈1800 그룹으로 나뉨). 절차:

**① 그룹 가설** — 잔기별 개별 fit kex(비교표의 `fit_kex`/individual `kex`)를 보고 뭉치는
값으로 잔기를 K개 그룹으로 나눈다.

**② 그룹별 global fit** — 각 그룹만 `flag: on`으로 둔 config를 만들어 따로 실행:
```bash
# groupA.conf : 그룹 A 잔기만 on   /   groupB.conf : 그룹 B 잔기만 on
sk-run groupA.conf
sk-run groupB.conf
```
각 그룹이 **`better` 전부 global** + **잭나이프 relative SE 작음**이면 그 그룹은 내부적으로
일관된 global fit이다.

**③ 구조 선택 (AICc)** — 세 가설의 총 AICc를 비교한다(같은 데이터 점 수 n):

| 구조 | chi² | 파라미터 k |
|------|------|-----------|
| 1-global | 전체 global fit chi² (큼) | 1 kex + 로컬 |
| **K-global** | Σ(그룹별 global chi²) | K kex + 로컬 |
| individual | Σ(잔기별 개별 chi²) | N kex + 로컬 |

`AICc = chi² + 2k + 2k(k+1)/(n−k−1)`, **가장 낮은 것이 최적 구조**. 여기서 `n`은 **전체
데이터 점 수**(모든 잔기 합, 그룹별 n이 아님) — 세 가설이 같은 데이터를 설명하므로 n은 동일.
각 그룹 config의 `Model comparison` 표에서 `global`의 `chi2`/`k`를 읽어 합하면 K-global의
chi²/k가 된다. K-그룹 데이터에서는 K-global이 1-global(과소적합)과 individual(과대적합)
둘 다 이긴다 — individual은 chi²가 조금 더 낮아도 교환 파라미터 N개 값을 치러 페널티를 받는다.

```
1-global     chi2=  2955.88  AICc=  3003.45
2-global     chi2=    91.69  AICc=   144.97   <-- WINNER
individual   chi2=    87.57  AICc=   165.99
```
→ "이 데이터에는 global fit이 2개 있다"가 통계적으로 확정.

> 8-A / 8-B 절차는 `demo/`(단일 outlier)와 `demo/multi/`(N개 individual, K개 그룹)에서
> 자동 재현된다(§10).

---

## 9. stdout JSON 연동 (프로그램 처리)

`sk-run` stdout에 라벨된 JSON 블록:

```
##### model_comparison
{ ... compareModelsAIC 결과 (per_residue 포함) ... }
##### jackknife
{ ... jackknifeGlobal 결과 ... }
#####
[ ... 잔기별 최종 fit 값 (reportAllValues) ... ]
```

- `##### model_comparison` / `##### jackknife` 블록은 **strict JSON** (비유한값은 `null`).
- **마지막 `#####` 다음 줄이 잔기 결과 리스트** — 기존 파서 호환을 위해 항상 마지막.
  프로그램은 "마지막 `#####` 뒤 줄"을 읽어야 안전(`split('#####')[-1]`).
- `compare_aic`/`jackknife`가 꺼져 있으면 해당 블록은 출력되지 않는다.

---

## 10. 빠른 검증 (선택)

기능이 ground-truth를 맞게 판정하는지 확인:

```bash
bash demo/run_demo.sh          # global / individual / mixed (단일 outlier)
bash demo/multi/run_multi.sh   # N개 individual (m1..m4), K개 global 그룹 (g2, g3)
```

- `demo/run_demo.sh` — 전역/개별/혼합 3개 시나리오 생성·실행·평가 (exit 0 = 전부 통과).
- `demo/multi/run_multi.sh` — 8-A(individual 1~4개 식별·복원), 8-B(2·3개 global 그룹의
  AICc 구조 선택) 자동 재현.
- 각 `demo/README.md`, `demo/multi/README.md`에 시나리오·기대 결과 설명.

---

## 11. 판정 규칙 요약

```
전체 판정 (whole dataset):
  delta AICc (individual - global) > 0  ->  GLOBAL
                                    < 0  ->  INDIVIDUAL

잔기별 (better 열):
  dAICc > 0  ->  global      |  dAICc < 0  ->  individual
  |z(kex)| > 2               ->  공유 kex와 불일치 (outlier 후보)

잭나이프:
  relative SE 작고 * 없음     ->  global fit 견고
  * 표시 잔기                 ->  영향력 큰 잔기

이질성 진단 (전체 individual일 때):
  chi2_glob 랭킹 상위 N개  ->  individual 잔기 (오염에 robust; §8-A)
    -> flag off -> 재적합 -> global 복원 확인
  fit kex가 K개 값에 뭉침   ->  그룹별 fit + AICc(1-global vs K-global vs individual)
    -> 최저 AICc가 K-global이면 "global fit K개" 구조 확정 (§8-B)
```

---

## 12. 트러블슈팅

| 증상 | 원인 / 대응 |
|------|-------------|
| `Two or more fields strengths are necessary for alpha calculation` | field 1개 → 비교/알파 불가. field ≥ 2 데이터 필요. |
| `Jackknife skipped: need at least 2 active residues.` | 활성 잔기 < 2. 잭나이프는 2개 이상 필요. |
| `covariance matrix is None` 경고 | fit이 singular/부족 제약. std 0으로 표시. 초기값·모델·데이터 점검. |
| `n is not an integer` | tcp × nu_CPMG가 정수 아님. `.dat`의 tcp/주파수 확인. |
| Matrix 모델이 느림 | 잭나이프는 N회 재적합. 큰 데이터는 Meiboom/London 또는 `init.mode: values`로 가속. |
| 전체 individual인데 원인 불명 | §8 절차로 `z(kex)`/`chi2_glob` 큰 잔기 outlier 의심 → 제외 재검증. |
