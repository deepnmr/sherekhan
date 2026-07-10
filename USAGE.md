# ShereKhan 1.2.0 사용법

NMR CPMG 분산 곡선 분석 도구.  
저자: Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)  
원본 코드: Adam Mazur, Bjoern Hammesfahr, Christian Griesinger, Donghan Lee, Martin Kollmar (Max-Planck-Institute for Biophysical Chemistry, 2012)  
Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)

---

## 빠른 시작 (Quick Start)

```bash
# 1. 합성 데이터 생성 (테스트용)
python sk_createSyntheticDataset.py

# 2. config 파일 생성
python sk_prepare.py synth-60.dat synth-90.dat > run.conf

# 3. config 수정 (exchange, model, residue flag 등 설정)
#    → run.conf 직접 편집

# 4. 분석 실행
python sk_run.py run.conf
```

---

## 스크립트 설명

### `sk_createSyntheticDataset.py` — 합성 데이터 생성

테스트용 CPMG 데이터셋을 생성한다.

```bash
python sk_createSyntheticDataset.py
```

**출력 파일:**
- `synth-60.dat` — 60.12 MHz 자기장 데이터
- `synth-90.dat` — 90.23 MHz 자기장 데이터

**파라미터 수정 (스크립트 내부):**

| 변수 | 기본값 | 설명 |
|------|--------|------|
| `fields` | `[60.12, 90.23]` | 1H 자기장 강도 (MHz) |
| `tcp` | `0.04` | CPMG pulse delay (s) |
| `kAB` | `20.0` | A→B 교환 속도 (s⁻¹) |
| `kBA` | `1000.0` | B→A 교환 속도 (s⁻¹) |
| `noiseRatio` | `0.02` | 노이즈 비율 (2%) |
| `res` | (잔기 목록) | 분석할 잔기와 화학적 이동 차이 (csd, ppm) |

---

### `sk_prepare.py` — Config 파일 생성

실험 데이터 파일(.dat)을 읽어 분석용 config(JSON)를 자동 생성한다.

```bash
python sk_prepare.py <field1.dat> <field2.dat> [field3.dat ...] > myproject.conf
```

**예시:**
```bash
python sk_prepare.py exp-600MHz.dat exp-900MHz.dat > project.conf
```

**출력:** JSON 형식의 config (stdout). 파일로 리다이렉트하여 저장 후 편집.

**자동 판별:** 각 잔기의 alpha 값을 계산하여 fast/slow exchange 체제를 추정하고 적합한 model을 제안한다.

---

### `sk_run.py` — 분석 실행

Config 파일을 읽어 CPMG 분산 곡선을 피팅하고 결과를 출력한다.

```bash
python sk_run.py <config.conf>
```

**예시:**
```bash
python sk_run.py project.conf
```

**출력 파일 (`Project Name` 기준):**
- `{name}.log` — 피팅 결과 텍스트 로그
- `{name}.pdf` — 분산 곡선 플롯

**stdout:** JSON 형식의 모든 잔기 피팅 값

---

## Config 파일 형식 (`.conf`)

```json
{
    "Project Name": "myproject",
    "comments": "",
    "datasets": [
        "exp-600MHz.dat",
        "exp-900MHz.dat"
    ],
    "exchange": "fast",
    "model": "Meiboom",
    "init": {
        "mode": "guess"
    },
    "residues": [
        {"name": "A10", "alpha": 1.44, "flag": "on"},
        {"name": "G15", "alpha": 0.83, "flag": "off"}
    ]
}
```

### 주요 필드

| 필드 | 값 | 설명 |
|------|-----|------|
| `exchange` | `"fast"` / `"slow"` | 교환 체제 |
| `model` | `"Meiboom"` / `"Matrix"` / `"London"` | 피팅 모델 |
| `init.mode` | `"guess"` / `"values"` | 초기값 설정 방식 |
| `residues[].flag` | `"on"` / `"off"` | 해당 잔기 분석 포함 여부 |
| `compare_aic` | `true` / `false` | global vs individual fit을 AIC로 비교 + 잔기별 우열 판정 (`sk_prepare` 기본 `true`, 키 없으면 `false`) |
| `jackknife` | `true` / `false` | leave-one-residue-out 잭나이프로 global fit 검증 (`sk_prepare` 기본 `true`, 키 없으면 `false`) |

### 초기값 직접 지정 (`init.mode = "values"`)

**Matrix / London 모델:**
```json
"init": {
    "mode": "values",
    "kex": 1000.0,
    "pB": 0.02,
    "csd": 2.5
}
```

**Meiboom 모델:**
```json
"init": {
    "mode": "values",
    "kex": 1000.0,
    "phi": 0.5
}
```

---

## 모델 선택 가이드

| 교환 체제 | 권장 모델 | 설명 |
|-----------|-----------|------|
| Fast | `Meiboom` | 빠른 교환 근사, 파라미터 적음 (φ, kex) |
| Fast | `Matrix` | 정확한 행렬 계산 (kAB, kBA, Δω) |
| Slow | `London` | 느린 교환 근사 |
| Slow | `Matrix` | 정확한 행렬 계산 |

> `alpha < 1` → slow exchange 가능성 높음  
> `alpha > 1` → fast exchange 가능성 높음

---

## Global vs Individual 피팅 비교 (AIC)

Config에 `"compare_aic": true`를 추가하면 두 모델을 비교한다.

- **Global fit** — 모든 잔기가 하나의 교환 과정(단일 kex, 또는 kAB/kBA)을 공유.
  잔기별 파라미터는 dd(csd 또는 phi)와 R2_0만.
- **Individual fit** — 각 잔기를 독립적으로 피팅. 잔기마다 자기 교환 속도를 가짐.

각 모델의 chi², 파라미터 수 k, 데이터 점 n으로 AIC를 계산한다:

```
AIC  = chi2 + 2k
AICc = AIC + 2k(k+1) / (n - k - 1)      # 소표본 보정
```

가중 잔차 (obs-calc)/sigma를 쓰므로 chi²가 -2 ln L 와 상수 차이. AICc가 낮은 모델이 선호되며,
Akaike weight로 상대적 지지도를 표시한다.

**출력:** `.log` 파일 끝과 stdout에 (1) 전체 데이터셋 비교 표, (2) 잔기별 individual 파라미터,
(3) **잔기별 우열 판정 표**.

```
Model               chi2      k      n          AIC         AICc
global            86.139     16    110      118.139      123.989
individual        83.084     20    110      123.084      132.523
---
delta AICc (individual - global):        8.534
Akaike weights (AICc): global=0.9862  individual=0.0138
Preferred model (whole dataset): global
```

### 잔기별 individual vs global 판정

전체 비교와 별도로, **잔기마다** 어느 쪽이 나은지 판정한다. 각 잔기의 global chi²
기여분(공유 kex 하에서의 그 잔기 잔차)과 individual chi²(자기 kex로 단독 피팅)을
잔기 단위 AICc로 비교한다. 파라미터 수 규약: individual은 자기 교환 파라미터(Meiboom 1개,
Matrix/London 2개)까지 비용을 지불하고, global은 공유 교환 속도를 무료로 빌린다 —
즉 "이 잔기가 자기만의 kex를 가질 값어치가 있는가?"를 묻는다.

```
 resId    chi2_glob   chi2_indiv        dAICc       z(kex)     better
   K1f       22.488       20.989        1.521        -1.43     global
   L2f       10.101        9.471        2.390         0.80     global
   ...
Residues preferring global: 5   preferring individual: 0
```

> `dAICc` = AICc(individual) − AICc(global): **> 0 → global**, **< 0 → individual** 선호.
> `z(kex)` = (kex_individual − kex_global) / σ(kex_individual): **|z| > 2**이면 그 잔기가
> 공유 kex와 통계적으로 불일치(잠재적 outlier).

---

## 잭나이프로 Global Fit 검증 (Jackknife)

Config에 `"jackknife": true`를 추가하면 global fit의 견고성을 검증한다.
활성 잔기를 하나씩 제외하고 global fit을 재실행(leave-one-residue-out)하여, 공유 교환
속도 kex가 특정 잔기에 얼마나 의존하는지 측정한다.

```
kex (full global fit):            1134.895
kex (jackknife mean):             1133.920
kex (bias-corrected):             1138.793
jackknife bias:                     -3.899
jackknife std error (kex):          67.180
relative std error:                  5.92%
---
 resId      kex(-res)      delta_kex        |z|
   K1f       1166.713        -31.818       0.47
   N4f       1179.261        -44.366       0.66
   ...
(* = influential: dropping this residue shifts kex by > 2 jackknife SE)
```

계산식 (표준 leave-one-out 잭나이프, N = 활성 잔기 수):

```
delta_kex_i = kex_full - kex_(-i)                        # 잔기 i 제외 시 kex 변화
bias        = (N-1) * (mean(kex_(-i)) - kex_full)
SE          = sqrt( (N-1)/N * sum_i (kex_(-i) - mean)^2 )
```

> **relative SE가 작고** `*`로 표시된 잔기가 없으면 global fit은 견고하다 —
> 어느 한 잔기도 공유 kex를 지배하지 않는다는 뜻.
> 특정 잔기가 `*`로 표시되면 그 잔기를 빼면 kex가 2 SE 이상 흔들린다 → outlier 후보.

```bash
# 실행 (config에 compare_aic / jackknife 설정 후 — sk_prepare는 기본 활성화)
python sk_run.py project.conf
```

**stdout JSON:** 잔기 결과(마지막 `#####` 줄)에 더해, 활성화 시
`##### model_comparison`과 `##### jackknife` 라벨의 JSON 블록이 앞쪽에 출력된다.
잔기 목록 줄은 마지막에 그대로 유지되어 기존 파서와 호환된다.

---

## 데이터 파일 형식 (`.dat`)

```
 60.12          ← 자기장 강도 (MHz, 1H)
0.040000        ← TCP (s)
#nu_cpmg(Hz)        R2(1/s)      Esd(R2)
# A10           ← 잔기 이름 (# 주석)
  50.000   22.415    0.448
 100.000   21.893    0.438
 200.000   20.912    0.418
 ...
# G15
  50.000   18.234    0.365
 ...
```

---

## 전체 검증 테스트

```bash
bash test_migration.sh        # 전체 파이프라인 검증
bash test/run_tests.sh        # fast/slow exchange 테스트 케이스
```

---

## 전형적인 워크플로우

```
실험 데이터 (.dat 파일)
        ↓
sk_prepare.py → config 초안 생성 (.conf)
        ↓
config 편집 (exchange, model, flag, 초기값)
        ↓
sk_run.py → 피팅 실행
        ↓
결과 확인 (.log, .pdf)
        ↓
필요시 config 수정 후 재실행
```
