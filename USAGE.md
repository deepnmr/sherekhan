# ShereKhan 1.2.0 레퍼런스

NMR CPMG 분산 곡선 분석 도구 — 스크립트·config·데이터 형식 레퍼런스.  
저자: Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)  
원본 코드: Adam Mazur, Bjoern Hammesfahr, Christian Griesinger, Donghan Lee, Martin Kollmar (Max-Planck-Institute for Biophysical Chemistry, 2012)  
Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)

> **판정·검증 절차** (individual vs global fit, AIC 비교, 잭나이프, 이질성 진단)는
> 단계별 매뉴얼 [`MANUAL.ko.md`](MANUAL.ko.md) · [English](MANUAL.en.md) 참고.
> 이 문서는 스크립트·config·데이터 **형식 레퍼런스**만 다룬다.

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
| `exchange` | `"fast"` / `"slow"` / `"undefined"` | 교환 체제. `"undefined"`는 항상 모델 `Matrix`와 짝. 애매한 다중 필드(fast/slow 동률) → 모든 잔기 자동 `flag: "off"`. 단일 필드 → 마찬가지 `"undefined"`/`Matrix`지만 모든 잔기 `flag: "on"` (alpha는 필드 ≥ 2 필요) |
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

## 판정·검증 (AIC 비교 · 잭나이프)

`compare_aic: true` / `jackknife: true`일 때 global vs individual 판정과 잭나이프 검증
결과가 `.log`·stdout에 출력된다. **읽는 법과 이질성 진단 절차는 단계별 매뉴얼**
[`MANUAL.ko.md`](MANUAL.ko.md) §5–8 참고.

stdout JSON 연동 (프로그램 처리) — 라벨된 블록:

```
##### model_comparison
{ ... compareModelsAIC 결과 (per_residue 포함) ... }
##### jackknife
{ ... jackknifeGlobal 결과 ... }
#####
[ ... 잔기별 최종 fit 값 ... ]
```

- `##### model_comparison` / `##### jackknife` 블록은 strict JSON (비유한값은 `null`).
- **마지막 `#####` 다음 줄이 잔기 결과 리스트** (기존 파서 호환, `split('#####')[-1]`).
- `compare_aic`/`jackknife`가 꺼져 있으면 해당 블록은 출력되지 않는다.

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
