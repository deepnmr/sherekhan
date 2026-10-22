# ShereKhan 1.3.0 레퍼런스

> 🌐 **한국어** · [English](USAGE.en.md)

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

옵션: `--seed N`, `--outdir DIR` ([CLI 레퍼런스](#cli-레퍼런스) 참고).

**파라미터 수정 (스크립트 내부):**

| 변수 | 기본값 | 설명 |
|------|--------|------|
| `fields` | `[60.12, 90.23]` | 관측 핵의 Larmor 주파수 (MHz) |
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
python sk_run.py <run.conf>
```

**예시:**
```bash
python sk_run.py project.conf
```

**출력 파일 (`Project Name` 기준):**
- `{name}.log` — 피팅 결과 텍스트 로그
- `{name}.pdf` — 분산 곡선 플롯

**stdout:** 배너, 로그 텍스트, 라벨된 JSON 줄. 마지막 줄이 잔기별 JSON ([stdout 구성](#sk-run-stdout-구성) 참고)

---

## CLI 레퍼런스

```
sk-run CONFIG [--quiet] [--version] [-h]
sk-prepare DAT [DAT ...] [--version] [-h]
sk-synth [--seed N] [--outdir DIR] [--version] [-h]
```

- `--quiet` (sk-run): 피팅 중 chi2 진행 줄을 출력하지 않는다.
- `--seed N` (sk-synth): 재현 가능한 노이즈용 seed (기본값: seed 없음).
- `--outdir DIR` (sk-synth): `synth-60.dat`/`synth-90.dat`를 쓸 디렉터리 (없으면 생성, 기본값: 현재 디렉터리).
- `--version`: `<prog> 1.3.0`을 출력하고 종료. `-h`/`--help`는 사용법 출력. 둘 다 파일을 쓰지 않는다.

`python sk_run.py`, `python sk_prepare.py`, `python sk_createSyntheticDataset.py` 형태도 같은 옵션을 받는다.

### 종료 코드

| 코드 | 의미 |
|------|------|
| `0` | 성공 |
| `1` | 입력 또는 실행 오류: stderr에 `<prog>: error: <message>` 한 줄 (파일·키·줄을 명시). traceback 없음. |
| `2` | 사용법 오류 (잘못되거나 빠진 인자) |

### sk-run stdout 구성

sk-run stdout은 하나의 JSON 문서가 아니다. 항상 이 순서로 나온다:

1. 배너 (버전, 인용 정보)
2. chi2 진행 줄 (`--quiet`이면 생략)
3. 로그 텍스트 (`{name}.log`와 같은 내용)
4. `Calculations finished successfully.`
5. 선택: `##### model_comparison` 다음 JSON 한 줄 (`compare_aic` 켜짐)
6. 선택: `##### jackknife` 다음 JSON 한 줄 (`jackknife` 켜짐)
7. `#####` 다음 잔기별 결과 JSON 한 줄. 항상 마지막 줄이다.

sk-prepare는 stdout에 config JSON만 쓴다. 경고는 stderr로 가므로 `sk-prepare a.dat > run.conf`는 항상 유효한 JSON이다.

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
| `exchange` | `"fast"` / `"slow"` / `"undefined"` | 교환 체제. `"undefined"`는 `Matrix` 또는 `London`으로 실행되며 `"slow"`처럼 보고된다. sk-prepare는 이때 `Matrix`를 제안. 애매한 다중 필드(fast/slow 동률) → 모든 잔기 자동 `flag: "off"`. 단일 필드 → 마찬가지 `"undefined"`/`Matrix`지만 모든 잔기 `flag: "on"` (alpha는 필드 ≥ 2 필요) |
| `model` | `"Meiboom"` / `"Matrix"` / `"London"` | 피팅 모델 |
| `init.mode` | `"guess"` / `"values"` | 초기값 설정 방식 |
| `residues[].flag` | `"on"` / `"off"` | 해당 잔기 분석 포함 여부 |
| `compare_aic` | `true` / `false` | global vs individual fit을 AIC로 비교 + 잔기별 우열 판정 (`sk_prepare` 기본 `true`, 키 없으면 `false`) |
| `jackknife` | `true` / `false` | leave-one-residue-out 잭나이프로 global fit 검증 (`sk_prepare` 기본 `true`, 키 없으면 `false`) |

### Config 규칙 (sk-run이 검증)

- `exchange`는 `fast`, `slow`, `undefined` 중 하나 (대소문자 구분). `undefined`는 `Matrix` 또는 `London`과 함께 쓸 수 있고 `slow`처럼 보고된다. Meiboom은 `fast` 필요.
- `model`은 `Matrix`, `London`, `Meiboom` 중 하나. `init.mode`는 `guess` 또는 `values`.
- `datasets`의 상대 경로가 현재 디렉터리에서 없으면 config 파일 디렉터리 기준으로 찾는다.
- `init.mode = "values"`이면 `kex`는 0보다 큰 유한한 수여야 한다. Matrix/London은 0 < pB < 1인 유한한 `pB`와 유한한 `csd`가, Meiboom은 유한한 `phi`가 필요하다.
- `residues`에 없는 잔기는 활성 상태로 남는다.
- 규칙 위반은 키를 명시한 `sk-run: error: ...` 한 줄과 함께 종료 코드 1.

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

## 모듈 구성

1.3.0부터 피팅 코드는 `cpmg/` 아래 네 모듈로 나뉜다.

- `cpmg/models.py`: R2 공식 (Matrix, London, Meiboom).
- `cpmg/fitting.py`: 피팅 엔진.
- `cpmg/analysis.py`: 모델 비교, 잭나이프, 파생 불확도.
- `cpmg/report.py`: log, PDF, JSON 출력. matplotlib은 그림을 그릴 때만 지연 import한다.

공개 진입점은 여전히 `cpmg.model_2state.CPMG_model`이라 기존 스크립트는 그대로 동작한다.

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
 60.12          ← field (MHz, 관측 핵)
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


- 1번째 줄 `field`: 관측 핵의 Larmor 주파수 (MHz). 예: 600 MHz 분광기의 15N은 60.8 MHz. 2번째 줄: tcp (s).
- 데이터 줄은 숫자 정확히 3개, `Esd(R2)`는 0보다 커야 한다. 잔기 라벨 줄은 `#`으로 시작.
- 규칙 위반 파일은 `<prog>: error: <file>:<line>: ...`과 종료 코드 1로 멈춘다.

---

## 전체 검증 테스트

```bash
pip install -e ".[test]"   # needs pip >= 21.3
pytest
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
