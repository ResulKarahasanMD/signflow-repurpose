# Pre-registration timestamp (Phase 0.6)

**Intent:** Externally verifiable freeze of the analysis plan **before Phase 2 data download**.

## Local git objects (both are intentional)

| Object | Hash / name | What it covers |
|---|---|---|
| **Plan-content freeze commit** | `2b1a39af2a2cb6ee9595bc1b8f251d4210a389ab` | Phase 0.5 docs bundle as first committed plan content |
| **Annotated tag** | `prereg-phase0.5` → `git rev-parse prereg-phase0.5^{}` | Tag tip may differ from the plan-content freeze commit; **this is expected** — the tag marks the pre-registration *record*, while the freeze commit marks the *plan text snapshot* |

```bash
git rev-parse 2b1a39af2a2cb6ee9595bc1b8f251d4210a389ab
git rev-parse prereg-phase0.5^{}
git show prereg-phase0.5
```

## External receipt (A6.3) — required before Phase 2

| Field | Value |
|---|---|
| **Public GitHub Gist** | https://gist.github.com/ResulKarahasanMD/f3d6af1db8b405f24b0026e1e13315ce |
| **Bundle SHA256** (Phase 0.5 core docs at receipt time) | `68d4eeaedffe2ea76806a08ac6118e709ad07992817371045ce2df7adbc97b63` |
| **Date (UTC)** | 2026-07-23 |

A follow-up gist will record the Phase 0.6 plan hash after this amendment commit.

**Manuscript statement (draft):**  
*The analysis plan was pre-registered prior to Phase 2 data acquisition via a public GitHub Gist receipt (https://gist.github.com/ResulKarahasanMD/f3d6af1db8b405f24b0026e1e13315ce) recording SHA-256 `68d4eeaedffe2ea76806a08ac6118e709ad07992817371045ce2df7adbc97b63`, with local git tag `prereg-phase0.5` and plan-content commit `2b1a39af…`. Phase 0.6 amendments were logged before Phase 1/2 execution.*
