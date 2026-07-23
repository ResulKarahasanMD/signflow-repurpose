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
| **Public GitHub Gist (Phase 0.5)** | https://gist.github.com/ResulKarahasanMD/f3d6af1db8b405f24b0026e1e13315ce |
| **Public GitHub Gist (Phase 0.6)** | https://gist.github.com/ResulKarahasanMD/958b8c80430732caff2bada19b5b2c6d |
| **Bundle SHA256** (Phase 0.5 core docs at first receipt) | `68d4eeaedffe2ea76806a08ac6118e709ad07992817371045ce2df7adbc97b63` |
| **Bundle SHA256** (all `docs/phase0/*.md` at Phase 0.6 commit) | `99789eb553084d57489d7bf87b437db80b83eb1843e764d342166723b83756a6` |
| **Phase 0.6 commit** | `d2e4bef71b7f924db761bbfb9c10603df82b17a2` |
| **Date (UTC)** | 2026-07-23 |

**Phase 2 download gate:** **CLEARED** — external receipts above are HTTP 200 as of 2026-07-23 (re-verified before Phase 2). Public GitHub Gists are the external timestamp (immutable content at gist URLs); local tag alone is not relied upon. OSF/Zenodo DOI remains optional enhancement.
