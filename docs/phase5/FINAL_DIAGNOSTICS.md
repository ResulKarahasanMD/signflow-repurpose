# Phase 5 — Final diagnostics (§A1–A4)

Four diagnostics that sharpen the central Phase 5 claim. §A1 is the decisive one.
All figures carry the chance baseline; all CIs are `n_boot=1000` bootstrap intervals.
Reproduced by:

- `scripts/run_phase5_a1_control_graphs.py` (A1, A2 Spearman)
- `scripts/run_phase5_a3_top50_profile.py` (A3)

Chance AUPRC on the benchmark scored pool = 29 / 2833 = **0.0102**.

---

## A1 — BENCHMARK control-graph AUPRC (decisive)

Control-graph scoring previously ran for TARGET only. This reports it for BENCHMARK,
holding the module and disease direction fixed so **only the network changes**. Each
control preserves structure that the real graph carries: `rewired` preserves the exact
degree sequence, `edge_shuffled` preserves the edge multiset, `sign_shuffled` preserves
topology and permutes edge signs/relations.

| Arm | AUPRC [95% CI] (chance 0.0102) | AUROC [95% CI] | Spearman(scores) vs real |
|---|---|---|---|
| signed, real graph | **0.035 [0.017, 0.069]** | **0.665 [0.537, 0.785]** | — |
| signed, **rewired** | 0.033 [0.016, 0.063] | 0.595 [0.446, 0.749] | **0.799** |
| signed, **edge_shuffled** | 0.027 [0.014, 0.052] | 0.641 [0.517, 0.765] | **0.764** |
| signed, **sign_shuffled** | 0.027 [0.015, 0.052] | 0.661 [0.543, 0.779] | **0.740** |

Detectable-threshold check: every control AUPRC CI **overlaps** the real-graph CI; the
observed real−control AUPRC gaps (0.0015 rewired, 0.0076 edge_shuffled, 0.0079
sign_shuffled) are an order of magnitude **below** the minimum detectable difference
(≈ 0.045). Drug scores correlate ρ = 0.74–0.80 with the real graph even after its
structure is destroyed.

**Verdict (pre-committed before running):** the real signed graph does **not** materially
exceed the structure-preserving controls, so the control-graph score **carries no network
information**. The modest discrimination that exists comes from the **drug target sets and
the disease signature**, not from propagation over the molecular network. This is the
stronger, more useful finding than "weak signal, no utility": the network component of the
method is inert on this benchmark. (`verdict.network_information_present = false`.)

---

## A2 — Invariance / instability panel

Three previously-scattered numbers plus the A1 Spearman quantification tell one
mechanistic story.

| Perturbation | Statistic | Value | Implication |
|---|---|---:|---|
| Substrate changed: causal → causal+STRING (~10× edges), TARGET track | Top-20 Jaccard | **0.905** | ranking nearly network-invariant |
| Edge set reduced 34.3% (unsigned_matched vs unsigned_full), BENCHMARK | ΔAUPRC | **≈ 0.0017** (0.0351 vs 0.0335) | connectivity loss costs ~nothing |
| Signature bootstrap H7 (TARGET) | mean top-20 Jaccard | **0.051** (thr 0.20) | ranking highly signature-unstable |
| Real vs **rewired** graph (BENCHMARK) | Spearman(scores) | **0.799** | scores network-invariant (A1) |
| Real vs **edge_shuffled** (BENCHMARK) | Spearman(scores) | **0.764** | scores network-invariant (A1) |
| Real vs **sign_shuffled** (BENCHMARK) | Spearman(scores) | **0.740** | sign carries no information |

For completeness, the TARGET-track control-graph Spearman (different seeds, injection
dominance) is lower and sign-sensitive: rewired ρ ≈ 0.33, edge_shuffled ρ ≈ 0.18,
sign_shuffled ρ ≈ −0.27.

**Statement:** the drug ranking is dominated by the **disease signature** and the
**drug-injection term**, and is close to invariant to the molecular network — invariant to
substrate choice (0.905), to a 34.3% connectivity loss (ΔAUPRC ≈ 0.002), and to graph
rewiring/edge-shuffling (score Spearman 0.74–0.80), while being highly **unstable** to the
disease signature (H7 Jaccard 0.051).

---

## A3 — Top-50 non-gold characterisation

Gold drugs are absent from the top 50 yet enriched 4.9× at k=100 and 5.5× at k=250.
Profile of what outranks them (signed-flow score; degrees on the real signed benchmark
graph; ATC from ChEMBL):

| Feature (median / mean) | Top-50 non-gold | Gold set | Full scored pool |
|---|---:|---:|---:|
| mapped targets (total) | 3 / 6.4 | 4 / 14.3 | 3 / 8.3 |
| targets in substrate | 2 / 5.2 | 4 / 11.3 | 3 / 5.9 |
| mean target degree | 13.5 / 22.0 | 19.3 / 25.7 | 10.7 / 20.4 |
| max target degree | 20 / 57.8 | 58 / 77.2 | 17 / 55.1 |

ATC level-1 (therapeutic area) of the top-50 non-gold: **Antineoplastic &
immunomodulating 17**, **Musculo-skeletal 11**, Nervous 5, Sensory 4, Blood 2,
Alimentary 2, Respiratory 1. (Gold set: Antineoplastic & immunomodulating 11, Alimentary
8, Dermatologicals 5, Systemic hormonal 5, …)

**Reading.** The top-50 non-gold are **not more promiscuous than gold** — gold actually
carries more targets (median 4 vs 3) and higher-degree targets (median max-degree 58 vs
20). They are, however, modestly **hub-enriched relative to the pool** (mean target degree
13.5 vs 10.7 median) and are dominated by **immunomodulatory / anti-inflammatory** ATC
classes. What outranks gold is therefore **signature alignment via the injection term**:
drugs whose targets sit on inflammatory hubs that the injection propagates into the
disease module score highly regardless of UC-approval status. This is consistent with the
A1 network-invariance finding and explains the shape of the enrichment curve (broad
signature-aligned enrichment at k=100/250, none at the very top where UC-specific gold
would need target-set specificity the score does not use).

---

## A4 — Applicability limit (Limitations text)

**Pentosan polysulfate sodium — the applicability limit of target-based network methods.**
Pentosan polysulfate sodium (CHEMBL1201516) is the **sole Open Targets APPROVAL-stage
agent for the target disease `MONDO_0018301`**. It is `not_in_scored_pool`: it has **no
protein-target mechanism** to inject into the signed-flow substrate. Stated plainly: **the
only approved drug for the target disease is unscoreable in principle by target-based
network methods.** This generalises to a whole class of agents — **polysaccharides and
other non-protein-target mechanisms, and physical/barrier agents** (pentosan polysulfate
is a semi-synthetic polysaccharide acting as a urothelial glycosaminoglycan-layer
replenisher). It is a limitation of the **method family**, not of this implementation: any
method that represents drugs solely by protein-target injection cannot rank such agents.

**Coverage caveat (attach to every enrichment figure).** Substrate coverage correlates
with approval status: gold coverage is **100%** (29/29) versus non-gold **88.5%**
(2804/3168). This inflates gold prevalence from **0.686%** on the acquired universe to
**1.024%** on the scored pool (**+49%** relative). Every enrichment/precision figure is
computed on the scored pool and therefore sits on the inflated 1.02% prevalence; the
chance baseline (0.0102) already reflects it, but the differential coverage means gold is
over-represented among scoreable drugs before any propagation occurs.

Soft positives at their ranks (TARGET, primary seeds × primary substrate; reported where
they land, never forced): cyclosporine rank 20, amitriptyline rank 502, mirabegron rank
2432; pentosan polysulfate `not_in_scored_pool` (above).
