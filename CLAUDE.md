# CLAUDE.md — bio-research (SignFlow-Repurpose)

## Önce bunu oku
Bu depoda **`README.md` ve `docs/phase0/` altındaki ön kayıtlı plan bağlayıcıdır.** Faz durumu, veri kapısı, ön kayıt (preregistration) ve değişiklik kaydı orada; burada tekrarlamıyorum. Plan dışına çıkan analiz yapma — `docs/phase0/ANALYSIS_PLAN.md` ön kayıtlıdır, sapma `ANALYSIS_PLAN_AMENDMENTS.md`'ye yazılır.

## Ne bu
IC/BPS (asıl uygulama) için kontrollü işaretli-ağ önceliklendirmesi + çok kanıtlı üçgenleme; ülseratif kolit metot kıyas kolu. Paket adı `signflow-repurpose`, CLI girişi `signflow` (`src/signflow/cli.py`).

Python ≥3.12, uv (`uv.lock`, `.venv`). pydeseq2 / decoupler / gseapy / torch / torch-geometric / rpy2.

## Komutlar
```bash
uv sync
uv run pytest                      # pythonpath = ["src"]
uv run python scripts/run_smoke.py
uv run python scripts/run_phase5.py
```

`scripts/run_phase*.py` faz koşucularıdır; sırayı atlama, her faz kendi kapısından geçer.

## Dikkat
- **Hipotez üretimi, tedavi önerisi değil.** Skorları PK/PD, nedensellik veya klinik etkinlik gibi sunma.
- Etiketler Open Targets allowlist'inden ve eşik ≥0.30 ile dondurulmuştur (`docs/phase0/LABEL_PROVENANCE.md`). DE sonuçlarını etiket olarak kullanma — döngüsel olur.
- Etiket kapısı: IC n=38 → birincil analiz **gözetimsiz**. Bunu gözetimliye çevirme.
- `provenance/`, `runs/`, `logs/` üretilmiş kayıtlar — elle düzenleme, silmeden önce sor.
- Veri depoda yok (`data/` yalnız `.gitkeep`); faz koşucuları indirir. Yeni büyük çıktıları `runs/` altında tarihli dizine yaz.
- Ana dalda doğrudan çalışma, izinsiz commit/push yok.
