# Preprint build and updates

The Markdown manuscript is the editable source. The PDF uses the shared MiAnemone/MiUrchin layout. Scientific inputs remain in the repository root.

Install `requirements-preprint.txt`, pandoc and Chrome; run `python publications/plot_paper.py`, `python publications/build_preprint.py`, `python publications/check_preprint.py`, `python scripts/check_publication.py` and `python scripts/preprint_integrity.py`. GitHub Actions rebuilds, checks and publishes a versioned PDF and primer FASTA after publication inputs change.
