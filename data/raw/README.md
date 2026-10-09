# Raw source files

`code/01_prepare_data.py` downloads the public OWID Energy dataset and codebook
to this directory and verifies them against the frozen SHA-256 checksums before
analysis. Raw CSV files are intentionally excluded from Git; the derived
analytical panel and all reported outputs are versioned in the repository.
