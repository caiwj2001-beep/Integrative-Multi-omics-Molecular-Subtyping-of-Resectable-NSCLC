#!/usr/bin/env python3
"""config.py — shared configuration for the Python data-preparation steps (portable).

Edit PROJECT_ROOT (or set env NSCLC_PROJECT_ROOT). Raw-cohort directories
ALCHEMIST_DIR / MSK_DIR point to the downloaded source data (see README section 1).
"""
import os

PROJECT_ROOT = os.environ.get(
    "NSCLC_PROJECT_ROOT",
    r"C:\Users\Administrator\科研文件\NSCLC多组学分子分型")

# Raw source-data directories (downloaded from GDC / cBioPortal; not redistributed)
ALCHEMIST_DIR = os.environ.get(
    "NSCLC_ALCHEMIST_DIR",
    r"C:\Users\Administrator\科研文件\ALCHEMIST肺癌试验")
MSK_DIR = os.environ.get(
    "NSCLC_MSK_DIR",
    r"C:\Users\Administrator\科研文件\smarca4\data\supplementary\msk_impact_50k_2026")

# Derived project paths
DATA    = os.path.join(PROJECT_ROOT, "data")
MATDIR  = os.path.join(DATA, "matrices")
TCGADIR = os.path.join(DATA, "tcga")
MSKOUT  = os.path.join(DATA, "msk")
RESULTS = os.path.join(PROJECT_ROOT, "results")
WORK    = os.path.join(ALCHEMIST_DIR, "work")

for d in (DATA, MATDIR, TCGADIR, MSKOUT, RESULTS):
    os.makedirs(d, exist_ok=True)

GDC_FILES_API = "https://api.gdc.cancer.gov/files"
CBIO_API = "https://www.cbioportal.org/api"

if __name__ == "__main__":
    print("PROJECT_ROOT:", PROJECT_ROOT)
    print("ALCHEMIST_DIR:", ALCHEMIST_DIR)
    print("MSK_DIR:", MSK_DIR)
