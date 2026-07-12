#!/usr/bin/env python3
"""
02_fetch_tcga.py — Retrieve TCGA-LUAD/LUSC data via the cBioPortal API and
assemble the TCGA analysis dataset (genomic + clinical) plus the signature-gene
expression matrix used for transcriptomic subtype transfer.

Outputs (under DATA/tcga/):
  tcga_analysis_dataset.csv      (drivers, TMB, FGA, aneuploidy, stage, OS, DFS)
  tcga_signature_expr_matrix.csv (built later by 12_tcga_build_signature -> here
                                   we fetch generic driver/clinical; signature
                                   expression is fetched in step 14's helper)

Note: study clinical files (stage, OS, DFS) are read from local cBioPortal
clinical_patient files if present; otherwise fetch via API. Driver mutations,
CNA, TMB, FGA are fetched via API. Requires network access to cbioportal.org.
"""
import os, json, subprocess, pickle
import pandas as pd, numpy as np
from config import TCGADIR, CBIO_API

STUDIES = ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"]
DRIVERS = ["ARID1A","ATM","BRAF","CDKN2A","EGFR","KEAP1","KRAS","MET","NF1","NFE2L2",
           "PIK3CA","PTEN","RB1","RBM10","SMARCA4","STK11","TP53"]

def curl_json(url, payload=None):
    args = ["curl","-s","--max-time","120","-H","accept: application/json"]
    if payload is not None:
        pf = os.path.join(TCGADIR,"_payload.json"); json.dump(payload, open(pf,"w"))
        args += ["-H","Content-Type: application/json","-d","@"+pf]
    args.append(url)
    r = subprocess.run(args, capture_output=True, text=True, timeout=150)
    try: return json.loads(r.stdout)
    except Exception as e:
        print("  parse error:", e, r.stdout[:200]); return None

def entrez_ids(genes):
    r = subprocess.run(["curl","-s","--max-time","60","-H","Content-Type: application/json",
                        "-H","accept: application/json","-d",json.dumps(genes),
                        f"{CBIO_API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL"],
                       capture_output=True, text=True, timeout=90)
    return {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in json.loads(r.stdout)}

def sample_ids(study, profile_suffix="_all"):
    r = subprocess.run(["curl","-s","--max-time","60","-H","accept: application/json",
                        f"{CBIO_API}/sample-lists/{study}{profile_suffix}/sample-ids"],
                       capture_output=True, text=True, timeout=90)
    return json.loads(r.stdout)

def fetch_mut_cna(study, eids):
    muts = curl_json(f"{CBIO_API}/molecular-profiles/{study}_mutations/mutations/fetch",
                     {"entrezGeneIds":eids,"sampleListId":f"{study}_all"})
    cnas = curl_json(f"{CBIO_API}/molecular-profiles/{study}_gistic/discrete-copy-number/fetch"
                     f"?discreteCopyNumberEventType=ALL",
                     {"entrezGeneIds":eids,"sampleListId":f"{study}_all"})
    return muts, cnas

def fetch_clin_genomic(study):
    sids = sample_ids(study)
    ident = [{"entityId":s,"studyId":study} for s in sids]
    return curl_json(f"{CBIO_API}/clinical-data/fetch?clinicalDataType=SAMPLE",
                     {"attributeIds":["TMB_NONSYNONYMOUS","FRACTION_GENOME_ALTERED",
                                      "ANEUPLOIDY_SCORE","MUTATION_COUNT"],
                      "identifiers":ident})

def to_wide(data):
    df = pd.DataFrame([{"patient":d["patientId"],"attr":d["clinicalAttributeId"],"val":d["value"]}
                       for d in data])
    return df.pivot_table(index="patient",columns="attr",values="val",aggfunc="first").reset_index()

def main():
    emap = entrez_ids(DRIVERS); e2s = {str(v):k for k,v in emap.items()}
    eids = list(emap.values())
    all_rows=[]
    for study in STUDIES:
        print(f"[TCGA] {study} ...")
        muts, cnas = fetch_mut_cna(study, eids)
        gen = fetch_clin_genomic(study)
        # mutation binary
        mdf = pd.DataFrame([{"patient":m["patientId"],"gene":e2s.get(str(m["entrezGeneId"]))}
                            for m in muts if str(m["entrezGeneId"]) in e2s])
        mbin = pd.crosstab(mdf.patient, mdf.gene).clip(upper=1)
        for g in DRIVERS:
            if g not in mbin.columns: mbin[g]=0
        mbin = mbin[DRIVERS].add_prefix("mut_").reset_index()
        gwide = to_wide(gen)
        d = gwide.merge(mbin, on="patient", how="left")
        for c in [x for x in d.columns if x.startswith("mut_")]:
            d[c]=d[c].fillna(0).astype(int)
        # clinical stage/OS/DFS from local cBioPortal patient file if present, else API
        clin_path = os.path.join(TCGADIR, f"{study}_clinical_patient.txt")
        if os.path.exists(clin_path):
            cp = pd.read_csv(clin_path, sep="\t", comment="#", low_memory=False)
            cp = cp.rename(columns={"PATIENT_ID":"patient"})
            keep = ["patient","AJCC_PATHOLOGIC_TUMOR_STAGE","OS_MONTHS","OS_STATUS",
                    "DFS_MONTHS","DFS_STATUS"]
            cp = cp[[c for c in keep if c in cp.columns]]
            d = d.merge(cp, on="patient", how="left")
        d["cohort"] = study.split("_")[0].upper()
        d["histology"] = d["cohort"]
        all_rows.append(d)
    tcga = pd.concat(all_rows, ignore_index=True)
    for c in ["TMB_NONSYNONYMOUS","FRACTION_GENOME_ALTERED","ANEUPLOIDY_SCORE",
              "MUTATION_COUNT","OS_MONTHS","DFS_MONTHS"]:
        if c in tcga.columns: tcga[c]=pd.to_numeric(tcga[c],errors="coerce")
    if "AJCC_PATHOLOGIC_TUMOR_STAGE" in tcga.columns:
        tcga["stage_simple"]=tcga.AJCC_PATHOLOGIC_TUMOR_STAGE.str.extract(r"STAGE\s+(IV|III|II|I)")
    tcga.to_csv(os.path.join(TCGADIR,"tcga_analysis_dataset.csv"), index=False)
    print(f"[TCGA] analysis dataset: {len(tcga)} patients -> {TCGADIR}")
    # clean temp
    tmp=os.path.join(TCGADIR,"_payload.json")
    if os.path.exists(tmp): os.remove(tmp)
    print("=== 02_fetch_tcga.py complete ===")
    print("NOTE: signature-gene expression matrix (tcga_signature_expr_matrix.csv) is "
          "produced by the R step 12/14 helper after the ALCHEMIST signature is defined.")

if __name__ == "__main__":
    main()
