#!/usr/bin/env python3
"""
02b_fetch_tcga_signature_expr.py — Fetch TCGA-LUAD/LUSC mRNA expression for the
subtype signature genes and build the matrix used by the transcriptomic-axis
transfer (R step 14).

Run AFTER R step 12 (12_tcga_build_signature.R), which writes
results/signature_genes.csv (the 154 signature-gene list).

Output: data/tcga/tcga_signature_expr_matrix.csv  (genes x TCGA patients, RSEM)

Requires network access to cbioportal.org.
"""
import os, json, subprocess
import pandas as pd
from config import TCGADIR, RESULTS, CBIO_API

STUDIES = ["luad_tcga_pan_can_atlas_2018", "lusc_tcga_pan_can_atlas_2018"]

def entrez_ids(genes):
    r = subprocess.run(["curl","-s","--max-time","60","-H","Content-Type: application/json",
                        "-H","accept: application/json","-d",json.dumps(genes),
                        f"{CBIO_API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL"],
                       capture_output=True, text=True, timeout=90)
    return {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in json.loads(r.stdout)}

def fetch_expr(study, eids):
    prof = f"{study}_rna_seq_v2_mrna"
    payload = {"entrezGeneIds": eids, "sampleListId": prof}
    pf = os.path.join(TCGADIR,"_ep.json"); json.dump(payload, open(pf,"w"))
    r = subprocess.run(["curl","-s","--max-time","120","-H","Content-Type: application/json",
                        "-H","accept: application/json","-d","@"+pf,
                        f"{CBIO_API}/molecular-profiles/{prof}/molecular-data/fetch"],
                       capture_output=True, text=True, timeout=150)
    try: return json.loads(r.stdout)
    except Exception as e:
        print("  parse error:", e, r.stdout[:200]); return None

def build(data, e2s):
    rows = [{"patient":x["patientId"],"entrez":str(x["entrezGeneId"]),"val":x["value"]} for x in data]
    df = pd.DataFrame(rows); df["gene"]=df.entrez.map(e2s); df=df.dropna(subset=["gene"])
    return df.pivot_table(index="gene",columns="patient",values="val",aggfunc="mean")

def main():
    sig = pd.read_csv(os.path.join(RESULTS,"signature_genes.csv"))
    genes = sig["gene"].tolist()
    print(f"[sig-expr] {len(genes)} signature genes")
    emap = entrez_ids(genes); e2s = {str(v):k for k,v in emap.items()}
    eids = list(emap.values())
    mats=[]
    for study in STUDIES:
        data = fetch_expr(study, eids)
        m = build(data, e2s)
        print(f"  {study}: {m.shape[0]} genes x {m.shape[1]} samples")
        mats.append(m)
    common = mats[0].index.intersection(mats[1].index)
    out = pd.concat([mats[0].loc[common], mats[1].loc[common]], axis=1)
    out.index.name="gene"
    out.to_csv(os.path.join(TCGADIR,"tcga_signature_expr_matrix.csv"))
    tmp=os.path.join(TCGADIR,"_ep.json")
    if os.path.exists(tmp): os.remove(tmp)
    print(f"[sig-expr] saved tcga_signature_expr_matrix.csv: {out.shape[0]} genes x {out.shape[1]} samples")

if __name__ == "__main__":
    main()
