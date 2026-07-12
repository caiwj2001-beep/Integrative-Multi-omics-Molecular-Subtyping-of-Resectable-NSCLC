#!/usr/bin/env python3
"""
01_prepare_data.py — Build analysis matrices from raw ALCHEMIST + MSK-IMPACT data.

Steps:
  A. Map ALCHEMIST GDC assay files (named by aliquot UUID) to patient barcodes
     via the GDC Files API (queried on file_name).
  B. Parse ALCHEMIST clinical JSON -> master clinical table (treatment arm,
     histology, sex, age).
  C. Build per-patient matrices: expression (TPM), mutation (binary drivers +
     merged MAF + TMB), gene-level CNA, miRNA (RPM); tag multi-omics availability.
  D. Build MSK-IMPACT NSCLC cohort: clinical (OS), mutation binary, gene/arm CNA.

Outputs (under DATA/):
  matrices/expr_tpm.parquet, expr_counts.parquet
  matrices/mutation_binary.tsv.gz, alchemist_merged.maf.gz, tmb.csv
  matrices/cna_gene.tsv.gz, mirna_rpm.tsv.gz
  alchemist_master_clinical.csv
  msk/msk_nsclc_clinical.csv, msk/msk_mutation_binary.csv,
  msk/msk_gene_cna.csv, msk/msk_armlevel_cna.csv

Requires: pandas, numpy, pyarrow; network access to api.gdc.cancer.gov.
"""
import os, json, glob, gzip, time, subprocess
from collections import defaultdict, Counter
import pandas as pd, numpy as np
from config import ALCHEMIST_DIR, MSK_DIR, DATA, MATDIR, MSKOUT, WORK, GDC_FILES_API

os.makedirs(WORK, exist_ok=True)
NONSYN = {"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
          "In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site",
          "Nonstop_Mutation","Splice_Region"}
DRIVERS = ["TP53","EGFR","KRAS","STK11","KEAP1","SMARCA4","NF1","BRAF","MET","RBM10",
           "PIK3CA","RB1","CDKN2A","ATM","ARID1A","SETD2","PTEN","NFE2L2","MGA","U2AF1",
           "CTNNB1","SMAD4","FGFR1","ERBB2","ALK","ROS1","RET"]

# ---------- A. GDC file -> patient mapping ----------
def build_gdc_mapping():
    omics = {"expression":"expression","maf":"maf","copynumber_gene":"copynumber_gene",
             "copynumber_allele":"copynumber_allele","mirna":"mirna","isoform":"isoform"}
    all_fn = {}
    for otype, sub in omics.items():
        d = os.path.join(ALCHEMIST_DIR, "data", sub)
        if os.path.isdir(d):
            for f in os.listdir(d):
                all_fn[f] = otype
    filenames = list(all_fn)
    print(f"[A] {len(filenames)} assay files; querying GDC ...")
    mapping = {}
    B = 150
    for i in range(0, len(filenames), B):
        batch = filenames[i:i+B]
        payload = json.dumps({
            "filters":{"op":"in","content":{"field":"file_name","value":batch}},
            "fields":"file_name,cases.submitter_id,associated_entities.entity_submitter_id,"
                     "associated_entities.case_id,data_type,experimental_strategy",
            "format":"json","size":str(len(batch)+10)})
        r = subprocess.run(["curl","-s","--max-time","90","-H","Content-Type: application/json",
                            "-d",payload,GDC_FILES_API], capture_output=True, text=True, timeout=120)
        for h in json.loads(r.stdout).get("data",{}).get("hits",[]):
            fn = h.get("file_name"); cases = h.get("cases",[{}]); ae = h.get("associated_entities",[{}])
            mapping[fn] = {"patient": cases[0].get("submitter_id") if cases else None,
                           "entity": ae[0].get("entity_submitter_id") if ae else None,
                           "omics": all_fn[fn]}
        time.sleep(0.3)
    json.dump(mapping, open(os.path.join(WORK,"file_patient_mapping.json"),"w"))
    print(f"[A] mapped {len(mapping)}/{len(filenames)} files")
    return mapping

# ---------- B. ALCHEMIST clinical ----------
def parse_clinical():
    files = glob.glob(os.path.join(ALCHEMIST_DIR,"data","clinical","*.json"))
    rows = {}
    for f in files:
        try:
            d = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        cd = d.get("ClinicalData",{}); sd = cd.get("SubjectData",{})
        rec = {"submitter_id": sd.get("submitter_id","")}
        for ev in sd.get("StudyEventData",[]):
            for form in ev.get("FormData",[]):
                for ig in form.get("ItemGroupData",[]):
                    for it in ig.get("ItemData",[]):
                        if it.get("FieldDisplayName"): rec[it["FieldDisplayName"]] = it.get("@Value","")
        for e in cd.get("GDC Clinical Data Entities",[]):
            t = e.get("type","")
            for k,v in e.items():
                if not isinstance(v,(dict,list)) and k not in ("type","projects"):
                    rec[f"{t}.{k}"] = v
        pid = rec.get("submitter_id") or os.path.basename(f).replace(".json","")
        rows[pid] = rec
    return rows

def pick_file(cands):
    # prefer primary tumor (TTP), then recurrent (TTR), then any (normal-entity MAF/CNA still hold tumor calls)
    def code(info):
        e = (info["entity"] or "").split("-"); return e[2] if len(e)>2 else ""
    for pref in ("TTP","TTR","NB"):
        for fn,c in cands:
            if c.startswith(pref): return fn
    return cands[0][0] if cands else ""

def build_alchemist_master(mapping, clinical):
    pat_cands = defaultdict(lambda: defaultdict(list))
    for fn,info in mapping.items():
        e = (info["entity"] or "").split("-"); c = e[2] if len(e)>2 else ""
        pat_cands[info["patient"]][info["omics"]].append((fn,c))
    pat_files = {p:{o:pick_file(cs) for o,cs in od.items()} for p,od in pat_cands.items()}
    rows=[]
    for pid,clin in clinical.items():
        of = pat_files.get(pid,{})
        db = clin.get("days_to_birth","")
        age = round(int(db)/-365.25,1) if str(db).lstrip("-").isdigit() else None
        h = clin.get("histologic_subtype","")
        hg = ("LUAD" if ("Adenocarcinoma" in h and "situ" not in h) else
              "LUSC" if "Squamous" in h else
              "Adenosquamous" if "Adenosquamous" in h else "Other")
        rows.append({"patient_id":pid,"treatment_arm":clin.get("ALCH_treatment_arm",""),
                     "histology":h,"histology_group":hg,"sex":clin.get("demographic.sex_at_birth",""),
                     "age_years":age,"primary_diagnosis":clin.get("diagnosis.primary_diagnosis",""),
                     "has_expr":"expression" in of,"has_maf":"maf" in of,
                     "has_cna_gene":"copynumber_gene" in of,"has_mirna":"mirna" in of,
                     "has_isoform":"isoform" in of,
                     "file_expr":of.get("expression",""),"file_maf":of.get("maf",""),
                     "file_cna_gene":of.get("copynumber_gene",""),"file_mirna":of.get("mirna","")})
    df = pd.DataFrame(rows)
    df["core3"] = df.has_expr & df.has_maf & df.has_cna_gene
    df.to_csv(os.path.join(DATA,"alchemist_master_clinical.csv"), index=False, encoding="utf-8-sig")
    print(f"[B] master clinical: {len(df)} patients; core3={int(df.core3.sum())}")
    return df

# ---------- C. matrices ----------
def build_expression(df):
    EXPR = os.path.join(ALCHEMIST_DIR,"data","expression")
    pts = df[df.has_expr][["patient_id","file_expr"]].values
    d0 = pd.read_csv(os.path.join(EXPR,pts[0][1]), sep="\t", skiprows=1)
    d0 = d0[d0.gene_type=="protein_coding"]
    gid = d0.gene_id.values; gmap = dict(zip(d0.gene_id,d0.gene_name))
    tpm, cnt = {}, {}
    for pid,fn in pts:
        d = pd.read_csv(os.path.join(EXPR,fn), sep="\t", skiprows=1,
                        usecols=["gene_id","tpm_unstranded","unstranded"]).set_index("gene_id")
        tpm[pid]=d.tpm_unstranded.reindex(gid).values; cnt[pid]=d.unstranded.reindex(gid).values
    tpm=pd.DataFrame(tpm,index=gid); cnt=pd.DataFrame(cnt,index=gid)
    tpm.index.name=cnt.index.name="gene_id"
    tpm.insert(0,"gene_name",[gmap[g] for g in tpm.index]); cnt.insert(0,"gene_name",[gmap[g] for g in cnt.index])
    try:
        tpm.to_parquet(os.path.join(MATDIR,"expr_tpm.parquet"))
        cnt.to_parquet(os.path.join(MATDIR,"expr_counts.parquet"))
    except Exception as e:
        print("  parquet unavailable, writing tsv.gz only:", e)
    tpm.to_csv(os.path.join(MATDIR,"expr_tpm.tsv.gz"), sep="\t", compression="gzip")
    cnt.to_csv(os.path.join(MATDIR,"expr_counts.tsv.gz"), sep="\t", compression="gzip")
    print(f"[C] expression: {tpm.shape[0]} genes x {tpm.shape[1]-1} samples")

def build_mutation(df):
    MAF = os.path.join(ALCHEMIST_DIR,"data","maf")
    pts = df[df.has_maf][["patient_id","file_maf"]].values
    key = ["Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification",
           "Variant_Type","Reference_Allele","Tumor_Seq_Allele2","HGVSp_Short"]
    rows=[]; gp=defaultdict(set); tmb=Counter()
    for pid,fn in pts:
        try:
            with gzip.open(os.path.join(MAF,fn),"rt") as f:
                hdr=None; idx={}
                for line in f:
                    if line.startswith("#"): continue
                    if hdr is None:
                        hdr=line.rstrip("\n").split("\t"); idx={c:hdr.index(c) for c in key if c in hdr}
                        vci=hdr.index("Variant_Classification"); gi=hdr.index("Hugo_Symbol"); continue
                    c=line.rstrip("\n").split("\t")
                    if len(c)<len(hdr): continue
                    row={k:c[idx[k]] for k in idx}; row["patient_id"]=pid; row["Tumor_Sample_Barcode"]=pid
                    rows.append(row)
                    if c[vci] in NONSYN: gp[c[gi]].add(pid); tmb[pid]+=1
        except Exception as e:
            print("  MAF error", pid, e)
    pd.DataFrame(rows).to_csv(os.path.join(MATDIR,"alchemist_merged.maf.gz"),
                              sep="\t", index=False, compression="gzip")
    n=len(pts); top=[g for g in gp if len(gp[g])>=n*0.02]
    mat=pd.DataFrame(0,index=top,columns=[p for p,_ in pts])
    for g in top:
        for p in gp[g]:
            if p in mat.columns: mat.loc[g,p]=1
    mat.index.name="gene"; mat.to_csv(os.path.join(MATDIR,"mutation_binary.tsv.gz"),sep="\t",compression="gzip")
    tv=np.array([tmb.get(p,0) for p,_ in pts])
    pd.DataFrame({"patient_id":[p for p,_ in pts],"nonsyn_count":tv,
                  "tmb_per_mb":tv/38.0}).to_csv(os.path.join(MATDIR,"tmb.csv"),index=False)
    print(f"[C] mutation: {mat.shape[0]} genes x {mat.shape[1]} samples; TMB written")

def build_cna(df):
    CNA = os.path.join(ALCHEMIST_DIR,"data","copynumber_gene")
    pts = df[df.has_cna_gene][["patient_id","file_cna_gene"]].values
    d0 = pd.read_csv(os.path.join(CNA,pts[0][1]), sep="\t")
    gid=d0.gene_id.values; gmap=dict(zip(d0.gene_id,d0.gene_name)); data={}
    for pid,fn in pts:
        d=pd.read_csv(os.path.join(CNA,fn),sep="\t",usecols=["gene_id","copy_number"]).set_index("gene_id")
        data[pid]=d.copy_number.reindex(gid).values
    m=pd.DataFrame(data,index=gid); m.index.name="gene_id"
    m.insert(0,"gene_name",[gmap[g] for g in m.index])
    m.to_csv(os.path.join(MATDIR,"cna_gene.tsv.gz"),sep="\t",compression="gzip")
    print(f"[C] CNA: {m.shape[0]} genes x {m.shape[1]-1} samples")

def build_mirna(df):
    MIR = os.path.join(ALCHEMIST_DIR,"data","mirna")
    pts = df[df.has_mirna][["patient_id","file_mirna"]].values
    d0 = pd.read_csv(os.path.join(MIR,pts[0][1]), sep="\t")
    mid=d0.miRNA_ID.values; data={}
    for pid,fn in pts:
        d=pd.read_csv(os.path.join(MIR,fn),sep="\t",
                      usecols=["miRNA_ID","reads_per_million_miRNA_mapped"]).set_index("miRNA_ID")
        data[pid]=d.reads_per_million_miRNA_mapped.reindex(mid).values
    m=pd.DataFrame(data,index=mid); m.index.name="miRNA_ID"
    m.to_csv(os.path.join(MATDIR,"mirna_rpm.tsv.gz"),sep="\t",compression="gzip")
    print(f"[C] miRNA: {m.shape[0]} miRNAs x {m.shape[1]} samples")

# ---------- D. MSK-IMPACT ----------
def build_msk():
    samp = pd.read_csv(os.path.join(MSK_DIR,"data_clinical_sample.txt"),sep="\t",comment="#",low_memory=False)
    pat  = pd.read_csv(os.path.join(MSK_DIR,"data_clinical_patient.txt"),sep="\t",comment="#",low_memory=False)
    nsclc = samp[samp.CANCER_TYPE.str.contains("Non-Small Cell",case=False,na=False)].copy()
    keep = nsclc.CANCER_TYPE_DETAILED.str.contains("Adenocarcinoma|Squamous|Non-Small Cell|Adenosquamous",case=False,na=False)
    nsclc = nsclc[keep].copy()
    msk = nsclc.merge(pat[["PATIENT_ID","OS_STATUS","OS_MONTHS","SEX","AGE_AT_DX"]],on="PATIENT_ID",how="left")
    hg=lambda x:("LUAD" if "Adenocarcinoma" in str(x) else "LUSC" if "Squamous" in str(x)
                 else "Adenosquamous" if "Adenosquamous" in str(x) else "NSCLC_NOS")
    msk["histology_group"]=msk.CANCER_TYPE_DETAILED.apply(hg)
    msk["OS_event"]=(msk.OS_STATUS=="DECEASED").astype(int)
    msk.loc[msk.OS_STATUS.isna(),"OS_event"]=np.nan
    msk["OS_months"]=pd.to_numeric(msk.OS_MONTHS,errors="coerce")
    msk["is_primary"]=(msk.SAMPLE_TYPE=="Primary").astype(int)
    mu=msk.sort_values(["PATIENT_ID","is_primary"],ascending=[True,False]).drop_duplicates("PATIENT_ID")
    cols=["SAMPLE_ID","PATIENT_ID","SAMPLE_TYPE","histology_group","CANCER_TYPE_DETAILED","GENE_PANEL",
          "TMB_SCORE","MSI_TYPE","FACETS_WGD","FACETS_PURITY","SEX","AGE_AT_DX","OS_months","OS_event","OS_STATUS"]
    mu[cols].to_csv(os.path.join(MSKOUT,"msk_nsclc_clinical.csv"),index=False)
    sids=set(mu.SAMPLE_ID)
    # mutation binary (drivers)
    muts=pd.read_csv(os.path.join(MSK_DIR,"data_mutations.txt"),sep="\t",low_memory=False,
                     usecols=["Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode"])
    muts=muts[muts.Variant_Classification.isin(NONSYN) & muts.Tumor_Sample_Barcode.isin(sids)]
    mb=muts[muts.Hugo_Symbol.isin(DRIVERS)].groupby(["Hugo_Symbol","Tumor_Sample_Barcode"]).size().unstack(fill_value=0)
    mb=(mb>0).astype(int).reindex(columns=list(sids),fill_value=0)
    mb.to_csv(os.path.join(MSKOUT,"msk_mutation_binary.csv"))
    # gene CNA + arm CNA (NSCLC subset)
    with open(os.path.join(MSK_DIR,"data_cna.txt")) as f: hdr=f.readline().strip().split("\t")
    usecols=["Hugo_Symbol"]+[c for c in hdr if c in sids]
    pd.read_csv(os.path.join(MSK_DIR,"data_cna.txt"),sep="\t",usecols=usecols,low_memory=False)\
      .to_csv(os.path.join(MSKOUT,"msk_gene_cna.csv"),index=False)
    arm=pd.read_csv(os.path.join(MSK_DIR,"data_armlevel_cna.txt"),sep="\t",low_memory=False)
    arm[[c for c in arm.columns if c in ("NAME","ENTITY_STABLE_ID") or c in sids]]\
      .to_csv(os.path.join(MSKOUT,"msk_armlevel_cna.csv"),index=False)
    print(f"[D] MSK NSCLC: {len(mu)} patients; OS available={int(mu.OS_months.notna().sum())}")

def main():
    mapping = build_gdc_mapping()
    clinical = parse_clinical()
    df = build_alchemist_master(mapping, clinical)
    build_expression(df); build_mutation(df); build_cna(df); build_mirna(df)
    build_msk()
    print("=== 01_prepare_data.py complete ===")

if __name__ == "__main__":
    main()
