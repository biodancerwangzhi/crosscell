#!/usr/bin/env python3
"""
easySCF Round-trip Test -- Using the Paper's Own Test Data (GSE152048)

Paper: Zhang et al. 2024 (easySCF, Genome Analysis)
Test data: GSE152048 (Zhou et al. 2020, Nature Communications) -- Osteosarcoma scRNA-seq
Paper Section 3.6: "We tested this pathway using data from GSE152048"

Tests:
  A: R Seurat -> easySCFr saveH5 -> easySCFpy loadH5 -> AnnData
  B: Python AnnData -> easySCFpy saveH5 -> easySCFr loadH5 -> Seurat
  C: Full roundtrip R -> H5 -> Py -> H5 -> R

Usage: python3 /benchmark/scripts/test_easyscf_roundtrip.py [--skip-download]
"""
import json, os, subprocess, sys, time

GSE152048_SAMPLES = {
    "GSM4600182_BC5": {"desc": "GSE152048 BC5 primary osteoblastic OS", "prefix": "GSM4600182_BC5"},
    "GSM4600189_BC22": {"desc": "GSE152048 BC22 primary chondroblastic OS", "prefix": "GSM4600189_BC22"},
}
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152048/suppl"
PAPER_DIR = "/tmp/easyscf_roundtrip/gse152048"
TMP = "/tmp/easyscf_roundtrip"

BENCH = [
    {"id":"v5_pbmc3k_raw","input":"/benchmark/data/generated/seurat_v5_pbmc3k_raw.rds","type":"rds","desc":"V5 2.7k"},
    {"id":"v4_pbmc3k_proc","input":"/benchmark/data/generated/seurat_v4_pbmc3k_processed.rds","type":"rds","desc":"V4 proc 2.7k"},
    {"id":"v5_ifnb_raw","input":"/benchmark/data/generated/seurat_v5_ifnb_raw.rds","type":"rds","desc":"V5 14k"},
    {"id":"v4_cbmc_raw","input":"/benchmark/data/generated/seurat_v4_cbmc_raw.rds","type":"rds","desc":"V4 CITE 8.6k"},
    {"id":"v5_bmcite_raw","input":"/benchmark/data/generated/seurat_v5_bmcite_raw.rds","type":"rds","desc":"V5 30.7k"},
    {"id":"cxg_pbmc15k","input":"/benchmark/data/generated/cellxgene_pbmc_15k.h5ad","type":"h5ad","desc":"CxG 14.8k"},
    {"id":"scanpy_pbmc3k","input":"/benchmark/data/generated/scanpy_pbmc3k.h5ad","type":"h5ad","desc":"Scanpy 2.7k"},
    {"id":"scvelo_panc","input":"/benchmark/data/generated/scvelo_pancreas.h5ad","type":"h5ad","desc":"scVelo 3.7k"},
]

def run(cmd, timeout=600):
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode == 0, r.stdout, r.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return False, "", "TIMEOUT", timeout


def dl_sample(sid, info):
    """Download one GSE152048 sample from GEO, process to Seurat RDS."""
    sdir = os.path.join(PAPER_DIR, sid)
    rds = os.path.join(PAPER_DIR, f"{sid}.rds")
    if os.path.exists(rds):
        print(f"    RDS exists: {rds}"); return rds
    os.makedirs(sdir, exist_ok=True)
    pfx = info["prefix"]
    fnames = [f"{pfx}_barcodes.tsv.gz", f"{pfx}_features.tsv.gz", f"{pfx}_matrix.mtx.gz"]
    print(f"    Downloading {sid} from GEO...")
    for fn in fnames:
        lp = os.path.join(sdir, fn)
        if os.path.exists(lp) and os.path.getsize(lp) > 0:
            print(f"      Exists: {fn}"); continue
        url = f"{GEO_FTP}/{fn}"
        print(f"      {fn}...")
        ok, _, e, t = run(["wget", "-q", "-O", lp, url], 300)
        if not ok: ok, _, e, t = run(["curl", "-sL", "-o", lp, url], 300)
        if not ok: print(f"      FAIL: {e[:200]}"); return None
        print(f"      OK ({os.path.getsize(lp)/1048576:.1f}MB {t:.0f}s)")
    # Symlinks for Read10X
    ld = os.path.join(sdir, "filtered"); os.makedirs(ld, exist_ok=True)
    for s, d in [(f"{pfx}_barcodes.tsv.gz","barcodes.tsv.gz"),
                 (f"{pfx}_features.tsv.gz","features.tsv.gz"),
                 (f"{pfx}_matrix.mtx.gz","matrix.mtx.gz")]:
        dp = os.path.join(ld, d)
        if os.path.exists(dp): os.remove(dp)
        os.symlink(os.path.join(sdir, s), dp)
    # Seurat processing
    print(f"    Seurat processing...")
    rs = f'''library(Seurat)
d<-Read10X("{ld}"); o<-CreateSeuratObject(counts=d,project="{sid}",min.cells=3,min.features=200)
o[["pmt"]]<-PercentageFeatureSet(o,pattern="^MT-")
o<-subset(o,subset=nFeature_RNA>300&pmt<15)
o<-NormalizeData(o); o<-FindVariableFeatures(o,nfeatures=2000)
o<-ScaleData(o); o<-RunPCA(o,npcs=30)
o<-FindNeighbors(o,dims=1:20); o<-FindClusters(o,resolution=0.5)
o<-RunUMAP(o,dims=1:20)
cat("Seurat:",ncol(o),"x",nrow(o),"\\n"); saveRDS(o,"{rds}")'''
    ok, out, err, t = run(["Rscript", "-e", rs], 600)
    if not ok: print(f"    FAIL: {(out+err)[:400]}"); return None
    print(f"    Done ({t:.0f}s) {out.strip()}")
    return rds

def mk_h5ad(sid, info):
    """Create H5AD from 10x data for Py->R test."""
    hp = os.path.join(PAPER_DIR, f"{sid}.h5ad")
    if os.path.exists(hp): print(f"    H5AD exists"); return hp
    ld = os.path.join(PAPER_DIR, sid, "filtered")
    if not os.path.exists(ld): return None
    print(f"    Creating H5AD...")
    ps = f'''import scanpy as sc
a=sc.read_10x_mtx("{ld}",var_names="gene_symbols",cache=False)
sc.pp.filter_cells(a,min_genes=300); sc.pp.filter_genes(a,min_cells=3)
a.var["mt"]=a.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(a,qc_vars=["mt"],inplace=True)
a=a[a.obs.pct_counts_mt<15,:].copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a,n_top_genes=2000)
sc.pp.pca(a,n_comps=30); sc.pp.neighbors(a,n_pcs=20)
sc.tl.umap(a); sc.tl.leiden(a,resolution=0.5)
print(f"AnnData: {{a.n_obs}}x{{a.n_vars}}"); a.write("{hp}")'''
    ok, out, err, t = run(["python3", "-c", ps], 600)
    if ok: print(f"    Done ({t:.0f}s) {out.strip()}"); return hp
    print(f"    FAIL: {(out+err)[:300]}"); return None

def prep_paper(skip):
    """Download and prepare GSE152048 paper data."""
    cases = []
    if skip:
        print("\n  Skip download (--skip-download)")
        for sid, info in GSE152048_SAMPLES.items():
            r = os.path.join(PAPER_DIR, f"{sid}.rds")
            if os.path.exists(r):
                cases.append({"id":f"paper_{sid}","input":r,"type":"rds",
                              "desc":f"[PAPER] {info['desc']}","paper":True})
            h = os.path.join(PAPER_DIR, f"{sid}.h5ad")
            if os.path.exists(h):
                cases.append({"id":f"paper_{sid}_h5ad","input":h,"type":"h5ad",
                              "desc":f"[PAPER] {info['desc']} H5AD","paper":True})
        return cases
    print("\n" + "="*70)
    print("  Downloading GSE152048 (Zhou et al. 2020, Nat Commun)")
    print("="*70)
    os.makedirs(PAPER_DIR, exist_ok=True)
    for sid, info in GSE152048_SAMPLES.items():
        print(f"\n  [{sid}] {info['desc']}")
        r = dl_sample(sid, info)
        if r: cases.append({"id":f"paper_{sid}","input":r,"type":"rds",
                            "desc":f"[PAPER] {info['desc']}","paper":True})
        h = mk_h5ad(sid, info)
        if h: cases.append({"id":f"paper_{sid}_h5ad","input":h,"type":"h5ad",
                            "desc":f"[PAPER] {info['desc']} H5AD","paper":True})
    return cases


def test_a(case):
    """Test A: R -> H5 -> Python"""
    h5 = os.path.join(TMP, f"{case['id']}_a.h5")
    vf = os.path.join(TMP, f"{case['id']}_va.json")
    tag = "[PAPER]" if case.get("paper") else "[BENCH]"
    print(f"\n  [A] R->H5->Py: {case['id']} {tag}")
    print(f"    {case['desc']}")
    # R saveH5
    rs = f'''library(easySCFr);library(Seurat)
o<-readRDS("{case['input']}");o<-UpdateSeuratObject(o)
cat("Seurat:",ncol(o),"x",nrow(o),"\\n")
saveH5(o,"{h5}");cat("saved\\n")'''
    ok,out,err,t = run(["Rscript","-e",rs])
    print(f"    R saveH5: {'OK' if ok else 'FAIL'} ({t:.1f}s)")
    if not ok:
        print(f"    {(out+err)[:300]}")
        return {"test":"A","id":case["id"],"status":"failed","step":"r_save",
                "error":(out+err)[:800],"paper":case.get("paper",False)}
    # Py loadH5
    ps = f'''import json;from easySCFpy import loadH5
try:
    a=loadH5("{h5}")
    d={{"n_obs":a.n_obs,"n_vars":a.n_vars,"obsm":list(a.obsm.keys()) if a.obsm else [],"ok":True}}
except Exception as e: d={{"ok":False,"err":str(e)}}
with open("{vf}","w") as f: json.dump(d,f)
print(json.dumps(d))'''
    ok2,out2,err2,t2 = run(["python3","-c",ps])
    print(f"    Py loadH5: {'OK' if ok2 else 'FAIL'} ({t2:.1f}s)")
    if not ok2:
        print(f"    {(out2+err2)[:300]}")
        return {"test":"A","id":case["id"],"status":"failed","step":"py_load",
                "error":(out2+err2)[:800],"paper":case.get("paper",False)}
    try:
        with open(vf) as f: d=json.load(f)
        if not d.get("ok"):
            return {"test":"A","id":case["id"],"status":"failed","step":"py_verify",
                    "error":d.get("err",""),"paper":case.get("paper",False)}
        print(f"    AnnData: {d['n_obs']}x{d['n_vars']}")
        return {"test":"A","id":case["id"],"status":"success","info":d,
                "t_r":round(t,2),"t_py":round(t2,2),"paper":case.get("paper",False)}
    except Exception as e:
        return {"test":"A","id":case["id"],"status":"failed","step":"verify",
                "error":str(e),"paper":case.get("paper",False)}

def test_b(case):
    """Test B: Python -> H5 -> R"""
    h5 = os.path.join(TMP, f"{case['id']}_b.h5")
    vf = os.path.join(TMP, f"{case['id']}_vb.json")
    tag = "[PAPER]" if case.get("paper") else "[BENCH]"
    print(f"\n  [B] Py->H5->R: {case['id']} {tag}")
    print(f"    {case['desc']}")
    # Py saveH5
    ps = f'''import scanpy as sc;from easySCFpy import saveH5
a=sc.read_h5ad("{case['input']}")
print(f"AnnData: {{a.n_obs}}x{{a.n_vars}}")
saveH5(a,"{h5}");print("saved")'''
    ok,out,err,t = run(["python3","-c",ps])
    print(f"    Py saveH5: {'OK' if ok else 'FAIL'} ({t:.1f}s)")
    if not ok:
        print(f"    {(out+err)[:300]}")
        return {"test":"B","id":case["id"],"status":"failed","step":"py_save",
                "error":(out+err)[:800],"paper":case.get("paper",False)}
    # R loadH5
    rs = f'''library(easySCFr);library(Seurat);library(jsonlite)
tryCatch({{
  o<-loadH5("{h5}")
  d<-list(n_cells=ncol(o),n_genes=nrow(o),status="ok")
  writeLines(toJSON(d,auto_unbox=TRUE),"{vf}")
  cat("Seurat:",ncol(o),"x",nrow(o),"\\n")
}},error=function(e){{
  d<-list(status="fail",error=conditionMessage(e))
  writeLines(toJSON(d,auto_unbox=TRUE),"{vf}")
  cat("Err:",conditionMessage(e),"\\n")
}})'''
    ok2,out2,err2,t2 = run(["Rscript","-e",rs])
    print(f"    R loadH5: {'OK' if ok2 else 'FAIL'} ({t2:.1f}s)")
    try:
        with open(vf) as f: d=json.load(f)
        if d.get("status")!="ok":
            em = d.get("error","unknown")
            print(f"    FAIL: {em[:200]}")
            return {"test":"B","id":case["id"],"status":"failed","step":"r_load",
                    "error":em,"paper":case.get("paper",False)}
        print(f"    Seurat: {d['n_cells']}x{d['n_genes']}")
        return {"test":"B","id":case["id"],"status":"success","info":d,
                "t_py":round(t,2),"t_r":round(t2,2),"paper":case.get("paper",False)}
    except Exception as e:
        if not ok2:
            return {"test":"B","id":case["id"],"status":"failed","step":"r_load",
                    "error":(out2+err2)[:800],"paper":case.get("paper",False)}
        return {"test":"B","id":case["id"],"status":"failed","step":"verify",
                "error":str(e),"paper":case.get("paper",False)}


def test_c(case):
    """Test C: R -> H5 -> Py -> H5 -> R (full roundtrip)"""
    h1 = os.path.join(TMP, f"{case['id']}_c1.h5")
    h2 = os.path.join(TMP, f"{case['id']}_c2.h5")
    vf = os.path.join(TMP, f"{case['id']}_vc.json")
    tag = "[PAPER]" if case.get("paper") else "[BENCH]"
    print(f"\n  [C] R->H5->Py->H5->R: {case['id']} {tag}")
    print(f"    {case['desc']}")
    # Step1: R saveH5
    r1 = f'''library(easySCFr);library(Seurat)
o<-readRDS("{case['input']}");o<-UpdateSeuratObject(o)
cat("Orig:",ncol(o),"x",nrow(o),"\\n")
saveH5(o,"{h1}");cat("s1 done\\n")'''
    ok1,o1,e1,t1 = run(["Rscript","-e",r1])
    print(f"    S1 R saveH5:   {'OK' if ok1 else 'FAIL'} ({t1:.1f}s)")
    if not ok1:
        return {"test":"C","id":case["id"],"status":"failed","step":1,
                "error":(o1+e1)[:500],"paper":case.get("paper",False)}
    # Step2: Py load+save
    p2 = f'''from easySCFpy import loadH5,saveH5
a=loadH5("{h1}");print(f"Py: {{a.n_obs}}x{{a.n_vars}}")
saveH5(a,"{h2}");print("s2 done")'''
    ok2,o2,e2,t2 = run(["python3","-c",p2])
    print(f"    S2 Py load+save: {'OK' if ok2 else 'FAIL'} ({t2:.1f}s)")
    if not ok2:
        return {"test":"C","id":case["id"],"status":"failed","step":2,
                "error":(o2+e2)[:500],"paper":case.get("paper",False)}
    # Step3: R load+compare
    r3 = f'''library(easySCFr);library(Seurat);library(jsonlite)
tryCatch({{
  orig<-readRDS("{case['input']}");orig<-UpdateSeuratObject(orig)
  rt<-loadH5("{h2}")
  oc<-ncol(orig);og<-nrow(orig);rc<-ncol(rt);rg<-nrow(rt)
  n<-min(100,og,rg);m<-min(100,oc,rc)
  ov<-as.vector(as.matrix(GetAssayData(orig,layer="counts")[1:n,1:m]))
  rv<-as.vector(as.matrix(GetAssayData(rt,layer="counts")[1:n,1:m]))
  cr<-if(length(ov)==length(rv)) cor(ov,rv) else NA
  d<-list(status="ok",oc=oc,og=og,rc=rc,rg=rg,cm=(oc==rc),gm=(og==rg),corr=cr)
  writeLines(toJSON(d,auto_unbox=TRUE),"{vf}")
  cat("Orig:",oc,"x",og,"RT:",rc,"x",rg,"Corr:",cr,"\\n")
}},error=function(e){{
  d<-list(status="fail",error=conditionMessage(e))
  writeLines(toJSON(d,auto_unbox=TRUE),"{vf}")
  cat("Err:",conditionMessage(e),"\\n")
}})'''
    ok3,o3,e3,t3 = run(["Rscript","-e",r3])
    print(f"    S3 R load+cmp:   {'OK' if ok3 else 'FAIL'} ({t3:.1f}s)")
    try:
        with open(vf) as f: d=json.load(f)
        if d.get("status")!="ok":
            print(f"    FAIL: {d.get('error','')[:200]}")
            return {"test":"C","id":case["id"],"status":"failed","step":3,
                    "error":d.get("error",""),"paper":case.get("paper",False)}
        st = "success" if d.get("cm") and d.get("gm") else "mismatch"
        print(f"    {'PASS' if st=='success' else 'MISMATCH'} "
              f"Orig:{d['oc']}x{d['og']} RT:{d['rc']}x{d['rg']} Corr:{d.get('corr')}")
        return {"test":"C","id":case["id"],"status":st,"info":d,
                "times":{"s1":round(t1,2),"s2":round(t2,2),"s3":round(t3,2)},
                "paper":case.get("paper",False)}
    except Exception as e:
        if not ok3:
            return {"test":"C","id":case["id"],"status":"failed","step":3,
                    "error":(o3+e3)[:500],"paper":case.get("paper",False)}
        return {"test":"C","id":case["id"],"status":"failed","step":"verify",
                "error":str(e),"paper":case.get("paper",False)}


def main():
    skip = "--skip-download" in sys.argv
    print("="*70)
    print("  easySCF Round-trip Test")
    print("  Paper data: GSE152048 (Zhou et al. 2020)")
    print("="*70)
    os.makedirs(TMP, exist_ok=True)

    paper = prep_paper(skip)
    print(f"\n  Paper data: {len(paper)} cases")

    all_c = paper + BENCH
    rds_c = [c for c in all_c if c["type"]=="rds"]
    h5ad_c = [c for c in all_c if c["type"]=="h5ad"]
    res = []

    print("\n"+"="*70)
    print("  Test A: R -> H5 -> Python")
    print("="*70)
    for c in rds_c: res.append(test_a(c))

    print("\n"+"="*70)
    print("  Test B: Python -> H5 -> R (expected FAIL)")
    print("="*70)
    for c in h5ad_c: res.append(test_b(c))

    print("\n"+"="*70)
    print("  Test C: Full roundtrip R -> H5 -> Py -> H5 -> R")
    print("="*70)
    for c in rds_c: res.append(test_c(c))

    # Summary
    print("\n"+"="*70)
    print("  SUMMARY")
    print("="*70)
    pr = [r for r in res if r.get("paper")]
    br = [r for r in res if not r.get("paper")]
    for lbl, sub in [("Paper (GSE152048)", pr), ("Benchmark", br)]:
        if not sub: continue
        print(f"\n  --- {lbl} ---")
        for r in sub:
            s = r["status"]
            sym = {"success":"PASS","failed":"FAIL","mismatch":"MISMATCH"}.get(s,"?")
            step = f" (step {r.get('step','?')})" if s=="failed" else ""
            print(f"  {sym} [{r['test']}] {r['id']}: {s}{step}")
            if s=="failed": print(f"       {r.get('error','')[:200]}")

    for lbl, sub in [("Paper",pr),("Bench",br),("Total",res)]:
        if not sub: continue
        ok=sum(1 for r in sub if r["status"]=="success")
        fl=sum(1 for r in sub if r["status"]=="failed")
        mm=sum(1 for r in sub if r["status"]=="mismatch")
        print(f"\n  {lbl}: {len(sub)} tests, {ok} pass, {fl} fail, {mm} mismatch")

    for tt in "ABC":
        sub=[r for r in res if r["test"]==tt]
        if sub:
            ok=sum(1 for r in sub if r["status"]=="success")
            print(f"  Test {tt}: {ok}/{len(sub)} pass")

    out = "/benchmark/results/easyscf_roundtrip.json"
    with open(out,"w") as f: json.dump(res,f,indent=2,ensure_ascii=False)
    print(f"\n  Saved: {out}")

    # Key findings
    ta=[r for r in res if r["test"]=="A"]
    tb=[r for r in res if r["test"]=="B"]
    tc=[r for r in res if r["test"]=="C"]
    print("\n"+"="*70)
    print("  KEY FINDINGS")
    print("="*70)
    print(f"  R->H5->Py (A):        {sum(1 for r in ta if r['status']=='success')}/{len(ta)}")
    print(f"  Py->H5->R (B):        {sum(1 for r in tb if r['status']=='success')}/{len(tb)}")
    print(f"  R->H5->Py->H5->R (C): {sum(1 for r in tc if r['status']=='success')}/{len(tc)}")
    if tb and sum(1 for r in tb if r["status"]=="success")==0:
        print("\n  CONCLUSION: easySCF Py->R is BROKEN")
        print("  easySCFpy uses anndata write_elem for var encoding,")
        print("  but easySCFr expects var/rawvar/_index as simple HDF5 dataset.")
    if tc and sum(1 for r in tc if r["status"]=="success")==0:
        print("\n  CONCLUSION: Full roundtrip FAILS at Py->H5->R leg")
        print("  easySCF is one-way only: R -> H5 -> Python")

if __name__ == "__main__":
    main()
