"""Compare a phu avger run with CheckAMG annotate's shipped v1.2 outputs.

Proteins are matched by (contig, gene ordinal): both tools call genes with
pyrodigal-gv, and CheckAMG names proteins <contig>_<n> while phu writes the
contig and gene_number columns. The match rate is reported so a numbering
drift is visible.
"""
import csv, re, sys
from collections import Counter

phu_dir, ca_dir = sys.argv[1], sys.argv[2]

def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

# CheckAMG predictions: key by (contig, ordinal).
ca_pred = {}
for r in load(f"{ca_dir}/avg_predictions.tsv"):
    m = re.match(r"^(.*)_(\d+)$", r["Protein"])
    ca_pred[(r["Contig"], int(m.group(2)))] = r
ca_all = {}
for r in load(f"{ca_dir}/detailed/all_proteins.tsv"):
    m = re.match(r"^(.*)_(\d+)$", r["Protein"])
    ca_all[(r["Contig"], int(m.group(2)))] = r

phu_rows = load(f"{phu_dir}/gene-annotations.tsv")
phu = {(r["contig"], int(r["gene_number"])): r for r in phu_rows}

shared = set(phu) & set(ca_all)
print(f"proteins: phu={len(phu)} checkamg={len(ca_all)} matched={len(shared)}")

CLASSMAP = {"putative_amg": "metabolic", "putative_apg": "physiological", "putative_areg": "regulatory"}
phu_calls = {k: CLASSMAP[r["avg_class"]] for k, r in phu.items() if r["classification_status"] == "classified"}
ca_calls = {k: r["Protein Classification"].strip().lower() for k, r in ca_all.items() if r["Protein Classification"].strip().lower() in ("metabolic", "physiological", "regulatory")}
print(f"calls: phu={len(phu_calls)} {dict(Counter(phu_calls.values()))}")
print(f"calls: checkamg={len(ca_calls)} {dict(Counter(ca_calls.values()))}")

both = {k for k in phu_calls if k in ca_calls}
agree = {k for k in both if phu_calls[k] == ca_calls[k]}
print(f"called by both={len(both)} same class={len(agree)} "
      f"class differs={sorted((phu_calls[k], ca_calls[k]) for k in both - agree)}")

def phu_status(k):
    r = phu.get(k)
    if r is None: return "not matched"
    s = r["classification_status"]
    if s == "classified": return f"phu calls {CLASSMAP[r['avg_class']]}"
    if s == "not_candidate":
        has = any(r[c] for c in ("pfam_id", "kofam_id", "phrog_id"))
        return "fails score rule" + ("" if has else ", no strict annotation")
    return s

only_ca = Counter(f"{ca_calls[k]}: {phu_status(k)}" for k in ca_calls if k not in phu_calls)
print("\nCheckAMG-only calls, by phu's view:")
for reason, n in only_ca.most_common(): print(f"  {n:3d}  {reason}")

only_phu = Counter()
for k in phu_calls:
    if k in ca_calls: continue
    r = ca_all.get(k)
    if r is None: only_phu[f"{phu_calls[k]}: not in CheckAMG output"] += 1; continue
    conf = r["Protein Viral Origin Confidence"]
    only_phu[f"{phu_calls[k]}: CheckAMG {r['Protein Classification'] or 'unclassified'} (viral-origin {conf})"] += 1
print("\nphu-only calls, by CheckAMG's view:")
for reason, n in only_phu.most_common(): print(f"  {n:3d}  {reason}")

# Annotation agreement on shared proteins, per database.
ca_ann = {}
for r in load(f"{ca_dir}/detailed/gene_annotations.tsv"):
    m = re.match(r"^(.*)_(\d+)$", r["Protein"])
    ca_ann[(r["Contig"], int(m.group(2)))] = r
for db, phu_col, ca_col, norm in (
    ("pfam", "pfam_id", "Pfam_hmm_id", lambda v: v.split(".")[0].upper()),
    ("kofam", "kofam_id", "KEGG_hmm_id", str.upper),
    ("phrog", "phrog_id", "PHROG_hmm_id", lambda v: v.lower().replace("phrog_", "")),
):
    pairs = [(norm(phu[k][phu_col]) if phu[k][phu_col] else "", norm(ca_ann[k][ca_col]) if ca_ann[k][ca_col] else "") for k in shared if k in ca_ann]
    both_ann = [(a, b) for a, b in pairs if a and b]
    same = sum(1 for a, b in both_ann if a == b)
    print(f"\n{db}: phu annotates {sum(1 for a,_ in pairs if a)}, checkamg {sum(1 for _,b in pairs if b)}, "
          f"both {len(both_ann)}, same model {same}")
