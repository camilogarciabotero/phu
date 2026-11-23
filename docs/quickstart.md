# Quickstart

A first run in a few minutes, on two small files, with no database download.

## 1. Get the example files

Two files, about 210 KB together:

- [`contigs.fna`](quickstart/contigs.fna): the genome of phage phiX174
  (`NC_001422.1`, 5,386 bp) and `decoy`, the same bases shuffled.
- [`major_head.hmm`](quickstart/major_head.hmm): a profile HMM of the phiX174
  major head protein (gene F).

```bash
mkdir phu-quickstart && cd phu-quickstart
base=https://camilogarciabotero.github.io/phu/quickstart
curl -O $base/contigs.fna -O $base/major_head.hmm
```

They were made by [`make_quickstart.py`](quickstart/make_quickstart.py) from
the GenBank record.

## 2. Screen the contigs with your own HMM

```bash
phu screen -i contigs.fna -o qs-out major_head.hmm
```

This takes a few seconds. phu predicts the proteins of both contigs, searches
them with the HMM and keeps the contigs with a hit.

## 3. Check the result

```bash
cat qs-out/matched-contig-ids.txt
```

```text
NC_001422.1
```

phiX174 is kept and the decoy is not: its shuffled bases encode no head
protein. The output folder holds:

| File | Content |
| --- | --- |
| `matched-contig-ids.txt` | IDs of the kept contigs |
| `matched-contigs.fna` | Their sequences |
| `major_head-hmmsearch-domains.tbl` | The hmmsearch domain table: the gene F hit scores about 1,000 bits |
| `.phu/query.json` | What was searched and with which settings |

## 4. Search a named family instead

Instead of your own HMM, name a family from the phu databases. Pfam
`PF02305` is the phiX174 capsid F protein family:

```bash
phu screen -i contigs.fna -o pf-out PF02305
cat pf-out/matched-contig-ids.txt   # NC_001422.1
```

The first run that names a database family installs the phu databases:
a **2.9 GB download** that needs about **20 GB free** while it installs and
keeps about **21 GB**. On a cluster, install them once on shared storage
first; see [Databases](databases.md).

## Next

- Your own contigs: [`screen`](commands/screen.md) for named families,
  [`jack`](commands/jack.md) to find relatives of proteins you have.
- Families by name: [`phrogs`](commands/phrogs.md), e.g.
  `phu screen -i contigs.fa --combine-mode all alias:mcp alias:terl`.
- Candidate auxiliary genes: [`avger`](commands/avger.md).
