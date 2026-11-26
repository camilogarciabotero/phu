"""Rebuild the quickstart files from the phiX174 GenBank record.

    curl -o phix.gb "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=gb&retmode=text"
    python make_quickstart.py phix.gb

Writes, next to this script:
- contigs.fna: phiX174 (NC_001422.1) and `decoy`, the same bases shuffled
  (seed 1), so one contig should be kept and one should not;
- major_head.hmm: a profile HMM of the phiX174 major head protein (gene F),
  built from that single sequence with PyHMMER.
"""

import random
import re
import sys
from pathlib import Path

import pyhmmer

record = Path(sys.argv[1]).read_text()
here = Path(__file__).parent

origin = record.split("\nORIGIN", 1)[1].split("//", 1)[0]
genome = re.sub(r"[^acgt]", "", origin).upper()
assert len(genome) == 5386, len(genome)

feature = record.split('/product="major head protein"', 1)[1]
protein = re.sub(r"\s", "", feature.split('/translation="', 1)[1].split('"', 1)[0])

bases = list(genome)
random.Random(1).shuffle(bases)


def wrap(seq: str) -> str:
    return "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))


(here / "contigs.fna").write_text(
    f">NC_001422.1 Escherichia phage phiX174\n{wrap(genome)}\n"
    f">decoy phiX174 bases shuffled\n{wrap(''.join(bases))}\n"
)

alphabet = pyhmmer.easel.Alphabet.amino()
seq = pyhmmer.easel.TextSequence(name=b"phiX174_F", sequence=protein).digitize(alphabet)
builder = pyhmmer.plan7.Builder(alphabet, seed=42)
hmm, _, _ = builder.build(seq, pyhmmer.plan7.Background(alphabet))
hmm.name = b"major_head"
hmm.description = b"phiX174 major head protein (gene F), from one sequence"
# pyhmmer records sys.argv and the current time in the model header. Both
# would put the maintainer's local path and the build date into a file that
# ships with the docs, and would change the file on every rebuild.
hmm.command_line = None
hmm.creation_time = None
with open(here / "major_head.hmm", "wb") as out:
    hmm.write(out)
print(f"genome {len(genome)} bp, protein F {len(protein)} aa")
