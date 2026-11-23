#!/usr/bin/env bash
# Run real phu commands from an installed wheel, the way users install it:
# `uv tool install`, with a minimal PATH and outside the source tree, so
# neither the checkout nor the dev venv can stand in for the package.
set -euo pipefail

wheel=$(ls "$PWD"/dist/*.whl)
quickstart="$PWD/docs/quickstart"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
export UV_TOOL_DIR="$work/tools" UV_TOOL_BIN_DIR="$work/bin"
uv tool install --python "${PYTHON_VERSION:-3.12}" "$wheel"
tool_python="$UV_TOOL_DIR/phu/bin/python"
phu="$UV_TOOL_BIN_DIR/phu"
cd "$work"

# Synthetic inputs: genomes for cluster; a contig with one planted gene, a
# contig without it, and an HMM built from that gene's protein for screen.
"$tool_python" - <<'PY'
import random

import pyhmmer

rng = random.Random(7)


def dna(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


a, b = dna(3000), dna(3000)
a_copy = "".join(c if rng.random() > 0.005 else rng.choice("ACGT") for c in a)
with open("genomes.fna", "w") as out:
    out.write(f">genome_a\n{a}\n>genome_a_copy\n{a_copy}\n>genome_b\n{b}\n")

bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code = {x + y + z: aas[16 * i + 4 * j + k]
        for i, x in enumerate(bases) for j, y in enumerate(bases)
        for k, z in enumerate(bases)}
sense = [c for c, aa in code.items() if aa != "*"]
gene = "ATG" + "".join(rng.choice(sense) for _ in range(299)) + "TAA"
protein = "".join(code[gene[i:i + 3]] for i in range(0, len(gene) - 3, 3))
with open("contigs.fna", "w") as out:
    out.write(f">hit\n{dna(400)}{gene}{dna(400)}\n>miss\n{dna(1700)}\n")

alphabet = pyhmmer.easel.Alphabet.amino()
seq = pyhmmer.easel.TextSequence(name=b"planted", sequence=protein).digitize(alphabet)
hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(seq, pyhmmer.plan7.Background(alphabet))
with open("planted.hmm", "wb") as out:
    hmm.write(out)
PY

export PATH=/usr/bin:/bin
"$phu" --version

"$phu" screen -i contigs.fna -o screen-out planted.hmm
test "$(cat screen-out/matched-contig-ids.txt)" = "hit"

# The documented quickstart (docs/quickstart.md) gives its documented result.
"$phu" screen -i "$quickstart/contigs.fna" -o qs-out "$quickstart/major_head.hmm"
test "$(cat qs-out/matched-contig-ids.txt)" = "NC_001422.1"

"$phu" cluster -i genomes.fna -o cluster-out --mode dereplication
test "$(grep -c "^genome_" cluster-out/cluster-representative-ids.txt)" -eq 2
echo "installed-wheel smoke runs passed"
