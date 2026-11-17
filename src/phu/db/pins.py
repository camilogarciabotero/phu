"""The phu-db bundle this phu release installs.

phu-db-0.4.0, Zenodo record 24296611 (https://doi.org/10.5281/zenodo.24296611),
built by phu-db-prep 0.4.0. The vendored phu/data/phu-db.json is the manifest
inside that tarball; its digest is pinned here too.
"""

from phu.db.bundle import BundlePins

PRODUCTION_PINS = BundlePins(
    url="https://zenodo.org/records/24296611/files/phu-db.tar.gz",
    sha256="93bd0b9d2f5efb85064b36a4a8d879f5b8689ec92e677d391be1c14fdc25e6e5",
    size=2885279568,
    bundle_version="phu-db-0.4.0",
    manifest_sha256="a1f3bfd79f26f87cbecde33010171cc2a0f016c438020722860f45479f4827fe",
    expected_members=frozenset(
        {
            "phu-db/phu-db.json",
            "phu-db/pfam/profiles.hmm",
            "phu-db/kofam/profiles.hmm",
            "phu-db/kofam/ko_list.tsv",
            "phu-db/dbcan/profiles.hmm",
            "phu-db/dbcan/family_substrate.tsv",
            "phu-db/dbcan/pul_rules.tsv",
            "phu-db/avg/v_scores.tsv",
            "phu-db/avg/positive_amg.tsv",
            "phu-db/avg/positive_apg.tsv",
            "phu-db/avg/positive_areg.tsv",
            "phu-db/avg/filters_amg.tsv",
            "phu-db/avg/filters_apg.tsv",
            "phu-db/avg/filters_areg.tsv",
            "phu-db/phrog/profiles.hmm",
            "phu-db/phrog/msa_a3m.tar.gz",
            "phu-db/phrog/annotations.tsv",
            "phu-db/ecod/profiles.hmm",
            "phu-db/ecod/msa_a3m.tar.gz",
            "phu-db/ecod/sequences.fa",
            "phu-db/ecod/families.tsv",
            "phu-db/ecod/domains.tsv",
        }
    ),
)
