"""The phu-db bundle this phu release installs.

phu-db-0.3.0, Zenodo record 24290671 (https://doi.org/10.5281/zenodo.24290671),
built by phu-db-prep 0.3.0. The vendored phu/data/phu-db.json is the manifest
inside that tarball; its digest is pinned here too.
"""

from phu.db.bundle import BundlePins

PRODUCTION_PINS = BundlePins(
    url="https://zenodo.org/records/24290671/files/phu-db.tar.gz",
    sha256="866153874ac37bfea9e88197f97baec9cd9b5b5b4fb1608b0aede60f047af641",
    size=2286259337,
    bundle_version="phu-db-0.3.0",
    manifest_sha256="fefa64508b1cea0d77e38604e67f6c9181dfc93fbe8dbaa1b8817377297c502d",
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
        }
    ),
)
