"""The phu-db bundle this phu release installs.

phu-db-0.1.1, Zenodo record 24121209 (https://doi.org/10.5281/zenodo.24121209),
built by phu-db-prep 0.1.1. The vendored phu/data/phu-db.json is the manifest
inside that tarball; its digest is pinned here too.
"""

from phu.db.bundle import BundlePins

PRODUCTION_PINS = BundlePins(
    url="https://zenodo.org/records/24121209/files/phu-db.tar.gz",
    sha256="7975c49ee78bb21b1e61b14c3f21554f3e2204a6fef58f5b039413152e6439b6",
    size=1852286645,
    bundle_version="phu-db-0.1.1",
    manifest_sha256="9fd6c63b0e77d7904cbc20e18af62620d3d79c686c37c64747d8d01818168be4",
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
        }
    ),
)
