"""The phu-db bundle this phu release installs.

phu-db-0.2.0, Zenodo record 24183949 (https://doi.org/10.5281/zenodo.24183949),
built by phu-db-prep 0.2.0. The vendored phu/data/phu-db.json is the manifest
inside that tarball; its digest is pinned here too.
"""

from phu.db.bundle import BundlePins

PRODUCTION_PINS = BundlePins(
    url="https://zenodo.org/records/24183949/files/phu-db.tar.gz",
    sha256="61db3c5b0840c37177adf05c76033715c8a7a709a7465a517349cd052704b5df",
    size=2242834752,
    bundle_version="phu-db-0.2.0",
    manifest_sha256="e735e55b2edb312194c589bdcdacd0fc0fe72a242996cdd4eb2022fe23374fd1",
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
            "phu-db/phrog/annotations.tsv",
        }
    ),
)
