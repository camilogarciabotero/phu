# Citation

`CITATION.cff` is the source of truth for PHU citation metadata. This
documentation describes the `0.10.0` release.

For a published study, cite the immutable PHU release used for the analysis
and record its version in the run manifest. Use the Zenodo concept DOI for
software-level references that are not tied to one release. The release DOI
must come from the corresponding release metadata and must not be inferred
from this development page.

PHU analyses should also cite the external tools and databases actually used,
including their versions or immutable releases. The relevant scientific
citations and database provenance belong in the run manifest alongside the PHU
version.

See the repository's [`CITATION.cff`](https://github.com/camilogarciabotero/phu/blob/main/CITATION.cff)
for the machine-readable author and software metadata.

For the databases, cite the `phu-db` bundle version recorded in `.phu/` and
the upstream databases you relied on; see
[Databases and Provenance](databases.md#citation).
