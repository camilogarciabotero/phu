from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union
from enum import Enum

import pandas as pd

class OutputFormat(str, Enum):
    csv = "csv"
    tsv = "tsv"


@dataclass
class TaxaConfig:
    """Configuration for taxa simplification."""
    input_file: Path
    output_file: Path
    add_lineage: bool = False
    lineage_col: str = "compact_lineage"
    sep: Optional[str] = None
    
    def plan(self) -> "TaxaPlan":
        """Create execution plan from configuration."""
        # Auto-detect separator if not provided
        separator = self.sep
        if separator is None:
            if str(self.input_file).lower().endswith(".tsv"):
                separator = "\t"
            else:
                separator = ","
        
        # Determine output format
        output_path_lower = str(self.output_file).lower()
        if output_path_lower.endswith(".tsv"):
            output_format = OutputFormat.tsv
        else:
            output_format = OutputFormat.csv
            
        return TaxaPlan(
            input_file=self.input_file,
            output_file=self.output_file,
            add_lineage=self.add_lineage,
            lineage_col=self.lineage_col,
            separator=separator,
            output_format=output_format,
        )


@dataclass
class TaxaPlan:
    """Execution plan for taxa simplification."""
    input_file: Path
    output_file: Path
    add_lineage: bool
    lineage_col: str
    separator: str
    output_format: OutputFormat
# --- Pattern helpers ---------------------------------------------------------

# Capture for the known ancestor/anchor name:
# anything except ":" or "|" to avoid colliding with our output separators.
ANCHOR = r"([^:\|]+)"

# Compile order: most specific → most general (first match wins).
TAXA_PATTERNS: Dict[str, List[Dict[str, re.Pattern]]] = {
    "kingdom": [
        {
            "pattern": re.compile(rf"novel_kingdom_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",  # placeholder, handled below (we use callable repl)
        },
    ],
    "phylum": [
        {
            "pattern": re.compile(rf"novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_phylum_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
    "class": [
        # Edge cases for class 0 chains - MUST be first to match before regular patterns
        {
            "pattern": re.compile(rf"novel_class_0_of_novel_phylum_0_of_novel_kingdom_(\d+)_of_{ANCHOR}"),
            "replacement": r"EDGECASE_WITH_ANCHOR",
        },
        {
            "pattern": re.compile(r"novel_class_0_of_novel_phylum_0_of_novel_kingdom"),
            "replacement": r"EDGECASE_NK0:NP0:NC0",
        },
        # Regular class patterns
        {
            "pattern": re.compile(rf"novel_class_(\d+)_of_novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_class_(\d+)_of_novel_phylum_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_class_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
    "order": [
        # Edge cases for order 0 chains - MUST be first to match before regular patterns
        {
            "pattern": re.compile(rf"novel_order_0_of_novel_class_0_of_novel_phylum_0_of_novel_kingdom_(\d+)_of_{ANCHOR}"),
            "replacement": r"EDGECASE_ORDER_WITH_ANCHOR",
        },
        {
            "pattern": re.compile(r"novel_order_0_of_novel_class_0_of_novel_phylum_0_of_novel_kingdom"),
            "replacement": r"EDGECASE_ORDER_NK0:NP0:NC0:NO0",
        },
        # Regular order patterns
        {
            "pattern": re.compile(
                rf"novel_order_(\d+)_of_novel_family_(\d+)_of_novel_class_(\d+)_of_novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"
            ),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_order_(\d+)_of_novel_family_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_order_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
    "family": [
        {
            "pattern": re.compile(
                rf"novel_family_(\d+)_of_novel_order_(\d+)_of_novel_class_(\d+)_of_novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"
            ),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_novel_order_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
    "subfamily": [
        {
            "pattern": re.compile(
                rf"novel_subfamily_(\d+)_of_novel_family_(\d+)_of_novel_order_(\d+)_of_novel_class_(\d+)_of_novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"
            ),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_subfamily_(\d+)_of_novel_family_(\d+)_of_novel_order_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_subfamily_(\d+)_of_novel_family_(\d+)$"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_subfamily_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        # helper fallbacks mirroring R
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_novel_order_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
    "genus": [
        {
            "pattern": re.compile(
                rf"novel_genus_(\d+)_of_novel_subfamily_(\d+)_of_novel_family_(\d+)_of_novel_order_(\d+)_of_novel_class_(\d+)_of_novel_phylum_(\d+)_of_novel_kingdom_(\d+)_of_{ANCHOR}"
            ),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(
                rf"novel_genus_(\d+)_of_novel_subfamily_(\d+)_of_novel_family_(\d+)_of_novel_order_(\d+)_of_{ANCHOR}"
            ),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_genus_(\d+)_of_novel_subfamily_(\d+)_of_novel_family_(\d+)$"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_genus_(\d+)_of_novel_subfamily_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_genus_(\d+)_of_novel_subfamily_(\d+)$"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_genus_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        # helper fallbacks
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_novel_order_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
        {
            "pattern": re.compile(rf"novel_family_(\d+)_of_{ANCHOR}"),
            "replacement": r"\1",
        },
    ],
}

# rank code tags
RANK_CODE = {
    "kingdom": "NK",
    "phylum": "NP",
    "class": "NC",
    "order": "NO",
    "family": "NF",
    "subfamily": "NSF",
    "genus": "NG",
}

def _apply_first_match(x: str, level: str) -> str:
    """Apply the first matching regex for a single taxon string at a given level."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    # Edgecase token direct return
    for rule in TAXA_PATTERNS.get(level, []):
        m = rule["pattern"].search(x)
        if not m:
            continue

        # Special hard-coded edge cases for 0 chains
        if m and "EDGECASE_NK0:NP0:NC0" in rule["replacement"]:
            return "NK0:NP0:NC0"
        elif m and "EDGECASE_WITH_ANCHOR" in rule["replacement"]:
            # Handle novel_class_0_of_novel_phylum_0_of_novel_kingdom_X_of_Anchor
            groups = m.groups()
            kingdom_id = groups[0]  # First captured group is kingdom ID
            anchor = groups[1]      # Second captured group is anchor
            return f"{anchor}:NK{kingdom_id}:NP0:NC0"
        elif m and "EDGECASE_ORDER_NK0:NP0:NC0:NO0" in rule["replacement"]:
            return "NK0:NP0:NC0:NO0"
        elif m and "EDGECASE_ORDER_WITH_ANCHOR" in rule["replacement"]:
            # Handle novel_order_0_of_novel_class_0_of_novel_phylum_0_of_novel_kingdom_X_of_Anchor
            groups = m.groups()
            kingdom_id = groups[0]  # First captured group is kingdom ID
            anchor = groups[1]      # Second captured group is anchor
            return f"{anchor}:NK{kingdom_id}:NP0:NC0:NO0"

        # We reconstruct the output using captured groups explicitly,
        # because we need both IDs and the anchor name, and the IDs
        # can be in multiple positions depending on the pattern.
        groups = m.groups()
        # Find the anchor: last group matching ANCHOR; IDs are digits.
        # We’ll pick the last non-digit group as anchor.
        anchor = None
        ids = []
        for g in groups:
            if g is None:
                continue
            if g.isdigit():
                ids.append(g)
            else:
                anchor = g

        code = []
        if anchor:
            code.append(anchor)

        # Map IDs in the order appropriate to each level’s pattern.
        # We don’t know which ID corresponds to which rank purely from regex
        # position across all patterns, so we infer from level and the
        # count of IDs captured. The patterns are ordered to keep this safe.
        tag = RANK_CODE[level]

        # For mixed patterns above, the *last* captured number is the target level’s ID.
        # Upstream lineage numbers come before it. We emit upstream tags when we
        # can infer them (based on count) in descending rank order.
        # This is a pragmatic, pattern-order-dependent reconstruction mirroring your R replacements.

        # Helper to append tag:value
        def add(tagcode: str, val: str):
            code.append(f"{tagcode}{val}")

        if level == "kingdom":
            # anchor:NK<id>
            if ids:
                add("NK", ids[0])

        elif level == "phylum":
            # anchor:NK<k>:NP<p>  or   anchor:NP<p>
            if len(ids) == 2:
                add("NK", ids[1])  # kingdom id
                add("NP", ids[0])  # phylum id
            elif len(ids) == 1:
                add("NP", ids[0])

        elif level == "class":
            # edge case handled earlier
            # anchor:NK<k>:NP<p>:NC<c>   or   anchor:NP<p>:NC<c>   or anchor:NC<c>
            if len(ids) == 3:
                add("NK", ids[2])
                add("NP", ids[1])
                add("NC", ids[0])
            elif len(ids) == 2:
                add("NP", ids[1])
                add("NC", ids[0])
            elif len(ids) == 1:
                add("NC", ids[0])

        elif level == "order":
            # anchor:NK<k>:NP<p>:NC<c>:NF<f>:NO<o>  (full)
            # or anchor:NF<f>:NO<o> ; or anchor:NO<o>
            if len(ids) == 5:
                add("NK", ids[4])
                add("NP", ids[3])
                add("NC", ids[2])
                add("NF", ids[1])
                add("NO", ids[0])
            elif len(ids) == 2:
                add("NF", ids[1])
                add("NO", ids[0])
            elif len(ids) == 1:
                add("NO", ids[0])

        elif level == "family":
            # anchor:NK<k>:NP<p>:NC<c>:NO<o>:NF<f>
            # or anchor:NO<o>:NF<f> ; or anchor:NF<f>
            if len(ids) == 5:
                add("NK", ids[4])
                add("NP", ids[3])
                add("NC", ids[2])
                add("NO", ids[1])
                add("NF", ids[0])
            elif len(ids) == 2:
                add("NO", ids[1])
                add("NF", ids[0])
            elif len(ids) == 1:
                add("NF", ids[0])

        elif level == "subfamily":
            # anchor:NK<k>:NP<p>:NC<c>:NO<o>:NF<f>:NSF<sf>
            # or anchor:NO<o>:NF<f>:NSF<sf>
            # or NF<f>:NSF<sf> (no anchor)
            # or anchor:NSF<sf>
            # or family fallbacks
            if len(ids) == 6 and anchor:
                add("NK", ids[5])
                add("NP", ids[4])
                add("NC", ids[3])
                add("NO", ids[2])
                add("NF", ids[1])
                add("NSF", ids[0])
            elif len(ids) == 3 and anchor:
                add("NO", ids[2])
                add("NF", ids[1])
                add("NSF", ids[0])
            elif len(ids) == 2 and not anchor:
                # pattern: novel_subfamily_<sf>_of_novel_family_<f>
                code = [f"NF{ids[1]}", f"NSF{ids[0]}"]
            elif len(ids) == 1 and anchor:
                add("NSF", ids[0])
            elif len(ids) == 2 and anchor:
                # fallback family patterns captured under subfamily group
                add("NO", ids[1])  # from fallback pattern
                add("NF", ids[0])

        elif level == "genus":
            # anchor:NK<k>:NP<p>:NC<c>:NO<o>:NF<f>:NSF<sf>:NG<g>
            # or anchor:NO<o>:NF<f>:NSF<sf>:NG<g>
            # or NF<f>:NSF<sf>:NG<g> (no anchor)
            # or anchor:NSF<sf>:NG<g>
            # or NSF<sf>:NG<g> (no anchor)
            # or anchor:NG<g>
            # plus family fallbacks
            if len(ids) == 7 and anchor:
                add("NK", ids[6])
                add("NP", ids[5])
                add("NC", ids[4])
                add("NO", ids[3])
                add("NF", ids[2])
                add("NSF", ids[1])
                add("NG", ids[0])
            elif len(ids) == 4 and anchor:
                add("NO", ids[3])
                add("NF", ids[2])
                add("NSF", ids[1])
                add("NG", ids[0])
            elif len(ids) == 3 and not anchor:
                code = [f"NF{ids[2]}", f"NSF{ids[1]}", f"NG{ids[0]}"]
            elif len(ids) == 2 and anchor:
                add("NSF", ids[1])
                add("NG", ids[0])
            elif len(ids) == 2 and not anchor:
                code = [f"NSF{ids[1]}", f"NG{ids[0]}"]
            elif len(ids) == 1 and anchor:
                add("NG", ids[0])
            elif len(ids) == 2 and anchor:
                # fallback family patterns (under genus rules)
                add("NO", ids[1])
                add("NF", ids[0])

        # Join
        return ":".join(code) if code else x

    return x  # no match for this level


def simplify_single_taxon(x: Optional[str], level: str) -> Optional[str]:
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    # Split by || (multiple candidates)
    parts = str(x).split("||")
    out_parts = []
    for p in parts:
        p = p.strip()
        out_parts.append(_apply_first_match(p, level))
    return "||".join("" if v is None else v for v in out_parts)


def simplify_series(series: pd.Series, level: str) -> pd.Series:
    return series.astype("string").map(lambda v: simplify_single_taxon(v, level))


# --- Main execution function ------------------------------------------------

PRED_LEVELS = [
    "realm",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "subfamily",
    "genus",
]


def _simplify_taxa(cfg: TaxaConfig) -> TaxaPlan:
    """
    Simplify vContact taxonomy predictions into compact lineage codes.
    """
    if not cfg.input_file.exists():
        raise FileNotFoundError(f"Input file not found: {cfg.input_file}")

    plan = cfg.plan()
    
    print(f"Reading taxonomy data from {plan.input_file}...")
    
    try:
        df = pd.read_csv(
            plan.input_file, 
            sep=plan.separator, 
            dtype=str, 
            na_values=["", "NA", "NaN"]
        )
    except Exception as e:
        raise RuntimeError(f"Error reading {plan.input_file}: {e}")

    # Clean column names
    df.columns = _clean_cols(list(df.columns))
    level_map = _detect_level_map(df)

    if not level_map:
        print("No *_prediction columns found. Writing unchanged data...")
        _write_df(df, plan.output_file, plan.output_format)
        return plan

    print(f"Found {len(level_map)} taxonomy prediction columns to simplify...")

    # Apply simplification
    for col, lvl in level_map.items():
        print(f"Simplifying {col} ({lvl} level)...")
        df[col] = _simplify_series(df[col], lvl)

    # Optional lineage column from deepest available rank
    if plan.add_lineage:
        print(f"Adding {plan.lineage_col} column...")
        ordered = [c for c in [ 
            "genus_prediction", "subfamily_prediction", "family_prediction", "order_prediction",
            "class_prediction", "phylum_prediction", "kingdom_prediction", "realm_prediction" 
        ] if c in df.columns]
        
        def deepest(row):
            for c in ordered:
                val = row.get(c)
                if isinstance(val, str) and val.strip():
                    return val
            return pd.NA
            
        df[plan.lineage_col] = df.apply(deepest, axis=1)

    print(f"Writing results to {plan.output_file}...")
    _write_df(df, plan.output_file, plan.output_format)

    # QA summary
    print("\nQA Summary:")
    for col in level_map:
        remaining_novel = df[col].astype("string").str.contains(r"novel_", na=False).sum()
        print(f"  {col}: {remaining_novel} remaining 'novel_' strings")

    return plan


def _clean_cols(cols: List[str]) -> List[str]:
    """Basic snake_case & trimming; avoids janitor dependency."""
    out = []
    for c in cols:
        c2 = c.strip().replace(" ", "_").replace("-", "_")
        c2 = re.sub(r"__+", "_", c2)
        out.append(c2.lower())
    return out


def _detect_level_map(df: pd.DataFrame) -> Dict[str, str]:
    """Detect *_prediction columns and map them to taxonomic levels."""
    want = {f"{lvl}_prediction": lvl for lvl in PRED_LEVELS}
    return {col: lvl for col, lvl in want.items() if col in df.columns}


def _simplify_series(series: pd.Series, level: str) -> pd.Series:
    """Apply simplification to a pandas Series for a given taxonomic level."""
    return series.astype("string").map(lambda v: _simplify_single_taxon(v, level))


def _simplify_single_taxon(x: Optional[str], level: str) -> Optional[str]:
    """Simplify a single taxon string, handling multiple candidates separated by ||."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    # Split by || (multiple candidates)
    parts = str(x).split("||")
    out_parts = []
    for p in parts:
        p = p.strip()
        out_parts.append(_apply_first_match(p, level))
    return "||".join("" if v is None else v for v in out_parts)


def _write_df(df: pd.DataFrame, path: Path, output_format: OutputFormat) -> None:
    """Write dataframe to file in the specified format."""
    if output_format == OutputFormat.tsv:
        df.to_csv(path, sep="\t", index=False)
    else:  # csv
        df.to_csv(path, index=False)