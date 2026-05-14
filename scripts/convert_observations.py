"""
convert_observations.py
-----------------------
Converts raw field observation CSVs (produced by geojson-to-csv.py) into the
EpiphytesDatabase format used by the R processing pipeline.

Each observation POINT in the GeoJSON may describe multiple individuals at
different heights. This script expands those into one row per individual.

Usage:
    python convert_observations.py --input ObservationData.csv --output ProcessedObservations.csv

Optional flags:
    --source        Value for the Source column        (default: "Field Observations")
    --site          Value for the Area_or_Site column  (default: "Unknown Site")
    --exp-start     Expedition start year              (default: current year)
    --exp-end       Expedition end year                (default: current year)

Output columns (in order):
    Source, Area_or_Site, lat, lon, Elevation_m, FieldID, Abundance,
    Height_m, CanopyHeight_m, pictures, note, AI_ID, FinalID, Genus, species

Description field formats handled (all field sites):
    Maquipucuna / Mashpi style (single-line, newlines collapsed by geojson-to-csv):
        "1292m Maxillariinae ?  6.4m  Canopy 12m Foto 96"
        "1351m Ind 5.1m foto 178 Ind 5.4m foto 179 Canopy 16.8m"
        "1326m Morfoespecie 1 fotos 109-113 2.5-3.5m Morfoespecie 2 fotos 114-117 3.9-4.2m Canopy 6.8m"

    MindoTarabita style (English / mixed):
        "Sudamelycaste maybe fimbriata  at 6.2m one big individual"
        "Xylobium sp1 At 3.1 m big individual At 5.5m smaller individual"

    MiradorMindo style (simplest):
        "1707m Maxillaria Flor blanca  6m"
"""

import argparse
import ast
import re
from datetime import datetime

import pandas as pd


# ── KNOWN TAXA ───────────────────────────────────────────────────────────────
KNOWN_TAXA = {
    "maxillaria", "maxillariinae", "brassea", "pleurothallis",
    "sudamelycaste", "sudamerlycaste", "xylobium", "morfoespecie",
    "orchidaceae", "anguloa", "acianthera", "ida", "mormolyca",
}

# ── MONTHS ───────────────────────────────────────────────────────────────

MONTH_ES = {
    "enero": 1, "febrero": 2, "marzo": 3, "abril": 4,
    "mayo": 5, "junio": 6, "julio": 7, "agosto": 8,
    "septiembre": 9, "octubre": 10, "noviembre": 11, "diciembre": 12
}

MONTH_EN = {
    "january": 1, "february": 2, "march": 3, "april": 4,
    "may": 5, "june": 6, "july": 7, "august": 8,
    "september": 9, "october": 10, "november": 11, "december": 12
}

# ── HELPERS ────────────────────────────────────────────────────────────────────

def _parse_datetime(name_str):
    """
    Parses datetime from the GeoJSON 'name' field.
    Handles two formats:
      Spanish: "27 de marzo de 2026, 11:52"
      English: "6 March 2026 at 14:04"
    Returns a datetime object or None.
    """
    s = str(name_str).strip()

    # Spanish: "27 de marzo de 2026, 11:52"
    m = re.match(
        r'(\d{1,2})\s+de\s+(\w+)\s+de\s+(\d{4}),\s+(\d{1,2}):(\d{2})',
        s, re.I
    )
    if m:
        day, month_str, year, hour, minute = m.groups()
        month = MONTH_ES.get(month_str.lower())
        if month:
            return datetime(int(year), month, int(day), int(hour), int(minute))

    # English: "6 March 2026 at 14:04"
    m = re.match(
        r'(\d{1,2})\s+(\w+)\s+(\d{4})\s+at\s+(\d{1,2}):(\d{2})',
        s, re.I
    )
    if m:
        day, month_str, year, hour, minute = m.groups()
        month = MONTH_EN.get(month_str.lower())
        if month:
            return datetime(int(year), month, int(day), int(hour), int(minute))

    return None

def _parse_elevation(line):
    """
    Extracts elevation from the start of a description.
    Must be a whole number >= 100 to distinguish from individual heights.
    """
    m = re.match(r'^(\d+)\s*m\b', line.strip(), re.I)
    if m:
        val = float(m.group(1))
        if val >= 100:
            return val
    return None


def _parse_canopy(line):
    """Extracts canopy/tree/arbol height from anywhere in a line."""
    m = re.search(r'(?:canopy|tree|arbol)\s*[~+]?\s*(\d+(?:\.\d+)?)\s*m', line, re.I)
    return float(m.group(1)) if m else None


def _extract_genus(line):
    """
    Finds the first known genus word in a line and returns it with any
    qualifier (sp., sp1, ?, maybe, cf.) that immediately follows.
    """
    for word in re.findall(r'[A-Za-z][a-z]+', line):
        if word.lower() in KNOWN_TAXA:
            idx = line.lower().find(word.lower())
            snippet = line[idx:].split()
            genus = snippet[0]
            extras = []
            for w in snippet[1:]:
                if re.match(r'^(sp\.?|sp\d+|\?|maybe|cf\.?)$', w, re.I):
                    extras.append(w)
                else:
                    break
            return (genus + " " + " ".join(extras)).strip() if extras else genus
    return None


def _is_genus_line(line):
    """True if the first word of the line is a known genus (multiline mode)."""
    words = line.strip().split()
    if not words:
        return False
    return words[0].lower().rstrip("?,.") in KNOWN_TAXA


def _parse_individuals_multiline(line):
    """
    Individual parser for multiline descriptions (one pattern per line).
    Used when geojson-to-csv preserves newlines (older GeoJSON exports).
    """
    inds = []
    line = line.strip()

    has_canopy = bool(re.search(r'(?:canopy|tree|arbol)\b', line, re.I))
    has_ind    = bool(re.search(
        r'(?:[Ii]nd\b|foto\s+\d|\bat\s+\d|\d+(?:\.\d+)?\s*m\s+foto)', line, re.I
    ))
    if has_canopy and not has_ind:
        return inds

    # 1. "Ind Xm [foto N]" / "x2 Ind Xm"
    for m in re.finditer(
        r'(?:x(\d+)\s+)?[Ii]nd(?:\s+\w+)?\s+(\d+(?:\.\d+)?)\s*m(?:\s+foto\s+([\d\-]+))?',
        line
    ):
        inds.append({'count': int(m.group(1)) if m.group(1) else 1,
                     'height_m': float(m.group(2)), 'photo': m.group(3), 'note': None})
    if inds: return inds

    # 2. "Foto(s) N a Xm"
    for m in re.finditer(
        r'[Ff]otos?\s+([\d\-]+)(?:\s*\((\d+)\s*ind\))?\s+a\s+(\d+(?:\.\d+)?)\s*m', line
    ):
        inds.append({'count': int(m.group(2)) if m.group(2) else 1,
                     'height_m': float(m.group(3)), 'photo': m.group(1), 'note': None})
    if inds: return inds

    # 3. "Xm foto N"
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s*m\s+foto\s+([\d\-]+)', line):
        inds.append({'count': 1, 'height_m': float(m.group(1)),
                     'photo': m.group(2), 'note': None})
    if inds: return inds

    # 4. "N ind foto N-N a Xm"
    for m in re.finditer(
        r'(\d+)\s+ind\s+foto\s+([\d\-]+)\s+a\s+(\d+(?:\.\d+)?)\s*m', line, re.I
    ):
        inds.append({'count': int(m.group(1)), 'height_m': float(m.group(3)),
                     'photo': m.group(2), 'note': None})
    if inds: return inds

    # 5. English "at Xm"
    for m in re.finditer(r'\bat\s+(\d+(?:\.\d+)?)\s*m\b', line, re.I):
        inds.append({'count': 1, 'height_m': float(m.group(1)),
                     'photo': None, 'note': line})
    if inds: return inds

    # 6. Height range "Xm - Ym"
    m = re.match(r'^(\d+(?:\.\d+)?)\s*m?\s*[-–]\s*(\d+(?:\.\d+)?)\s*m\s*$', line.strip())
    if m:
        h_min, h_max = float(m.group(1)), float(m.group(2))
        return [{'count': 1, 'height_m': round((h_min + h_max) / 2, 2),
                 'photo': None, 'note': f"range {h_min}-{h_max}m"}]

    # 7. Bare height
    m = re.match(r'^(\d+(?:\.\d+)?)\s*m\s*$', line.strip())
    if m:
        inds.append({'count': 1, 'height_m': float(m.group(1)),
                     'photo': None, 'note': None})
    return inds


def _parse_individuals_singleline(line, elevation=None, canopy=None):
    """
    Individual parser for collapsed single-line descriptions.
    Strips elevation prefix and canopy segment before matching,
    then tries patterns in order of specificity.
    """
    inds = []

    # Strip elevation prefix
    clean = re.sub(r'^\d+\s*m\b\s*', '', line, count=1) if elevation else line
    # Strip canopy segment
    clean = re.sub(
        r'(?:canopy|tree|arbol)\s*[~+]?\s*\d+(?:\.\d+)?\s*m\b[^,\n]*?(?=\s+[A-Z]|\s*$)',
        '', clean, flags=re.I
    ).strip()

    # P1: "Ind Xm [foto N]" / "x2 Ind Xm"
    for m in re.finditer(
        r'(?:x(\d+)\s+)?[Ii]nd(?:\s+\w+)?\s+(\d+(?:\.\d+)?)\s*m(?:\s+foto\s+([\d\-]+))?',
        clean
    ):
        inds.append({'count': int(m.group(1)) if m.group(1) else 1,
                     'height_m': float(m.group(2)), 'photo': m.group(3), 'note': None})
    if inds: return inds

    # P2: "foto N a Xm" / "fotos N-N (X ind) a Xm"
    for m in re.finditer(
        r'[Ff]otos?\s+([\d\-]+)(?:\s*\((\d+)\s*ind\))?\s+a\s+(\d+(?:\.\d+)?)\s*m',
        clean
    ):
        inds.append({'count': int(m.group(2)) if m.group(2) else 1,
                     'height_m': float(m.group(3)), 'photo': m.group(1), 'note': None})
    if inds: return inds

    # P3: "Xm foto N"
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s*m\s+foto\s+([\d\-]+)', clean):
        inds.append({'count': 1, 'height_m': float(m.group(1)),
                     'photo': m.group(2), 'note': None})
    if inds: return inds

    # P4: "N ind foto N-N a Xm"
    for m in re.finditer(
        r'(\d+)\s+ind\s+foto\s+([\d\-]+)\s+a\s+(\d+(?:\.\d+)?)\s*m', clean, re.I
    ):
        inds.append({'count': int(m.group(1)), 'height_m': float(m.group(3)),
                     'photo': m.group(2), 'note': None})
    if inds: return inds

    # P5: "foto N Xm a Ym" (range after photo, no 'a' before range)
    for m in re.finditer(
        r'[Ff]otos?\s+([\d\-]+)(?:\s+\w+)?\s+(\d+(?:\.\d+)?)\s*m\s+a\s+(\d+(?:\.\d+)?)\s*m',
        clean
    ):
        h_min, h_max = float(m.group(2)), float(m.group(3))
        inds.append({'count': 1, 'height_m': round((h_min + h_max) / 2, 2),
                     'photo': m.group(1), 'note': f"range {h_min}-{h_max}m"})
    if inds: return inds

    # P6: "fotos N-N Xm-Ym" or "fotos N-N [word] Xm-Ym" (photo range + height range)
    for m in re.finditer(
        r'[Ff]otos?\s+([\d\-]+)(?:\s+\w+)?\s+(\d+(?:\.\d+)?)\s*m?\s*[-–]\s*(\d+(?:\.\d+)?)\s*m',
        clean
    ):
        h_min, h_max = float(m.group(2)), float(m.group(3))
        inds.append({'count': 1, 'height_m': round((h_min + h_max) / 2, 2),
                     'photo': m.group(1), 'note': f"range {h_min}-{h_max}m"})
    if inds: return inds

    # P7: "foto N Xm Ym ..." (photo then bare heights)
    for m in re.finditer(
        r'[Ff]otos?\s+([\d\-]+)\s+((?:\d+(?:\.\d+)?\s*m\s*)+)', clean
    ):
        photo = m.group(1)
        for hm in re.finditer(r'(\d+(?:\.\d+)?)\s*m', m.group(2)):
            val = float(hm.group(1))
            if val < 100:
                inds.append({'count': 1, 'height_m': val,
                             'photo': photo, 'note': None})
    if inds: return inds

    # P8: English "at Xm"
    for m in re.finditer(r'\bat\s+(\d+(?:\.\d+)?)\s*m\b', clean, re.I):
        inds.append({'count': 1, 'height_m': float(m.group(1)),
                     'photo': None, 'note': clean.strip()})
    if inds: return inds

    # P9: Height range "Xm - Ym" embedded in line
    m = re.search(r'(\d+(?:\.\d+)?)\s*m?\s*[-–]\s*(\d+(?:\.\d+)?)\s*m\b', clean)
    if m:
        h_min, h_max = float(m.group(1)), float(m.group(2))
        if h_min < 100 and h_max < 100:
            inds.append({'count': 1, 'height_m': round((h_min + h_max) / 2, 2),
                         'photo': None, 'note': f"range {h_min}-{h_max}m"})
            return inds

    # P10: Any bare Xm < 100 embedded in line
    for m in re.finditer(r'(?<!\d)(\d+(?:\.\d+)?)\s*m\b', clean):
        val = float(m.group(1))
        if val < 100:
            inds.append({'count': 1, 'height_m': val, 'photo': None, 'note': None})

    return inds


# ── DESCRIPTION PARSER ─────────────────────────────────────────────────────────

def parse_description(desc):
    """
    Parses a full field-note description into structured fields.
    Handles both multiline GeoJSON and single-line collapsed versions.

    Returns a dict:
        elevation_m  : float or None
        genus        : str or None
        canopy_m     : float or None
        individuals  : list of {height_m, count, photo, note}
    """
    if not desc or not desc.strip():
        return None

    result = {'elevation_m': None, 'genus': None, 'canopy_m': None, 'individuals': []}

    lines = [l.strip() for l in desc.strip().splitlines() if l.strip()]

    if len(lines) > 1:
        # ── Multiline ─────────────────────────────────────────────────────
        for i, line in enumerate(lines):
            if i == 0:
                elev = _parse_elevation(line)
                if elev is not None:
                    result['elevation_m'] = elev
                    continue
            canopy = _parse_canopy(line)
            if canopy is not None:
                result['canopy_m'] = canopy
                continue
            if re.match(r'^canopy\s*$', line, re.I):
                continue
            if _is_genus_line(line) and result['genus'] is None:
                if not re.search(r'\d+(?:\.\d+)?\s*m\b', line):
                    result['genus'] = line
                    continue
            result['individuals'].extend(_parse_individuals_multiline(line))

    else:
        # ── Single line (collapsed by geojson-to-csv) ──────────────────────
        line = desc.strip()

        result['elevation_m'] = _parse_elevation(line)
        result['canopy_m']    = _parse_canopy(line)
        result['genus']       = _extract_genus(line)
        result['individuals'] = _parse_individuals_singleline(
            line,
            elevation = result['elevation_m'],
            canopy    = result['canopy_m']
        )

    return result


# ── COORDINATE PARSER ──────────────────────────────────────────────────────────

def parse_coords(coord_str):
    """Converts GeoJSON '[lon, lat]' string to (lon, lat) floats."""
    coords = ast.literal_eval(coord_str)
    return float(coords[0]), float(coords[1])


# ── MAIN CONVERSION ────────────────────────────────────────────────────────────

# Column order for the output CSV — matches the agreed schema exactly.
OUTPUT_COLUMNS = [
    'Source', 'Area_or_Site', 'lat', 'lon', 'Elevation_m',
    'FieldID', 'Abundance', 'Height_m', 'CanopyHeight_m',
    'pictures', 'note', 'AI_ID', 'FinalID', 'Genus', 'species', 'datetime'
]


def convert(input_csv, output_csv, source, site):
    """
    Reads the raw GeoJSON-derived CSV, parses every observation point, and
    writes one output row per individual (expanding multi-individual points).
    """
    df_raw = pd.read_csv(input_csv)

    obs = df_raw[
        (df_raw['g__type'] == 'Point') & df_raw['description'].notna()
    ].copy()

    print(f"Input rows:                       {len(df_raw)}")
    print(f"Observation points with desc:     {len(obs)}")
    print(f"Skipped (tracks / empty markers): {len(df_raw) - len(obs)}")

    rows = []
    for _, obs_row in obs.iterrows():
        parsed = parse_description(obs_row['description'])
        if parsed is None:
            continue

        lon, lat = parse_coords(obs_row['g__coordinates'])
        obs_dt = _parse_datetime(obs_row.get('name', ''))

        base = {
            'Source':         source,
            'Area_or_Site':   site,
            'lat':            lat,
            'lon':            lon,
            'Elevation_m':    parsed['elevation_m'],
            'FieldID':        None,   # left blank — fill in post-hoc if needed
            'Abundance':      None,
            'Height_m':       None,
            'CanopyHeight_m': parsed['canopy_m'],
            'pictures':       None,
            'note':           None,
            'AI_ID':          None,   # filled during identification workflow
            'FinalID':        None,   # filled after expert review
            'Genus':          parsed['genus'],
            'species':        None,
            'datetime':       obs_dt.strftime("%Y-%m-%d %H:%M:%S") if obs_dt else None,
        }

        if parsed['individuals']:
            for ind in parsed['individuals']:
                r = base.copy()
                r['Height_m']  = ind['height_m']
                r['Abundance'] = ind['count']
                r['pictures']  = ind.get('photo')
                r['note']      = ind.get('note') or (
                    f"foto {ind['photo']}" if ind.get('photo') else None
                )
                rows.append(r)
        else:
            rows.append(base)

    df_out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    # Ensure pictures is always stored as string (prevents type clash in R bind_rows)
    df_out['pictures'] = df_out['pictures'].astype(str).where(df_out['pictures'].notna())

    df_out.to_csv(output_csv, index=False, na_rep='')
    print(f"\nSaved {len(df_out)} rows -> {output_csv}")
    print(f"  With Elevation_m:    {df_out['Elevation_m'].notna().sum()}")
    print(f"  With Height_m:       {df_out['Height_m'].notna().sum()}")
    print(f"  With CanopyHeight_m: {df_out['CanopyHeight_m'].notna().sum()}")
    print(f"  With Genus:          {df_out['Genus'].notna().sum()}")
    print(f"  With pictures:       {df_out['pictures'].notna().sum()}")
    return df_out


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    current_year = datetime.now().year

    parser = argparse.ArgumentParser(
        description="Convert raw field observation CSV to EpiphytesDatabase format."
    )
    parser.add_argument("--input",     required=True,
                        help="Path to input CSV (from geojson-to-csv.py)")
    parser.add_argument("--output",    required=True,
                        help="Path for output CSV")
    parser.add_argument("--source",    default="Field Observations")
    parser.add_argument("--site",      default="Unknown Site")

    args = parser.parse_args()

    convert(
        input_csv  = args.input,
        output_csv = args.output,
        source     = args.source,
        site       = args.site
    )


if __name__ == "__main__":
    main()