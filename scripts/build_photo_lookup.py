"""
build_photo_lookup.py
---------------------
Builds a photo identification lookup table from the raw GeoJSON-derived CSVs.

Works directly from the raw CSVs (output of geojson-to-csv.py) rather than
the processed EpiphytesDatabase CSVs, because the raw descriptions contain
all photo references in their original form.

Produces photo_lookup.csv with one row per photo, columns aligned to the
combined.csv schema. AI_ID, FinalID, Genus, species, and note are left blank
for you to fill in during the identification workflow.

Usage:
    # Single site:
    python build_photo_lookup.py \\
        --input  geojson_to_csv/csv/Maquipucuna.csv \\
        --sites  Maquipucuna \\
        --output data/photo_lookup.csv

    # All sites at once:
    python build_photo_lookup.py \\
        --input  geojson_to_csv/csv/Maquipucuna.csv \\
                 geojson_to_csv/csv/Mashpi.csv \\
                 geojson_to_csv/csv/MindoTarabita.csv \\
        --sites  Maquipucuna Mashpi MindoTarabita \\
        --output data/photo_lookup.csv \\
        --sd-card /Volumes/YOUR_SD_CARD/

Optional flags:
    --offset    Offset from field note number to camera file number (default: 9000)
                field note "foto 96" + 9000 = _MG_9096.CR2
    --sd-card   Path to SD card folder — checks which files actually exist
"""

import argparse
import ast
import os
import re

import pandas as pd


# ── FOTO REFERENCE EXTRACTOR ───────────────────────────────────────────────────

def extract_foto_refs(desc):
    """
    Extracts all photo numbers from a raw description string.
    Handles all field note formats found across sites:

        "Foto 96"                      -> [96]
        "fotos 109-113"                -> [109,110,111,112,113]
        "133-135 fotos"                -> [133,134,135]
        "152-156 Maxillaria"           -> [152,153,154,155,156]  (bare range)
        "Foto 139 a 2m Foto 140 a 1m"  -> [139,140]
        "Ind 5.1m foto 178"            -> [178]
        "foto 283-285 a 9.1m"          -> [283,284,285]

    Returns a list of ints in order of appearance (deduplicated).
    """
    if not desc or not isinstance(desc, str):
        return []

    refs = []

    # Pattern 1: "foto/fotos N[-N]" — number NOT followed by decimal point
    for m in re.finditer(r'fotos?\s+(\d+)(?:-(\d+))?(?!\.\d)', desc, re.I):
        start = int(m.group(1))
        end   = int(m.group(2)) if m.group(2) else start
        refs.extend(range(start, end + 1))

    # Pattern 2: "N-N fotos?" — range BEFORE the keyword
    for m in re.finditer(r'(\d+)-(\d+)\s+fotos?', desc, re.I):
        refs.extend(range(int(m.group(1)), int(m.group(2)) + 1))

    # Pattern 3: standalone number ranges with no 'm' suffix, in photo-number
    # range (50-500). Catches "152-156 Maxillaria" style bare photo ranges.
    for m in re.finditer(r'(?<!\d)(\d{2,4})-(\d{2,4})(?!\s*m\b)(?!\.\d)', desc, re.I):
        s, e = int(m.group(1)), int(m.group(2))
        if 50 <= s <= 500 and 50 <= e <= 500 and e > s:
            refs.extend(range(s, e + 1))

    # Deduplicate preserving order of first appearance
    seen, result = set(), []
    for r in refs:
        if r not in seen:
            seen.add(r)
            result.append(r)
    return result


def foto_to_filename(foto_number, offset=9000):
    """Converts a field note foto number to the camera filename."""
    return f"_MG_{foto_number + offset}.CR2"


def parse_coords(coord_str):
    """Converts GeoJSON '[lon, lat]' string to (lon, lat) floats."""
    try:
        coords = ast.literal_eval(coord_str)
        return float(coords[0]), float(coords[1])
    except Exception:
        return None, None


# ── MAIN ───────────────────────────────────────────────────────────────────────

# Column order matches combined.csv schema exactly.
# FieldID is omitted here — it is not derivable from raw GeoJSON CSVs.
OUTPUT_COLUMNS = [
    'Source', 'Area_or_Site', 'lat', 'lon',
    'foto_number', 'filename', 'file_exists',
    'pictures', 'note',
    'AI_ID', 'FinalID', 'Genus', 'species',
    'description',   # kept for reference during ID workflow
]


def build_lookup(input_csvs, site_names, output_csv, offset=9000, sd_card_path=None):

    all_rows = []

    for csv_path, site_name in zip(input_csvs, site_names):
        df = pd.read_csv(csv_path)

        # Only Point geometries with a description
        obs = df[(df['g__type'] == 'Point') & df['description'].notna()].copy()
        print(f"\n{site_name}: {len(obs)} observation points")

        for _, row in obs.iterrows():
            desc     = row['description']
            lon, lat = parse_coords(row['g__coordinates'])
            refs     = extract_foto_refs(desc)

            if refs:
                for foto_num in refs:
                    filename = foto_to_filename(foto_num, offset)
                    if sd_card_path:
                        file_exists = os.path.exists(os.path.join(sd_card_path, filename))
                    else:
                        file_exists = None

                    all_rows.append({
                        'Source':       'Field Observations',
                        'Area_or_Site': site_name,
                        'lat':          lat,
                        'lon':          lon,
                        'foto_number':  foto_num,
                        'filename':     filename,
                        'file_exists':  file_exists,
                        'pictures':     str(foto_num),  # mirrors combined.csv pictures column
                        'note':         None,
                        'AI_ID':        None,   # fill during identification workflow
                        'FinalID':      None,   # fill after expert review
                        'Genus':        None,   # fill during identification workflow
                        'species':      None,   # fill during identification workflow
                        'description':  desc,
                    })
            else:
                # No photo reference — include so nothing gets lost
                all_rows.append({
                    'Source':       'Field Observations',
                    'Area_or_Site': site_name,
                    'lat':          lat,
                    'lon':          lon,
                    'foto_number':  None,
                    'filename':     None,
                    'file_exists':  None,
                    'pictures':     None,
                    'note':         None,
                    'AI_ID':        None,
                    'FinalID':      None,
                    'Genus':        None,
                    'species':      None,
                    'description':  desc,
                })

    df_out = pd.DataFrame(all_rows, columns=OUTPUT_COLUMNS)
    df_out.to_csv(output_csv, index=False)

    total      = len(df_out)
    with_photo = df_out['foto_number'].notna().sum()
    no_photo   = df_out['foto_number'].isna().sum()

    print(f"\n{'='*50}")
    print(f"Total rows:            {total}")
    print(f"  With photo ref:      {with_photo}")
    print(f"  Without photo ref:   {no_photo}")
    if sd_card_path:
        found   = (df_out['file_exists'] == True).sum()
        missing = (df_out['file_exists'] == False).sum()
        print(f"  Files on SD found:   {found}")
        print(f"  Files on SD missing: {missing}")
    print(f"\nSaved → {output_csv}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build photo lookup table from raw GeoJSON-derived CSVs."
    )
    parser.add_argument("--input",    required=True, nargs="+",
                        help="One or more raw CSV files (from geojson-to-csv.py)")
    parser.add_argument("--sites",    required=True, nargs="+",
                        help="Site name for each input file (same order)")
    parser.add_argument("--output",   required=True,
                        help="Path for output photo_lookup.csv")
    parser.add_argument("--offset",   default=9000, type=int,
                        help="Offset from field note number to file number (default: 9000)")
    parser.add_argument("--sd-card",  default=None,
                        help="Path to SD card folder to verify files exist")

    args = parser.parse_args()

    if len(args.input) != len(args.sites):
        parser.error("--input and --sites must have the same number of entries")

    build_lookup(
        input_csvs   = args.input,
        site_names   = args.sites,
        output_csv   = args.output,
        offset       = args.offset,
        sd_card_path = args.sd_card,
    )


if __name__ == "__main__":
    main()