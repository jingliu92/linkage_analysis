# espK/espV/espN–Serotype Linkage Analysis
⭐ Rationale
The espK, espV, and espN effector genes are strongly associated with enterohemorrhagic E. coli (EHEC), with approximately 99% of EHEC genomes carrying at least one of these markers. However, these esp genes are also present in a subset of enteropathogenic E. coli (EPEC) strains, which confounds the use of esp genotypes as simple, stand-alone indicators for distinguishing EHEC from EPEC. To clarify the lineage specificity of these genes, it is essential to determine whether esp-positive genotypes occur preferentially within particular EPEC serotypes and whether esp-negative EPEC strains are restricted to distinct serotype lineages.

By resolving the serotype distribution of esp-positive and esp-negative EPEC isolates, we can assess whether specific serotypes reliably predict esp gene carriage. If such linkage is observed, serotype markers may be used to identify esp-positive EPEC pathovars with higher accuracy, improving virulence prediction and strain classification in surveillance and genomic diagnostics.

To accomplish this, we performed a three-step analysis:

- Classify the distribution of espK/espV/espN within EHEC and EPEC, defining esp-positive and esp-negative subpopulations.

- Extract genome assemblies corresponding to each genotype category (esp+ EHEC, esp+ EPEC, esp– EPEC).

- Perform high-resolution serotyping using ECTyper to determine whether specific serotypes are enriched for esp-positive or esp-negative profiles.
# Step1: Classify the distribution of espK/espV/espN within EHEC and EPEC, defining esp-positive and esp-negative subpopulations.
## Create marker gene file
nano all.gene.fasta
## make db
```
makeblastdb -in all_markers.fasta -dbtype nucl -out all_markers
```
## Blast EHEC
```
mkdir -p blast_out_EHEC
mkdir -p blast_hits_only_EHEC

for f in $(find /home/jing/E.coli/blast_results/EHEC_assemblies -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out_EHEC/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query all_markers.fasta \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=90 && cov>=90) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only_EHEC/
      echo "✅ Hits found in $folder — copied to blast_hits_only_EHEC/"
  else
      echo "❌ No hits for $folder"
  fi
done
```

## Blast EPEC
```
mkdir -p blast_out_EPEC
mkdir -p blast_hits_only_EPEC

for f in $(find /home/jing/E.coli/blast_results/EPEC_assemblies -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out_EPEC/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query all_markers.fasta \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=90 && cov>=90) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only_EPEC/
      echo "✅ Hits found in $folder — copied to blast_hits_only_EPEC/"
  else
      echo "❌ No hits for $folder"
  fi
done
```

## Blast STEC
```
mkdir -p blast_out_STEC
mkdir -p blast_hits_only_STEC

for f in $(find /home/jing/E.coli/blast_results/STEC_assemblies -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out_STEC/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query all_markers.fasta \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=90 && cov>=90) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only_STEC/
      echo "✅ Hits found in $folder — copied to blast_hits_only_STEC/"
  else
      echo "❌ No hits for $folder"
  fi
done
```
## Make presence and absence martrix
```
#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd

def read_markers_from_fasta(fasta_path):
    """
    Read marker IDs from FASTA headers.
    Assumes headers like: >markerID something something
    """
    markers = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                marker_id = line[1:].strip().split()[0]
                markers.append(marker_id)
    return markers

def main(blast_dir, fasta_path, out_tsv):
    # 1. Get list of all markers from FASTA (even if some never hit)
    markers = read_markers_from_fasta(fasta_path)
    markers = sorted(set(markers))

    # 2. Find all BLAST result files
    blast_files = sorted(glob.glob(os.path.join(blast_dir, "*_hits.tsv")))
    if not blast_files:
        raise SystemExit(f"No *_hits.tsv files found in {blast_dir}")

    # 3. Build presence/absence dict: { sample_id: {marker: 0/1} }
    presence = {}

    for bf in blast_files:
        # Extract sample ID from file name: e.g. folder_base_hits.tsv → folder_base
        fname = os.path.basename(bf)
        sample_id = fname.replace("_hits.tsv", "")

        # Initialize all markers as absent (0)
        presence[sample_id] = {m: 0 for m in markers}

        with open(bf) as f:
            for line in f:
                parts = line.strip().split("\t")
                if not parts or len(parts) < 1:
                    continue
                qseqid = parts[0]  # marker ID from BLAST qseqid
                if qseqid in presence[sample_id]:
                    presence[sample_id][qseqid] = 1

    # 4. Convert to DataFrame: rows = samples, cols = markers
    df = pd.DataFrame.from_dict(presence, orient="index")
    df.index.name = "Sample"
    df = df[sorted(df.columns)]  # consistently ordered markers

    # 5. Write to TSV
    df.to_csv(out_tsv, sep="\t")

    print(f"✅ Matrix saved to: {out_tsv}")
    print(f"   Samples: {df.shape[0]}")
    print(f"   Markers: {df.shape[1]}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make marker presence/absence matrix from BLAST hits."
    )
    parser.add_argument(
        "-b", "--blast_dir",
        required=True,
        help="Directory with *_hits.tsv files (e.g. blast_out or blast_hits_only)"
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="FASTA file with all markers (e.g. all_markers.fasta)"
    )
    parser.add_argument(
        "-o", "--out",
        default="marker_presence_absence.tsv",
        help="Output TSV file (default: marker_presence_absence.tsv)"
    )

    args = parser.parse_args()
    main(args.blast_dir, args.fasta, args.out)
```
## Run the script
```
python make_marker_matrix.py \
  -b blast_out_EHEC \
  -f all_markers.fasta \
  -o marker_presence_absence_EHEC.tsv
```
```
python make_marker_matrix.py \
  -b blast_out_EPEC \
  -f all_markers.fasta \
  -o marker_presence_absence_EPEC.tsv

python make_marker_matrix.py \
  -b blast_out_STEC \
  -f all_markers.fasta \
  -o marker_presence_absence_STEC.tsv
```

## Count all espK/espV/espN combinations with in EHEC
count_EHEC_esp.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your EHEC-only table
df = pd.read_csv("marker_presence_absence_EHCE.tsv", sep="\t")

# 2. Boolean presence for each gene
K = df["espK"] == 1
V = df["espV"] == 1
N = df["espN"] == 1

total = len(df)

rows = []

def add(label, mask):
    count = int(mask.sum())
    pct = round(count / total * 100, 2) if total > 0 else 0.0
    rows.append({"Pattern": label, "Count": count, "Percent(%)": pct})

# 3. Individual genes
add("espK(+)", K)
add("espV(+)", V)
add("espN(+)", N)

# 4. Simple combinations (AND / OR)
add("espK(+) AND espV(+)", K & V)
add("espK(+) AND espN(+)", K & N)
add("espV(+) AND espN(+)", V & N)
add("espK(+) AND espV(+) AND espN(+)", K & V & N)
add("espK(+) OR espV(+) OR espN(+)", K | V | N)

# 5. 8 mutually exclusive patterns (K/V/N all combos)
onlyK =  K & ~V & ~N
onlyV = ~K &  V & ~N
onlyN = ~K & ~V &  N
KV    =  K &  V & ~N
KN    =  K & ~V &  N
VN    = ~K &  V &  N
KVN   =  K &  V &  N
none  = ~K & ~V & ~N

add("espK only (K+ V- N-)", onlyK)
add("espV only (K- V+ N-)", onlyV)
add("espN only (K- V- N+)", onlyN)
add("espK & espV (K+ V+ N-)", KV)
add("espK & espN (K+ V- N+)", KN)
add("espV & espN (K- V+ N+)", VN)
add("espK & espV & espN (K+ V+ N+)", KVN)
add("none esp (K- V- N-)", none)

# 6. Make DataFrame and save / print
summary = pd.DataFrame(rows)

print(f"Total EHEC isolates: {total}\n")
print(summary)

summary.to_csv("EHEC_esp_counts_percentages.csv", index=False)
print("\nSaved to: EHEC_esp_counts_percentages.csv")
```
<img width="445" height="295" alt="image" src="https://github.com/user-attachments/assets/f5dcd2bc-c708-4ac9-8cd3-309cff08b7ef" />


## Count all espK/espV/espN combinations with in EPEC
count_PHEC_esp.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your EPEC-only table
df = pd.read_csv("marker_presence_absence_EPEC.tsv", sep="\t")

# 2. Boolean presence for each gene
K = df["espK"] == 1
V = df["espV"] == 1
N = df["espN"] == 1

total = len(df)

rows = []

def add(label, mask):
    count = int(mask.sum())
    pct = round(count / total * 100, 2) if total > 0 else 0.0
    rows.append({"Pattern": label, "Count": count, "Percent(%)": pct})

# 3. Individual genes
add("espK(+)", K)
add("espV(+)", V)
add("espN(+)", N)

# 4. Simple combinations (AND / OR)
add("espK(+) AND espV(+)", K & V)
add("espK(+) AND espN(+)", K & N)
add("espV(+) AND espN(+)", V & N)
add("espK(+) AND espV(+) AND espN(+)", K & V & N)
add("espK(+) OR espV(+) OR espN(+)", K | V | N)

# 5. 8 mutually exclusive patterns (K/V/N all combos)
onlyK =  K & ~V & ~N
onlyV = ~K &  V & ~N
onlyN = ~K & ~V &  N
KV    =  K &  V & ~N
KN    =  K & ~V &  N
VN    = ~K &  V &  N
KVN   =  K &  V &  N
none  = ~K & ~V & ~N

add("espK only (K+ V- N-)", onlyK)
add("espV only (K- V+ N-)", onlyV)
add("espN only (K- V- N+)", onlyN)
add("espK & espV (K+ V+ N-)", KV)
add("espK & espN (K+ V- N+)", KN)
add("espV & espN (K- V+ N+)", VN)
add("espK & espV & espN (K+ V+ N+)", KVN)
add("none esp (K- V- N-)", none)

# 6. Make DataFrame and save / print
summary = pd.DataFrame(rows)

print(f"Total EPEC isolates: {total}\n")
print(summary)

summary.to_csv("EPEC_esp_counts_percentages.csv", index=False)
print("\nSaved to: EPEC_esp_counts_percentages.csv")
```
<img width="457" height="296" alt="image" src="https://github.com/user-attachments/assets/fb7e05df-f688-480f-bee4-8fdc0550f1ac" />

## Count all espK/espV/espN combinations with in STEC
count_STEC_esp.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your STEC-only table
df = pd.read_csv("marker_presence_absence_STEC.tsv", sep="\t")

# 2. Boolean presence for each gene
K = df["espK"] == 1
V = df["espV"] == 1
N = df["espN"] == 1

total = len(df)

rows = []

def add(label, mask):
    count = int(mask.sum())
    pct = round(count / total * 100, 2) if total > 0 else 0.0
    rows.append({"Pattern": label, "Count": count, "Percent(%)": pct})

# 3. Individual genes
add("espK(+)", K)
add("espV(+)", V)
add("espN(+)", N)

# 4. Simple combinations (AND / OR)
add("espK(+) AND espV(+)", K & V)
add("espK(+) AND espN(+)", K & N)
add("espV(+) AND espN(+)", V & N)
add("espK(+) AND espV(+) AND espN(+)", K & V & N)
add("espK(+) OR espV(+) OR espN(+)", K | V | N)

# 5. 8 mutually exclusive patterns (K/V/N all combos)
onlyK =  K & ~V & ~N
onlyV = ~K &  V & ~N
onlyN = ~K & ~V &  N
KV    =  K &  V & ~N
KN    =  K & ~V &  N
VN    = ~K &  V &  N
KVN   =  K &  V &  N
none  = ~K & ~V & ~N

add("espK only (K+ V- N-)", onlyK)
add("espV only (K- V+ N-)", onlyV)
add("espN only (K- V- N+)", onlyN)
add("espK & espV (K+ V+ N-)", KV)
add("espK & espN (K+ V- N+)", KN)
add("espV & espN (K- V+ N+)", VN)
add("espK & espV & espN (K+ V+ N+)", KVN)
add("none esp (K- V- N-)", none)

# 6. Make DataFrame and save / print
summary = pd.DataFrame(rows)

print(f"Total STEC isolates: {total}\n")
print(summary)

summary.to_csv("STEC_esp_counts_percentages.csv", index=False)
print("\nSaved to: STEC_esp_counts_percentages.csv")
```

## Step 2.1: Filter EHEC strains with espK/espV/espN positive

filter_esp_positive_EHEC.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your EHEC presence/absence table
df = pd.read_csv("marker_presence_absence_EHCE.tsv", sep="\t")

# 🔴 Change this if your ID column has a different name
id_col = "Sample"

# 2. Boolean mask for espK(+) OR espV(+) OR espN(+)
mask_any_esp = (df["espK"] == 1) | (df["espV"] == 1) | (df["espN"] == 1)

subset = df[mask_any_esp].copy()

print(f"Total EHEC isolates: {len(df)}")
print(f"EHEC with espK/V/N positive: {len(subset)}")

# 3. Save list of IDs (one per line) for looping
subset_ids = subset[id_col].drop_duplicates()
subset_ids.to_csv("serotype/EHEC_esp_any_positive_ids.txt",
                  index=False, header=False)

# Optional: save full metadata for those
subset.to_csv("serotype/EHEC_esp_any_positive_metadata.tsv",
              sep="\t", index=False)

print("Saved IDs to: EHEC_esp_any_positive_ids.txt")
print("Saved metadata to: EHEC_esp_any_positive_metadata.tsv")
```

```
cd /home/jing/E.coli/blast_results/linkage/serotype

awk -F'EHEC_assemblies_' '{print $2}' EHEC_esp_any_positive_ids.txt \
  > EHEC_esp_any_positive_ids.core.txt
```
## Serotyping with ectyper
### EHEC
run_ectyper_esp_positive_EHEC.sh
```
#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="/home/jing/E.coli/blast_results/linkage/serotype/EHEC_esp_any_positive_ids.core.txt"
ASSEMBLY_ROOT="/home/jing/E.coli/blast_results/EHEC_assemblies"
OUTDIR="/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out_EHEC"
mkdir -p "$OUTDIR"

while read -r core; do
    [ -z "$core" ] && continue

    echo "=== Processing core: $core"

    fna=$(find "$ASSEMBLY_ROOT" -type f -name "${core}.fna" | head -n 1)

    if [ -z "$fna" ]; then
        echo "    ❌ No .fna found for $core"
        continue
    fi

    sample_out="${OUTDIR}/${core}"
    mkdir -p "$sample_out"

    ectyper -i "$fna" -o "$sample_out"

    echo "    ✅ ECTyper done -> $sample_out"
done < "$IDS_FILE"
```
```
chmod +x run_ectyper_esp_positive_EHEC.sh
```
```
cd /home/jing/E.coli/blast_results/linkage/serotype
./run_ectyper_esp_positive_EHEC.sh
```


### EPEC
run_ectyper_esp_positive_EPEC.sh
```
#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="/home/jing/E.coli/blast_results/linkage/serotype/EPEC_esp_any_positive_ids.core.txt"
ASSEMBLY_ROOT="/home/jing/E.coli/blast_results/EPEC_assemblies"
OUTDIR="/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out_EPEC"
mkdir -p "$OUTDIR"

while read -r core; do
    [ -z "$core" ] && continue

    echo "=== Processing core: $core"

    fna=$(find "$ASSEMBLY_ROOT" -type f -name "${core}.fna" | head -n 1)

    if [ -z "$fna" ]; then
        echo "    ❌ No .fna found for $core"
        continue
    fi

    sample_out="${OUTDIR}/${core}"
    mkdir -p "$sample_out"

    ectyper -i "$fna" -o "$sample_out"

    echo "    ✅ ECTyper done -> $sample_out"
done < "$IDS_FILE"
```
```
chmod +x run_ectyper_esp_positive_EPEC.sh
```
```
cd /home/jing/E.coli/blast_results/linkage/serotype
./run_ectyper_esp_positive_EPEC.sh
```
## Summarize EPEC_ECTyper
summarize_ectyper_EPEC.py
```
#!/usr/bin/env python3
import os
import pandas as pd

BASE = "/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out_EPEC"

rows = []

for sample in sorted(os.listdir(BASE)):
    sample_dir = os.path.join(BASE, sample)
    if not os.path.isdir(sample_dir):
        continue

    out_tsv = os.path.join(sample_dir, "output.tsv")
    if not os.path.exists(out_tsv):
        print(f"⚠️  No output.tsv in {sample_dir}, skipping")
        continue

    df = pd.read_csv(out_tsv, sep="\t")

    if df.empty:
        print(f"⚠️  Empty output.tsv in {sample_dir}, skipping")
        continue

    row = df.iloc[0].to_dict()

    # Make sure sample name is included (folder name is safest)
    row["Sample_folder"] = sample
    # If Name column exists, keep it too
    # (E.g., Name = GCA_000167875.2_ASM16787v2_genomic)
    # You’ll see both in the final table
    rows.append(row)

# Build DataFrame
summary = pd.DataFrame(rows)

# If you only want a subset of columns, specify them here:
preferred_cols = [
    "Sample_folder",
    "Name",
    "Species",
    "O-type",
    "H-type",
    "Serotype",
    "QC",
    "Evidence",
    "Pathotype",
    "StxSubtypes"
]
# Keep only those that actually exist in the DataFrame
cols = [c for c in preferred_cols if c in summary.columns]
# Plus any remaining columns at the end (optional)
other_cols = [c for c in summary.columns if c not in cols]
summary = summary[cols + other_cols]

out_file = "/home/jing/E.coli/blast_results/linkage/serotype/EPEC_ECTyper_summary.tsv"
summary.to_csv(out_file, sep="\t", index=False)

print(f"\n✅ Saved summary to: {out_file}")
print("\nPreview:")
print(summary.head())
```
```
python3 /home/jing/E.coli/blast_results/linkage/serotype/summarize_ectyper_EPEC.py
```
### Merge esp markers with serotype

merge_EPEC_esp_only_serotype.py

```
#!/usr/bin/env python3
import pandas as pd
import re

# === Paths ===
matrix_path = "/home/jing/E.coli/blast_results/linkage/marker_presence_absence_EPEC.tsv"
ectyper_path = "/home/jing/E.coli/blast_results/linkage/serotype/EPEC_ECTyper_summary.tsv"
out_path = "/home/jing/E.coli/blast_results/linkage/EPEC_esp_with_serotype.tsv"

# Load datasets
matrix = pd.read_csv(matrix_path, sep="\t")
ect = pd.read_csv(ectyper_path, sep="\t")

# === 1. Detect ID columns ===
matrix_id_col = "Sample" if "Sample" in matrix.columns else matrix.columns[0]
ect_id_col = "Sample_folder" if "Sample_folder" in ect.columns else "Name"

# === 2. Clean IDs (remove prefixes, extract GCA_*/GCF_*) ===
def extract_accession(x):
    m = re.search(r'(GC[AF]_[0-9]+\.[0-9]+_[A-Za-z0-9]+)', str(x))
    return m.group(1) if m else str(x)

matrix["CleanID"] = matrix[matrix_id_col].apply(extract_accession)
ect["CleanID"] = ect[ect_id_col].apply(extract_accession)

# === 3. Extract ONLY espK, espV, espN ===
esp_cols = [col for col in matrix.columns if col.lower() in ["espk", "espv", "espn"]]

if len(esp_cols) == 0:
    raise ValueError("❌ ERROR: espK, espV, espN not found in matrix file.")

matrix_sub = matrix[["CleanID"] + esp_cols]

# === 4. Extract ONLY Serotype column ===
if "Serotype" not in ect.columns:
    raise ValueError("❌ ERROR: 'Serotype' column not found in ECTyper summary!")

ect_sub = ect[["CleanID", "Serotype"]]

# === 5. Merge ===
merged = matrix_sub.merge(ect_sub, on="CleanID", how="left")

# Rename CleanID to Sample
merged = merged.rename(columns={"CleanID": "Sample"})

# === 6. Save output ===
merged.to_csv(out_path, sep="\t", index=False)

print(f"\n✅ Saved clean table: {out_path}")
print("Preview:")
print(merged.head())
```
### assemblies that have at least one esp gene (espK OR espV OR espN = 1).
filter_EPEC_esp_positive_only.py
```
#!/usr/bin/env python3
import pandas as pd

in_path = "/home/jing/E.coli/blast_results/linkage/EPEC_esp_with_serotype.tsv"
out_path = "/home/jing/E.coli/blast_results/linkage/EPEC_espPOS_with_serotype.tsv"

df = pd.read_csv(in_path, sep="\t")

# Make sure esp columns exist
for col in ["espK", "espV", "espN"]:
    if col not in df.columns:
        raise ValueError(f"Column {col} not found in {in_path}")

# Filter: any esp gene present
mask = (df["espK"] == 1) | (df["espV"] == 1) | (df["espN"] == 1)
df_pos = df[mask].copy()

df_pos.to_csv(out_path, sep="\t", index=False)

print(f"Total isolates: {len(df)}")
print(f"esp+ isolates (espK/espV/espN): {len(df_pos)}")
print(f"Saved to: {out_path}")
print(df_pos.head())
```

### Calculate the distribution
```
#!/usr/bin/env python3
import pandas as pd

# === Load EHEC and EPEC esp-positive tables ===
ehec = pd.read_csv("EHEC_espPOS_with_serotype.tsv", sep="\t")
epec = pd.read_csv("EPEC_espPOS_with_serotype.tsv", sep="\t")

def make_summary(df, label):
    """Summarize serotype distribution + espK/V/N for a cohort."""

    # Group by serotype
    summary = df.groupby("Serotype").agg(
        N_total = ("Sample", "count"),
        espK_pos = ("espK", "sum"),
        espV_pos = ("espV", "sum"),
        espN_pos = ("espN", "sum"),
    ).reset_index()

    # Percentages
    summary["espK_percent"] = (summary["espK_pos"] / summary["N_total"] * 100).round(2)
    summary["espV_percent"] = (summary["espV_pos"] / summary["N_total"] * 100).round(2)
    summary["espN_percent"] = (summary["espN_pos"] / summary["N_total"] * 100).round(2)

    summary["Cohort"] = label
    return summary


# === Build summaries ===
ehec_sum = make_summary(ehec, "EHEC")
epec_sum = make_summary(epec, "EPEC")

# === Save results ===
ehec_sum.to_csv("EHEC_serotype_esp_distribution.tsv", sep="\t", index=False)
epec_sum.to_csv("EPEC_serotype_esp_distribution.tsv", sep="\t", index=False)

# Combined file
combined = pd.concat([ehec_sum, epec_sum], ignore_index=True)
combined.to_csv("EHEC_EPEC_serotype_esp_distribution.tsv", sep="\t", index=False)

print("\n=== EHEC Summary ===")
print(ehec_sum.head())

print("\n=== EPEC Summary ===")
print(epec_sum.head())

print("\nSaved:")
print(" - EHEC_serotype_esp_distribution.tsv")
print(" - EPEC_serotype_esp_distribution.tsv")
print(" - EHEC_EPEC_serotype_esp_distribution.tsv")
```
#### Plot distribution
plot_serotypes_horizontal_single_color.py
```
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

# ---------- SETTINGS ----------
EHEC_FILE = "EHEC_serotype_esp_distribution.tsv"
EPEC_FILE = "EPEC_serotype_esp_distribution.tsv"

# Sort options: "N_total", "espK_percent", "espV_percent", "espN_percent"
SORT_BY = "N_total"

TOP_N = 30          # show top N serotypes
BAR_COLOR = "tab:blue"   # <-- ONE COLOR FOR ALL BARS
# ------------------------------


def make_barh(ax, df, title, sort_by=SORT_BY):
    """Plot a horizontal barplot with a single color and bar labels."""

    df = df.copy()

    # Sort by selected metric
    if sort_by not in df.columns:
        print(f"⚠ SORT_BY '{sort_by}' not found; using 'N_total'")
        sort_by = "N_total"

    df = df.sort_values(sort_by, ascending=True).tail(TOP_N)

    y = df["Serotype"]
    x = df["N_total"]

    ax.barh(y, x, color=BAR_COLOR)

    # Title / labels
    ax.set_title(title)
    ax.set_xlabel("Number of esp+ isolates")
    ax.set_ylabel("Serotype")

    # Add labels showing espK/V/N %
    for i, row in enumerate(df.itertuples()):
        label = (
            f"{int(row.N_total)} | "
            f"K:{row.espK_percent:.1f}% "
            f"V:{row.espV_percent:.1f}% "
            f"N:{row.espN_percent:.1f}%"
        )
        ax.text(
            row.N_total, i,
            label,
            va="center", ha="left",
            fontsize=8
        )

    ax.margins(x=0.15)


def main():
    # Load tables
    ehec = pd.read_csv(EHEC_FILE, sep="\t")
    epec = pd.read_csv(EPEC_FILE, sep="\t")

    # Figure height auto-scales
    max_rows = max(len(ehec), len(epec), 10)
    fig, axes = plt.subplots(1, 2, figsize=(18, max(6, max_rows * 0.25)))

    # EHEC plot
    make_barh(
        axes[0],
        ehec,
        title=f"EHEC esp+ Serotype Distribution (sorted by {SORT_BY})"
    )

    # EPEC plot
    make_barh(
        axes[1],
        epec,
        title=f"EPEC esp+ Serotype Distribution (sorted by {SORT_BY})"
    )

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
```
<img width="2000" height="1035" alt="image" src="https://github.com/user-attachments/assets/b1293c7d-a20a-42a9-a146-ee654d70d63e" />


## Step 2.1: Filter EPEC strains with espK/espV/espN ALL negative

1. filter_esp_negative_EPEC.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your EPEC presence/absence table
df = pd.read_csv("marker_presence_absence_EPEC.tsv", sep="\t")

# 🔴 Change this if your ID column has a different name
id_col = "Sample"

# 2. Boolean mask for espK(-) AND espV(-) AND espN(-)
mask_all_esp_neg = (df["espK"] == 0) & (df["espV"] == 0) & (df["espN"] == 0)

subset = df[mask_all_esp_neg].copy()

print(f"Total EPEC isolates: {len(df)}")
print(f"EPEC with espK/V/N all negative: {len(subset)}")

# 3. Save list of IDs (one per line) for looping
subset_ids = subset[id_col].drop_duplicates()
subset_ids.to_csv("serotype/EPEC_esp_all_negative_ids.txt",
                  index=False, header=False)

# Optional: save full metadata for those
subset.to_csv("serotype/EPEC_esp_all_negative_metadata.tsv",
              sep="\t", index=False)

print("Saved IDs to: serotype/EPEC_esp_all_negative_ids.txt")
print("Saved metadata to: serotype/EPEC_esp_all_negative_metadata.tsv")
```
2. Strip the EPEC_assemblies_ prefix to match assembly filenames
```
cd /home/jing/E.coli/blast_results/linkage/serotype

awk -F'EPEC_assemblies_' '{print $2}' EPEC_esp_all_negative_ids.txt \
  > EPEC_esp_all_negative_ids.core.txt
```
3. Serotyping esp– EPEC with ECTyper
run_ectyper_esp_negative_EPEC.sh
```
#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="/home/jing/E.coli/blast_results/linkage/serotype/EPEC_esp_all_negative_ids.core.txt"
ASSEMBLY_ROOT="/home/jing/E.coli/blast_results/EPEC_assemblies"
OUTDIR="/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out_EPEC_espNEG"

mkdir -p "$OUTDIR"

while read -r core; do
    [ -z "$core" ] && continue

    echo "=== Processing core: $core"

    fna=$(find "$ASSEMBLY_ROOT" -type f -name "${core}.fna" | head -n 1)

    if [ -z "$fna" ]; then
        echo "    ❌ No .fna found for $core"
        continue
    fi

    sample_out="${OUTDIR}/${core}"
    mkdir -p "$sample_out"

    ectyper -i "$fna" -o "$sample_out"

    echo "    ✅ ECTyper done -> $sample_out"
done < "$IDS_FILE"
```

4. Summarize ECTyper outputs for EPEC esp– isolates
summarize_ectyper_EPEC_espNEG.py
```
#!/usr/bin/env python3
import os
import pandas as pd

BASE = "/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out_EPEC_espNEG"

rows = []

for sample in sorted(os.listdir(BASE)):
    sample_dir = os.path.join(BASE, sample)
    if not os.path.isdir(sample_dir):
        continue

    out_tsv = os.path.join(sample_dir, "output.tsv")
    if not os.path.exists(out_tsv):
        print(f"⚠️  No output.tsv in {sample_dir}, skipping")
        continue

    df = pd.read_csv(out_tsv, sep="\t")
    if df.empty:
        print(f"⚠️  Empty output.tsv in {sample_dir}, skipping")
        continue

    row = df.iloc[0].to_dict()
    row["Sample_folder"] = sample
    rows.append(row)

summary = pd.DataFrame(rows)

# Put Sample_folder first
cols = ["Sample_folder"] + [c for c in summary.columns if c != "Sample_folder"]
summary = summary[cols]

out_file = "/home/jing/E.coli/blast_results/linkage/serotype/EPEC_espNEG_ECTyper_summary.tsv"
summary.to_csv(out_file, sep="\t", index=False)

print(f"\n✅ Saved espNEG ECTyper summary to: {out_file}")
print("Preview:")
print(summary.head())
```

5. Merge esp– metadata with Serotype
merge_EPEC_espNEG_serotype_only.py
```
#!/usr/bin/env python3
import pandas as pd
import re

meta_path = "serotype/EPEC_esp_all_negative_metadata.tsv"
ect_path = "serotype/EPEC_espNEG_ECTyper_summary.tsv"
out_path = "EPEC_espNEG_with_serotype.tsv"

# Load
meta = pd.read_csv(meta_path, sep="\t")
ect = pd.read_csv(ect_path, sep="\t")

# ID columns
meta_id_col = "Sample" if "Sample" in meta.columns else meta.columns[0]
ect_id_col = "Sample_folder" if "Sample_folder" in ect.columns else "Name"

# Extract GCA/GCF core id
def extract_accession(x):
    m = re.search(r'(GC[AF]_[0-9]+\.[0-9]+_[A-Za-z0-9]+)', str(x))
    return m.group(1) if m else str(x)

meta["CleanID"] = meta[meta_id_col].apply(extract_accession)
ect["CleanID"] = ect[ect_id_col].apply(extract_accession)

# Just espK/espV/espN + CleanID
esp_cols = [c for c in meta.columns if c.lower() in ["espk", "espv", "espn"]]
if len(esp_cols) == 0:
    raise ValueError("No espK/espV/espN columns found in metadata.")

meta_sub = meta[["CleanID"] + esp_cols]

# Serotype only
if "Serotype" not in ect.columns:
    raise ValueError("No 'Serotype' column found in ECTyper summary.")

ect_sub = ect[["CleanID", "Serotype"]]

# Merge
merged = meta_sub.merge(ect_sub, on="CleanID", how="left")
merged = merged.rename(columns={"CleanID": "Sample"})

merged.to_csv(out_path, sep="\t", index=False)

print(f"\n✅ Saved espNEG markers+Serotype to: {out_path}")
print(merged.head())
```
6. Compare serotype distributions: esp+ vs esp– EPEC
EPEC_serotype_espPOS_vs_espNEG.py
```
#!/usr/bin/env python3
import pandas as pd

pos_path = "EPEC_espPOS_with_serotype.tsv"
neg_path = "EPEC_espNEG_with_serotype.tsv"
out_path = "EPEC_serotype_espPOS_vs_espNEG.tsv"

pos = pd.read_csv(pos_path, sep="\t").dropna(subset=["Serotype"])
neg = pd.read_csv(neg_path, sep="\t").dropna(subset=["Serotype"])

# Flag each set
pos["esp_status"] = "espPOS"
neg["esp_status"] = "espNEG"

# Count per serotype per status
pos_counts = (
    pos.groupby("Serotype")
       .size()
       .reset_index(name="espPOS_count")
)
neg_counts = (
    neg.groupby("Serotype")
       .size()
       .reset_index(name="espNEG_count")
)

# Combine
summary = pd.merge(pos_counts, neg_counts, on="Serotype", how="outer").fillna(0)

# Total and percentages within each serotype
summary["Total"] = summary["espPOS_count"] + summary["espNEG_count"]
summary["espPOS_percent"] = (summary["espPOS_count"] / summary["Total"] * 100).round(2)
summary["espNEG_percent"] = (summary["espNEG_count"] / summary["Total"] * 100).round(2)

summary = summary.sort_values("Total", ascending=False)

summary.to_csv(out_path, sep="\t", index=False)

print(f"\n✅ Saved comparison to: {out_path}")
print(summary.head())
```
7. Plot the distribution of both
plot_espPOS_vs_espNEG.py
```
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

TOP_N = 20

# Colors
COLOR_NEG = "#cfcfcf"        # light grey
COLOR_POS = "#4a90e2"        # steel blue
TEXT_DARK = "#333333"
TEXT_LIGHT = "white"


def inside_pretty_plot(input_tsv, cohort_label, out_prefix):

    df = pd.read_csv(input_tsv, sep="\t")
    df = df.sort_values("Total", ascending=True).tail(TOP_N).reset_index(drop=True)

    serotypes = df["Serotype"]
    neg_pct = df["espNEG_percent"]
    pos_pct = df["espPOS_percent"]
    neg_cnt = df["espNEG_count"]
    pos_cnt = df["espPOS_count"]

    # Base plot
    plt.figure(figsize=(13, max(6, len(df) * 0.45)))
    plt.style.use("seaborn-v0_8-muted")

    # Draw stacked bars
    plt.barh(serotypes, neg_pct, color=COLOR_NEG, edgecolor="white", label="espNEG%")
    plt.barh(serotypes, pos_pct, left=neg_pct, color=COLOR_POS, edgecolor="white", label="espPOS%")

    plt.title(
        f"{cohort_label}: espPOS% vs espNEG% per Serotype (Top {TOP_N})",
        fontsize=18, weight="bold", pad=20, color=TEXT_DARK
    )

    plt.xlabel("Percentage of isolates (%)", fontsize=14, weight="bold")
    plt.ylabel("Serotype", fontsize=14, weight="bold")
    plt.yticks(fontsize=12, color=TEXT_DARK)
    plt.grid(axis="x", color="#eeeeee", linewidth=1.2)

    xmax = max((neg_pct + pos_pct).max(), 100) * 1.02
    plt.xlim(0, xmax)

    # === LABELS INSIDE THE BAR ===
    for i, row in df.iterrows():

        # --- NEGATIVE section label ---
        if row["espNEG_percent"] > 5:
            label_neg = f"NEG {row['espNEG_percent']:.1f}% ({row['espNEG_count']})"
            plt.text(
                row["espNEG_percent"] * 0.5,
                i,
                label_neg,
                va="center",
                ha="center",
                fontsize=10,
                color=TEXT_DARK
            )

        # --- POSITIVE section label ---
        if row["espPOS_percent"] > 5:
            label_pos = f"POS {row['espPOS_percent']:.1f}% ({row['espPOS_count']})"
            plt.text(
                row["espNEG_percent"] + row["espPOS_percent"] * 0.5,
                i,
                label_pos,
                va="center",
                ha="center",
                fontsize=10,
                color=TEXT_LIGHT if row["espPOS_percent"] > 20 else TEXT_DARK
            )

        # If a bar segment is <5%, label OUTSIDE (to avoid overlap)
        if row["espPOS_percent"] <= 5:
            label_pos = f"POS {row['espPOS_percent']:.1f}% ({row['espPOS_count']})"
            plt.text(
                row["espNEG_percent"] + row["espPOS_percent"] + 2,
                i,
                label_pos,
                va="center",
                ha="left",
                fontsize=9,
                color=TEXT_DARK
            )

        if row["espNEG_percent"] <= 5:
            label_neg = f"NEG {row['espNEG_percent']:.1f}% ({row['espNEG_count']})"
            plt.text(
                0.5,
                i,
                label_neg,
                va="center",
                ha="left",
                fontsize=9,
                color=TEXT_DARK
            )

    # Legend
    plt.legend(
        frameon=True, fontsize=12, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=2
    )

    plt.tight_layout()

    # Save
    plt.savefig(f"{out_prefix}.png", dpi=400, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight")

    plt.show()
    print(f"Saved {out_prefix}.png and {out_prefix}.pdf")


def main():
    inside_pretty_plot(
        "EPEC_serotype_espPOS_vs_espNEG.tsv",
        cohort_label="EPEC",
        out_prefix="EPEC_espPOS_vs_espNEG_top20_inside"
    )

    inside_pretty_plot(
        "EHEC_serotype_espPOS_vs_espNEG.tsv",
        cohort_label="EHEC",
        out_prefix="EHEC_espPOS_vs_espNEG_top20_inside"
    )


if __name__ == "__main__":
    main()
```
<img width="1490" height="1125" alt="image" src="https://github.com/user-attachments/assets/c9b4d65a-7e62-496d-8d9d-7ae5ee580dd3" />


8. Plot Distribution of esp-Negative Serotypes Only
plot_espNEG_only.py
```
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

TOP_N = 20  # how many serotypes to show

COLOR_NEG = "#cfcfcf"          # grey
TEXT_DARK = "#333333"


def plot_neg_only(input_file, cohort_label, out_prefix):

    # Load espNEG table
    df = pd.read_csv(input_file, sep="\t")

    # Drop missing serotypes
    df = df.dropna(subset=["Serotype"])

    # Count serotypes
    sero_counts = (
        df.groupby("Serotype").size().reset_index(name="count")
    )

    # Add percentages
    total = sero_counts["count"].sum()
    sero_counts["percent"] = sero_counts["count"] / total * 100

    # Sort & keep top
    sero_counts = sero_counts.sort_values("count", ascending=True).tail(TOP_N)
    sero_counts = sero_counts.reset_index(drop=True)

    # Plot size
    plt.figure(figsize=(12, max(6, len(sero_counts) * 0.4)))
    plt.style.use("seaborn-v0_8-muted")

    # Horizontal barplot
    plt.barh(
        sero_counts["Serotype"],
        sero_counts["percent"],
        color=COLOR_NEG,
        edgecolor="white",
        label="espNEG only"
    )

    # Title & axis
    plt.title(
        f"{cohort_label}: espNEG-only Serotype Distribution (Top {TOP_N})",
        fontsize=18, weight="bold", pad=20, color=TEXT_DARK
    )
    plt.xlabel("Percentage of isolates (%)", fontsize=14, weight="bold")
    plt.ylabel("Serotype", fontsize=14, weight="bold")
    plt.yticks(fontsize=12, color=TEXT_DARK)
    plt.grid(axis="x", color="#eeeeee", linewidth=1.2)

    xmax = sero_counts["percent"].max() * 1.18
    plt.xlim(0, xmax)

    # Add labels inside the bar
    for i, row in sero_counts.iterrows():
        label = f"{row['percent']:.1f}%  ({row['count']})"
        plt.text(
            row["percent"] * 0.5,
            i,
            label,
            ha="center",
            va="center",
            fontsize=10,
            color=TEXT_DARK
        )

    plt.tight_layout()

    # Save
    plt.savefig(f"{out_prefix}.png", dpi=400, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight")

    plt.show()
    print(f"Saved {out_prefix}.png and {out_prefix}.pdf")


def main():
    # EPEC espNEG only
    plot_neg_only(
        "EPEC_espNEG_with_serotype.tsv",
        "EPEC",
        "EPEC_espNEG_only_distribution"
    )

    # EHEC espNEG only (optional — enable if needed)
    # plot_neg_only(
    #     "EHEC_espNEG_with_serotype.tsv",
    #     "EHEC",
    #     "EHEC_espNEG_only_distribution"
    # )


if __name__ == "__main__":
    main()
```
<img width="1181" height="799" alt="image" src="https://github.com/user-attachments/assets/410bbc6d-7d54-4fbf-8fd5-21e4bed7c406" />

