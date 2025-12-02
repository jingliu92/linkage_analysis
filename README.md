# linkage_analysis
Rationale: esp+ genotype is ~99% associated with EHEC but it also associated with EPEC pathovar.  This confounds the use of esp genotypes to screen for EHEC pathovars. Here this analysis will determine whether a particular serotype is associate specifically with esp positive and esp negative EPEC pathovars. If such linkage is found, serotype markers can be especially used to identify esp positve EPEC pathovars. 

1. EHEC with espK+, espV+, espN+

2. EPEC with espK+, espV+, espN+

3. EPEC with espK-, espV-, espN-

To do this, three steps are required:
step 1: Classify the distribution of esp+ genes within EHEC or EPEC
step 2: extract assemblies with/without espK+, espV+, espN+ for EHEC and EPEC
step 3: serotyping
## Create marker gene file
nano all.gene.fasta

## make db
```
makeblastdb -in all_markers.fasta -dbtype nucl -out all_markers
```
## Blast EHEC
```
mkdir -p blast_out
mkdir -p blast_hits_only

for f in $(find /home/jing/E.coli/blast_results/EHEC_assemblies -name "*.fna"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query all_markers.fasta \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=80 && cov>=80) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only/
      echo "✅ Hits found in $folder — copied to blast_hits_only/"
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
         | awk '{cov=($4/$5)*100; if($3>=80 && cov>=80) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only_EPEC/
      echo "✅ Hits found in $folder — copied to blast_hits_only_EPEC/"
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
  -o marker_presence_absence_EHCE.tsv
```
```
python make_marker_matrix.py \
  -b blast_out_EPEC \
  -f all_markers.fasta \
  -o marker_presence_absence_EPEC.tsv
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
<img width="424" height="267" alt="image" src="https://github.com/user-attachments/assets/bbf02961-c7ea-4737-b0d1-821386de40f6" />

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

print(f"Total EHEC isolates: {total}\n")
print(summary)

summary.to_csv("EPEC_esp_counts_percentages.csv", index=False)
print("\nSaved to: EPEC_esp_counts_percentages.csv")
```
<img width="464" height="294" alt="image" src="https://github.com/user-attachments/assets/99772169-4d96-48fe-a4c4-5e1c64dabe3d" />


Step 2: Filter EHEC strains with espK/espV/espN positive
```
cd /home/jing/E.coli/blast_results/linkage/serotype

awk -F'EHEC_assemblies_' '{print $2}' EHEC_esp_any_positive_ids.txt \
  > EHEC_esp_any_positive_ids.core.txt
```
filter_esp_positive_EHEC.py
```
#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="/home/jing/E.coli/blast_results/linkage/serotype/EHEC_esp_any_positive_ids.core.txt"
ASSEMBLY_ROOT="/home/jing/E.coli/blast_results/EHEC_assemblies"
OUTDIR="/home/jing/E.coli/blast_results/linkage/serotype/ectyper_out"
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
chmod +x run_ectyper_esp_positive.sh
```
```
cd /home/jing/E.coli/blast_results/linkage/serotype
./run_ectyper_esp_positive.sh
```

