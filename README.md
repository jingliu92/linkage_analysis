# linkage_analysis

## Create marker gene file
nano all.gene.fasta

## make db
```
makeblastdb -in all_markers.fasta -dbtype nucl -out all_markers
```
## Blast
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
## count all espK/espV/espN combinations with in EHEC
```
#!/usr/bin/env python3
import pandas as pd

# 1. Read matrix
df = pd.read_csv("marker_presence_absence.tsv", sep="\t")

genes = ["espK", "espV", "espN"]

# 2. Make a nice combo label like espK+_espV-_espN+
def make_label(row):
    parts = []
    for g in genes:
        val = row[g]
        sign = "+" if val == 1 else "-"
        parts.append(f"{g}{sign}")
    return "_".join(parts)

df["combo"] = df.apply(make_label, axis=1)

# 3. Overall counts of all combinations
combo_counts = df["combo"].value_counts().reset_index()
combo_counts.columns = ["Combination", "N_isolates"]
print("Overall espK/espV/espN combinations:")
print(combo_counts)

# 4. If you have a Pathovar column, stratify by EHEC vs others
if "Pathovar" in df.columns:
    table = (
        df.groupby(["Pathovar", "combo"])
          .size()
          .reset_index(name="N_isolates")
          .pivot(index="combo", columns="Pathovar", values="N_isolates")
          .fillna(0)
          .astype(int)
    )
    print("\nCombination counts by Pathovar:")
    print(table)

    # Save tables
    combo_counts.to_csv("espK_espV_espN_combos_overall.tsv", sep="\t", index=False)
    table.to_csv("espK_espV_espN_combos_by_pathovar.tsv", sep="\t")
else:
    combo_counts.to_csv("espK_espV_espN_combos_overall.tsv", sep="\t", index=False)
```
