#!/bin/bash
###############################################################################
# CarpeDeam Pipeline — Ancient Chicken Cloaca Metagenomics (aDNA)
#
# Stages:
#   ref_prep_targeted → Download microbial reference genomes
#   merge_reads       → Merge paired reads (leeHom)
#   map_damage        → Map merged reads; compute damage; create .prof
#   assemble_carpedeam→ Run CarpeDeam ancient_assemble
#
# Run order:
# sbatch assembly_carpedeam.sh ref_prep_targeted
# sbatch --array=0-1 assembly_carpedeam.sh merge_reads
# sbatch --array=0-1 assembly_carpedeam.sh map_damage
# sbatch --array=0-1 assembly_carpedeam.sh assemble_carpedeam
###############################################################################

#SBATCH --job-name="carpedeam"
#SBATCH --time=48:00:00
#SBATCH -A gbru_fy23_ancient_chickens
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --error=%x.%A_%a.err
#SBATCH --output=%x.%A_%a.out

set -euo pipefail
export LC_ALL=C
export LANG=C

############################
# GLOBAL CONFIGURATION
############################
STAGE="${1:-}"
BASE="/project/gbru_fy23_ancient_chickens/annette"

# Resource & Tool Config
THREADS="${SLURM_CPUS_PER_TASK:-8}"
DP_HEAP="${DP_HEAP:-16g}"
DP_LEN="${DP_LEN:-100}"
ADAPTER_F="${ADAPTER_F:-}"
ADAPTER_R="${ADAPTER_R:-}"
LEEHOM_OPTS=("--ancientdna" "--prob" "0.8")

# Input/Output Paths
RAW_READS="$BASE/reads"
FILTERED_READS="$BASE/shortread/filtered_reads"
OUT_ROOT="$BASE/carpedeam"
REF_DIR="$OUT_ROOT/ref"
MERGE_DIR="$OUT_ROOT/merged"
MAP_DIR="$OUT_ROOT/mapping"
DMG_DIR="$OUT_ROOT/damage"
ASM_DIR="$OUT_ROOT/assembly"
TMP_DIR="$OUT_ROOT/tmp"

# Reference Files
REF_FASTA="$REF_DIR/targeted_microbes.fna"
REF_MMI="$REF_DIR/targeted_microbes.mmi"

# Samples (Slurm array indexing 0..1)
R1_FILES=(
  "$RAW_READS/Pen1-A10_S10_L001_R1_001.fastq.gz"
  "$RAW_READS/Pen9-B10_S22_L001_R1_001.fastq.gz"
)
R2_FILES=(
  "$RAW_READS/Pen1-A10_S10_L001_R2_001.fastq.gz"
  "$RAW_READS/Pen9-B10_S22_L001_R2_001.fastq.gz"
)

mkdir -p "$OUT_ROOT" "$REF_DIR" "$MERGE_DIR" "$MAP_DIR" "$DMG_DIR" "$ASM_DIR" "$TMP_DIR"

###############################################################################
# SAMPLE SETUP (Array-based)
###############################################################################
if [[ "$STAGE" != "ref_prep_targeted" && -n "$STAGE" ]]; then
  idx="${SLURM_ARRAY_TASK_ID:-0}"
  RAW_R1="${R1_FILES[$idx]}"
  SAMPLE="$(basename "$RAW_R1" _R1_001.fastq.gz)"

  # Resolve raw vs filtered reads
  F1="$FILTERED_READS/${SAMPLE}.clean_final1.fq.gz"
  F2="$FILTERED_READS/${SAMPLE}.clean_final2.fq.gz"
  if [[ -r "$F1" && -r "$F2" ]]; then
    USE_R1="$F1"; USE_R2="$F2"
  else
    USE_R1="$RAW_R1"; USE_R2="${R2_FILES[$idx]}"
    echo "[WARN] Filtered reads not found for $SAMPLE. Using RAW." >&2
  fi

  # Common Stage Paths
  DMG_PREFIX="$DMG_DIR/${SAMPLE}_CarpeDeam_"
  F5P="${DMG_PREFIX}5p.prof"
  F3P="${DMG_PREFIX}3p.prof"

  # Merged reads preference (prefer min20 filtered if it exists)
  MERGED_BASE="$MERGE_DIR/${SAMPLE}.merged"
  if [[ -s "${MERGED_BASE}.min20.fq.gz" ]]; then
    INPUT_MERGED="${MERGED_BASE}.min20.fq.gz"
  else
    INPUT_MERGED="${MERGED_BASE}.fq.gz"
  fi

  # Common output paths
  BAM="$MAP_DIR/${SAMPLE}.merged.sorted.bam"
  DP_OUT="$DMG_DIR/${SAMPLE}_damageprofiler"
  ASM_OUT="$ASM_DIR/${SAMPLE}.carpedeam.fasta"
  TMP_SUB="$TMP_DIR/${SAMPLE}"
fi

###############################################################################
# STAGE DISPATCH
###############################################################################
case "$STAGE" in

ref_prep_targeted)
  TOP_TSV="${2:-$REF_DIR/top_taxa_targets.tsv}"
  [[ ! -s "$TOP_TSV" ]] && { echo "[ERROR] top_taxa_targets.tsv not found at: $TOP_TSV" >&2; exit 1; }

  echo "[REF] Targeted ref prep from: $TOP_TSV"
  mkdir -p "$REF_DIR/targeted_ncbi" "$REF_DIR/logs"
  TLIST="$REF_DIR/targeted_ncbi/targeted_taxa.list"
  : > "$TLIST"

  module purge || true
  module load minimap2 samtools || true

  awk -F'\t' 'NR>1 && $8=="yes" && ($3=="species" || $3=="genus") {print $1"\t"$2"\t"$3}' "$TOP_TSV" > "$TLIST"
  echo "[REF] Selected $(wc -l < "$TLIST") taxa for download."

  CONDA_SH="/software/el9/apps/miniconda/24.7.1-2/etc/profile.d/conda.sh"
  [[ -f "$CONDA_SH" ]] && . "$CONDA_SH"
  USER_ENV="$OUT_ROOT/conda_env"
  if [[ ! -d "$USER_ENV" ]]; then
    conda create -y -p "$USER_ENV" -c conda-forge -c bioconda ncbi-datasets-cli seqkit
  fi

  download_taxon() {
    local taxid="$1" tname="$2" rank="$3"
    local query="${tname//_/ }"
    local zipfile="$REF_DIR/targeted_ncbi/${taxid}_${tname}.zip"
    local outdir="$REF_DIR/targeted_ncbi/${taxid}_${tname}"

    if conda run -p "$USER_ENV" datasets download genome taxon "$query" --assembly-source refseq --reference --dehydrated --filename "$zipfile" >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then :;
    elif conda run -p "$USER_ENV" datasets download genome taxon "$query" --assembly-source refseq --dehydrated --filename "$zipfile" >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then :;
    else return 1; fi

    mkdir -p "$outdir"
    unzip -q "$zipfile" -d "$outdir"
    conda run -p "$USER_ENV" datasets rehydrate --directory "$outdir" >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1
    find "$outdir/ncbi_dataset/data" -name "*_genomic.fna" -exec cat {} + >> "$REF_FASTA"
    return 0
  }

  : > "$REF_FASTA"
  while IFS=$'\t' read -r taxid name rank; do download_taxon "$taxid" "$name" "$rank" || true; done < "$TLIST"
  conda run -p "$USER_ENV" seqkit rmdup -s "$REF_FASTA" > "${REF_FASTA}.tmp" && mv "${REF_FASTA}.tmp" "$REF_FASTA"

  echo "[REF] Building Minimap2 index..."
  minimap2 -d "$REF_MMI" "$REF_FASTA" -I 4G
  ;;

merge_reads)
  echo "[MERGE] SAMPLE=$SAMPLE"
  module load miniconda || true
  source activate leehom_env || true

  MERGER_BIN="leeHom"
  command -v leeHomMulti >/dev/null 2>&1 && MERGER_BIN="leeHomMulti"

  OUT_PREFIX="$MERGE_DIR/${SAMPLE}.leehom"
  ADAPTER_OPTS=(--auto)
  [[ -n "$ADAPTER_F" && -n "$ADAPTER_R" ]] && ADAPTER_OPTS=(-f "$ADAPTER_F" -s "$ADAPTER_R")

  $MERGER_BIN "${LEEHOM_OPTS[@]}" -t "$THREADS" "${ADAPTER_OPTS[@]}" \
    -fq1 "$USE_R1" -fq2 "$USE_R2" -fqo "$OUT_PREFIX"

  ln -sf "${OUT_PREFIX}.fq.gz"    "${MERGE_DIR}/${SAMPLE}.merged.fq.gz"
  ln -sf "${OUT_PREFIX}_r1.fq.gz" "${MERGE_DIR}/${SAMPLE}.unmerged_R1.fq.gz"
  ln -sf "${OUT_PREFIX}_r2.fq.gz" "${MERGE_DIR}/${SAMPLE}.unmerged_R2.fq.gz"

  if command -v seqkit >/dev/null 2>&1; then
    seqkit seq -m 20 "${MERGE_DIR}/${SAMPLE}.merged.fq.gz" | gzip > "${MERGE_DIR}/${SAMPLE}.merged.min20.fq.gz"
  fi

  if [[ ! -s "${MERGE_DIR}/${SAMPLE}.merged.fq.gz" ]]; then
    echo "[ERROR] leeHom failed to create ${MERGE_DIR}/${SAMPLE}.merged.fq.gz" >&2
    exit 1
  fi
  ;;

map_damage)
  echo "[MAP+DMG] SAMPLE=$SAMPLE"
  [[ -s "$REF_MMI" ]] || { echo "[ERROR] Index missing"; exit 1; }
  module load minimap2 samtools miniconda || true
  source activate carpedeam_env

  minimap2 -a -x sr --MD -t "$THREADS" "$REF_MMI" "$INPUT_MERGED" | samtools sort -@ "$THREADS" -o "$BAM"
  samtools index "$BAM"

  mkdir -p "$DP_OUT"
  damageprofiler -i "$BAM" -r "$REF_FASTA" -o "$DP_OUT" -l "$DP_LEN"

  # Convert DP to CarpeDeam .prof
  export FREQ5="$DP_OUT/5pCtoT_freq.txt"
  export FREQ3="$DP_OUT/3pGtoA_freq.txt"
  export F5P F3P
  python - <<'PY'
import os, sys
subs = ["A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"]
def read_vals(path):
    vals=[]
    with open(path,'r') as f:
        for line in f:
            if line.lower().startswith(('pos','position','#')) or not line.strip(): continue
            vals.append(line.strip().split('\t')[-1])
    return vals
def write_prof(path, vals, col):
    with open(path,'w') as f:
        f.write("\t".join(subs) + "\n")
        for v in vals:
            row=['0.0']*12; row[col]=v
            f.write("\t".join(row) + "\n")
v5 = read_vals(os.environ['FREQ5']); v3 = read_vals(os.environ['FREQ3'])
write_prof(os.environ['F5P'], v5, 5); write_prof(os.environ['F3P'], v3, 6)
PY

  # Streamlined Normalization
  [[ -s "$F5P" && -s "$F3P" ]] || { echo "[ERROR] .prof files missing"; exit 1; }
  sed -E 's/[[:space:]]+/\t/g; s/\t$//' -i "$F5P" "$F3P"
  ;;

assemble_carpedeam)
  echo "[ASM] SAMPLE=$SAMPLE"
  module load miniconda || true
  source activate carpedeam_env

  mkdir -p "$TMP_SUB"

  carpedeam ancient_assemble "$INPUT_MERGED" "$ASM_OUT" "$TMP_SUB" \
    --ancient-damage "$DMG_PREFIX" --num-iter-reads-only 5 --num-iterations 12 \
    --min-merge-seq-id 95 --min-cov-safe 2
  ;;

*)
  echo "Usage: $0 <stage>"
  exit 1
  ;;
esac
