#!/bin/bash
###############################################################################
# CarpeDeam Pipeline — Ancient Chicken Cloaca Metagenomics (aDNA)
#
# Stages:
#   ref_prep_targeted → Download microbial reference genomes based on previous SPAdes assembly and taxonomy
#   merge_reads       → Merge paired reads (BBMerge), optionally drop ultra-short merges
#   map_damage        → Map merged reads → BAM; DamageProfiler → CarpeDeam *.prof damage matrices
#   assemble_carpedeam→ Run CarpeDeam ancient_assemble with merged reads and damage matrix
#
# Notes:
# - CarpeDeam requires merged reads + a strict tab-separated damage matrix (two files: *_5p.prof, *_3p.prof).  # See README. [7](https://link.springer.com/content/pdf/10.1186/s13059-025-03839-5.pdf)
# - DamageProfiler computes 5′C→T and 3′G→A profiles on mapped reads; we convert these into CarpeDeam’s 12-column .prof format. [6](https://usdagcc-my.sharepoint.com/personal/annette_hynes_usda_gov/Documents/Microsoft%20Copilot%20Chat%20Files/Pen1-A10_S10_k2_khist.txt)
# - Bowtie2 is used for single-end mapping of merged reads; index built once in ref_prep. [4](https://www.biostars.org/p/208696/)[5](https://www.biostars.org/p/485916/)
# - BBMerge is used for overlap-based merging. [8](https://www.anaconda.com/docs/getting-started/working-with-conda/environments)[9](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html)
#
# Run order:
# sbatch assembly_carpedeam.sh ref_prep_targeted
# sbatch --array=0-1 assembly_carpedeam.sh merge_reads
# sbatch --array=0-1 assembly_carpedeam.sh map_damage
# sbatch --array=0-1 assembly_carpedeam.sh assemble_carpedeam
###############################################################################

############################
# SLURM CONFIGURATION
############################
#SBATCH --job-name="carpedeam"
#SBATCH --time=48:00:00
#SBATCH -A gbru_fy23_ancient_chickens
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --error=%x.%A_%a.err
#SBATCH --output=%x.%A_%a.out

set -euo pipefail

############################
# GLOBAL PATHS
############################
STAGE="${1:-}"

BASE="/project/gbru_fy23_ancient_chickens/annette"

# Input reads
RAW_READS="$BASE/reads"
FILTERED_READS="$BASE/shortread/filtered_reads"   # from your host_filter stage (nochicken/nohuman → clean_final*.fq.gz)

# Pipeline outputs
OUT_ROOT="$BASE/carpedeam"
REF_DIR="$OUT_ROOT/ref"
MERGE_DIR="$OUT_ROOT/merged"
MAP_DIR="$OUT_ROOT/mapping"
DMG_DIR="$OUT_ROOT/damage"
ASM_DIR="$OUT_ROOT/assembly"
TMP_DIR="$OUT_ROOT/tmp"

mkdir -p "$OUT_ROOT" "$REF_DIR" "$MERGE_DIR" "$MAP_DIR" "$DMG_DIR" "$ASM_DIR" "$TMP_DIR"

# Reference outputs
REF_FASTA="$REF_DIR/microbes_reference.fna"
REF_INDEX_PREFIX="$REF_DIR/microbes_index"

# Samples (Slurm array indexing 0..1)
R1_FILES=(
  "$RAW_READS/Pen1-A10_S10_L001_R1_001.fastq.gz"
  "$RAW_READS/Pen9-B10_S22_L001_R1_001.fastq.gz"
)
R2_FILES=(
  "$RAW_READS/Pen1-A10_S10_L001_R2_001.fastq.gz"
  "$RAW_READS/Pen9-B10_S22_L001_R2_001.fastq.gz"
)

THREADS="${SLURM_CPUS_PER_TASK:-8}"

###############################################################################
# Helper: resolve filtered vs raw paired reads for a SAMPLE
###############################################################################
reads_for_sample() {
  local sample="$1"
  local f1="$FILTERED_READS/${sample}.clean_final1.fq.gz"
  local f2="$FILTERED_READS/${sample}.clean_final2.fq.gz"

  if [[ -r "$f1" && -r "$f2" ]]; then
    # Print each path on its own line, exactly as two separate records
    printf '%s\n' "$f1" "$f2"
  else
    local idx="${SLURM_ARRAY_TASK_ID:-0}"
    local r1="${R1_FILES[$idx]}"
    local r2="${R2_FILES[$idx]}"
    printf '%s\n' "$r1" "$r2"
    echo "[WARN] Filtered reads not found for $sample. Using RAW." >&2
  fi
}

###############################################################################
# STAGE DISPATCH
###############################################################################
case "$STAGE" in

# ---------------------------------------------------------------------------
# REF PREP (TARGETED) — use top_taxa_targets.tsv → download genomes → build Minimap2 index
# ---------------------------------------------------------------------------
ref_prep_targeted)
  set -euo pipefail
  TOP_TSV="${2:-$REF_DIR/top_taxa_targets.tsv}"
  [[ ! -s "$TOP_TSV" ]] && { echo "[ERROR] top_taxa_targets.tsv not found at: $TOP_TSV" >&2; exit 1; }

  echo "[REF] Targeted ref prep from: $TOP_TSV"
  mkdir -p "$REF_DIR/targeted_ncbi" "$REF_DIR/logs"
  TLIST="$REF_DIR/targeted_ncbi/targeted_taxa.list"
  : > "$TLIST"

  # 0) Module hygiene: avoid python_3 conflicts, keep toolchain minimal
  module purge || true
  module load minimap2 || true
  module load samtools || true

  # 1) Select only species/genus rows with include=yes (skip higher ranks)
  awk -F'\t' 'NR>1 && $8=="yes" && ($3=="species" || $3=="genus") {print $1"\t"$2"\t"$3}' "$TOP_TSV" > "$TLIST"
  NSEL=$(wc -l < "$TLIST" || echo 0)
  [[ "$NSEL" -eq 0 ]] && { echo "[ERROR] No species/genus marked include=yes in $TOP_TSV"; exit 1; }
  echo "[REF] Selected $NSEL taxa (species/genus) for download."

  # 2) Use a user-writable conda env under the project; avoid system base
  CONDA_SH="/software/el9/apps/miniconda/24.7.1-2/etc/profile.d/conda.sh"
  if [[ -f "$CONDA_SH" ]]; then
    . "$CONDA_SH"
  else
    echo "[WARN] conda.sh not found; ensure 'conda' is in PATH."
  fi

  USER_ENV="$OUT_ROOT/conda_env"
  if [[ ! -d "$USER_ENV" ]]; then
    echo "[REF] Creating local env: $USER_ENV"
    conda create -y -p "$USER_ENV" -c conda-forge -c bioconda ncbi-datasets-cli seqkit || {
      echo "[ERROR] Could not create local Conda env at $USER_ENV"; exit 1; }
  else
    echo "[REF] Using existing env: $USER_ENV"
  fi

  # 3) Robust downloader: try RefSeq reps → RefSeq → GenBank reps → GenBank
  download_taxon() {
    local taxid="$1" tname="$2" rank="$3"
    local query="${tname//_/ }"   # underscores → spaces
    local outbase="$REF_DIR/targeted_ncbi/${taxid}_${tname}"
    local zipfile="$outbase.zip"
    local outdir="$outbase"

    # Try 1: RefSeq representatives (preferred; compact)
    if conda run -p "$USER_ENV" datasets download genome taxon "$query" \
         --assembly-source refseq --reference --dehydrated --filename "$zipfile" \
         >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then
      :
    # Try 2: RefSeq (no representative filter)
    elif conda run -p "$USER_ENV" datasets download genome taxon "$query" \
         --assembly-source refseq --dehydrated --filename "$zipfile" \
         >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then
      :
    # Try 3: GenBank representatives
    elif conda run -p "$USER_ENV" datasets download genome taxon "$query" \
         --assembly-source genbank --reference --dehydrated --filename "$zipfile" \
         >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then
      :
    # Try 4: GenBank (no representative filter)
    elif conda run -p "$USER_ENV" datasets download genome taxon "$query" \
         --assembly-source genbank --dehydrated --filename "$zipfile" \
         >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1; then
      :
    else
      echo "[WARN] No assemblies for $rank $tname (taxid=$taxid). Skipping."
      return 1
    fi

    [[ ! -s "$zipfile" ]] && { echo "[WARN] Empty package for $tname (taxid=$taxid)."; return 1; }

    mkdir -p "$outdir"
    unzip -q "$zipfile" -d "$outdir"
    conda run -p "$USER_ENV" datasets rehydrate --directory "$outdir" \
      >> "$REF_DIR/logs/datasets_${taxid}.log" 2>&1

    # Append *_genomic.fna files to combined FASTA (Datasets package path)
    local nfa=0
    while IFS= read -r -d '' fna; do
      cat "$fna" >> "$REF_DIR/targeted_microbes.fna"
      ((nfa++))
    done < <(find "$outdir/ncbi_dataset/data" -name "*_genomic.fna" -print0)

    if [[ "$nfa" -eq 0 ]]; then
      echo "[WARN] No genomic FASTA found for $tname (taxid=$taxid)."
      return 1
    else
      echo "[OK] $tname: appended $nfa FASTA(s)."
    fi
    return 0
  }

  # 4) Build combined FASTA
  : > "$REF_DIR/targeted_microbes.fna"
  echo "[REF] Downloading & concatenating genomic FASTAs..."
  while IFS=$'\t' read -r taxid name rank; do
    download_taxon "$taxid" "$name" "$rank" || true
  done < "$TLIST"

  # Optional: remove exact duplicate sequences (if seqkit is available)
  if conda run -p "$USER_ENV" which seqkit >/dev/null 2>&1; then
    echo "[REF] Removing exact duplicates with seqkit rmdup..."
    conda run -p "$USER_ENV" seqkit rmdup -s "$REF_DIR/targeted_microbes.fna" > "$REF_DIR/targeted_microbes.dedup.fna"
    mv -f "$REF_DIR/targeted_microbes.dedup.fna" "$REF_DIR/targeted_microbes.fna"
  fi

  NSEQ=$(grep -c "^>" "$REF_DIR/targeted_microbes.fna" || echo 0)
  [[ "$NSEQ" -eq 0 ]] && { echo "[ERROR] Combined FASTA is empty."; exit 1; }
  echo "[REF] Combined FASTA has $NSEQ sequences: $REF_DIR/targeted_microbes.fna"

  # 5) Build Minimap2 index (OOM-safe batching) [2](https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2)
  echo "[REF] Building Minimap2 index (batched, -I 4G)..."
  minimap2 -d "$REF_DIR/targeted_microbes.mmi" "$REF_DIR/targeted_microbes.fna" -I 4G
  echo "[REF] Index built: $REF_DIR/targeted_microbes.mmi"

  echo "[REF] Targeted reference prepared."
  ;;

#############################################################################
# MERGE READS — overlap-based merging with BBMerge (robust twin-file inputs)
#############################################################################
merge_reads)
  set -euo pipefail
  module load bbtools || true  # bbmerge.sh wrapper lives in the BBTools module
  mkdir -p "$MERGE_DIR"

  # SLURM array index and sample name
  idx="${SLURM_ARRAY_TASK_ID:-0}"
  R1="${R1_FILES[$idx]}"
  SAMPLE="$(basename "$R1" _R1_001.fastq.gz)"

  # Resolve filtered vs raw paired reads for this sample:
  # reads_for_sample prints R1 on line1 and R2 on line2. We read both lines safely.
  mapfile -t __PAIR < <(reads_for_sample "$SAMPLE")
  USE_R1="${__PAIR-}"
  USE_R2="${__PAIR[1]:-}"

  # Hard stop if either is empty or unreadable
  if [[ -z "$USE_R1" || -z "$USE_R2" || ! -r "$USE_R1" || ! -r "$USE_R2" ]]; then
    printf '[ERROR] Cannot read one or both input files:\n  R1=%qR2=%q\n' "$USE_R1" "$USE_R2" >&2
    exit 1
  fi

  printf '[MERGE] R1: %q\n[MERGE] R2: %q\n' "$USE_R1" "$USE_R2"

  # Outputs
  MERGED="$MERGE_DIR/${SAMPLE}.merged.fq.gz"
  UNM_R1="$MERGE_DIR/${SAMPLE}.unmerged_R1.fq.gz"
  UNM_R2="$MERGE_DIR/${SAMPLE}.unmerged_R2.fq.gz"
  IHIST="$MERGE_DIR/${SAMPLE}.insert_hist.txt"

  # Overlap-based merging (preset 'loose' as a flag; not 'strictness=loose')
  # For twin files use: in1=<R1> in2=<R2>. Preset flags are documented in BBMerge guide. [1](https://anaconda.org/conda-forge/ncbi-datasets-cli)
  bbmerge.sh \
    in1="$USE_R1" in2="$USE_R2" \
    out="$MERGED" outu1="$UNM_R1" outu2="$UNM_R2" ihist="$IHIST" \
    loose=t t="${THREADS}"

  # Optional: drop ultra-short merges (<20 bp). Helps downstream tools avoid 1–19 bp artifacts.
  if command -v seqkit >/dev/null 2>&1; then
    echo "[MERGE] Filtering merged reads <20bp..."
    seqkit seq -m 20 "$MERGED" > "$MERGE_DIR/${SAMPLE}.merged.min20.fq"
    gzip -f "$MERGE_DIR/${SAMPLE}.merged.min20.fq"
    MERGED="$MERGE_DIR/${SAMPLE}.merged.min20.fq.gz"
  fi

  # Minimal summary
  if [[ -s "$IHIST" ]]; then
    echo "[MERGE] Insert-size histogram (head):"
    head -n 5 "$IHIST" | sed 's/\t/  /g'
  fi
  echo "[MERGE] Done: $MERGED"
  ;;

#############################################
# MAP DAMAGE: mapping → DamageProfiler → .prof
#############################################
map_damage)
  set -euo pipefail

  # --- Config knobs (override at submit time) ---
  DP_HEAP="${DP_HEAP:-16g}"   # JVM heap for DamageProfiler
  DP_LEN="${DP_LEN:-100}"     # Number of bases for frequency computations (DP default is 100)
  THREADS="${SLURM_CPUS_PER_TASK:-8}"

  mkdir -p "$MAP_DIR" "$DMG_DIR" "$TMP_DIR"

  # --- Verify targeted reference & index ---
  REF_FASTA_CAND=("$REF_DIR/targeted_microbes.fna" "$REF_FASTA")
  REF_MMI="$REF_DIR/targeted_microbes.mmi"
  REF_FASTA_FOUND=""
  for f in "${REF_FASTA_CAND[@]}"; do
    [[ -s "$f" ]] && REF_FASTA_FOUND="$f" && break
  done
  [[ -z "$REF_FASTA_FOUND" ]] && { echo "[ERROR] Reference FASTA not found."; exit 1; }
  [[ ! -s "$REF_MMI" ]] && { echo "[ERROR] Minimap2 index missing."; exit 1; }

  # --- Environment ---
  module load minimap2   || true
  module load samtools   || true
  module load miniconda  || true
  source activate carpedeam_env  # must contain DamageProfiler + Java ≥11

  # --- Resolve sample + merged reads ---
  idx="${SLURM_ARRAY_TASK_ID:-0}"
  R1="${R1_FILES[$idx]}"
  SAMPLE="$(basename "$R1" _R1_001.fastq.gz)"
  MERGED_CAND=("$MERGE_DIR/${SAMPLE}.merged.min20.fq.gz" "$MERGE_DIR/${SAMPLE}.merged.fq.gz")
  MERGED=""
  for f in "${MERGED_CAND[@]}"; do
    [[ -s "$f" ]] && MERGED="$f" && break
  done
  [[ -z "$MERGED" ]] && { echo "[ERROR] Merged reads not found for $SAMPLE"; exit 1; }

  # --- Map with Minimap2 → position-sorted BAM ---
  BAM="$MAP_DIR/${SAMPLE}.merged.sorted.bam"
  echo "[MAP] Minimap2 mapping (-x sr) → sorted BAM..."
  minimap2 -a -x sr -t "$THREADS" "$REF_MMI" "$MERGED" \
    | samtools sort -@ "$THREADS" -o "$BAM"
  samtools index "$BAM"

  # --- Contigs with any mapped reads (species/header list for DP) ---
  SPECIES_FILE="$DMG_DIR/${SAMPLE}_species.txt"
  samtools idxstats "$BAM" \
    | awk -F'\t' '$3 > 0 {print $1}' \
    | sort -u > "$SPECIES_FILE"

  # --- Reduced FASTA (mapped contigs only) to minimize memory ---
  REDUCED_FASTA="$REF_DIR/targeted_microbes.hitonly.${SAMPLE}.fna"
  if [[ -s "$SPECIES_FILE" ]]; then
    samtools faidx "$REF_FASTA_FOUND" -r "$SPECIES_FILE" > "$REDUCED_FASTA"
    samtools faidx "$REDUCED_FASTA"
  else
    # No hits found: fall back to full FASTA but still indexed
    cp "$REF_FASTA_FOUND" "$REDUCED_FASTA"
    [[ ! -s "${REF_FASTA_FOUND}.fai" ]] && samtools faidx "$REF_FASTA_FOUND"
    samtools faidx "$REDUCED_FASTA"
  fi

  # --- DamageProfiler (always with reference + species file) ---
  echo "[DMG] Running DamageProfiler with reference (hit-only FASTA)..."
  DP_OUT="$DMG_DIR/${SAMPLE}_damageprofiler"
  mkdir -p "$DP_OUT"

  export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx${DP_HEAP} -Xms4g"
  if [[ -s "$SPECIES_FILE" ]]; then
    damageprofiler -i "$BAM" -r "$REDUCED_FASTA" -sf "$SPECIES_FILE" -o "$DP_OUT" -l "$DP_LEN"
  else
    damageprofiler -i "$BAM" -r "$REDUCED_FASTA" -o "$DP_OUT" -l "$DP_LEN"
  fi

  # Validate DP outputs exist (these are the ones we convert)
  [[ -s "$DP_OUT/5pCtoT_freq.txt" && -s "$DP_OUT/3pGtoA_freq.txt" ]] \
    || { echo "[ERROR] DamageProfiler did not produce freq tables."; exit 1; }

  conda deactivate || true

  # --- Convert DP outputs → CarpeDeam .prof ---
  echo "[DMG] Converting DamageProfiler output to CarpeDeam .prof files..."
  FREQ5="$DP_OUT/5pCtoT_freq.txt"
  FREQ3="$DP_OUT/3pGtoA_freq.txt"
  DMG_PREFIX="$DMG_DIR/${SAMPLE}_CarpeDeam"
  F5P="${DMG_PREFIX}_5p.prof"
  F3P="${DMG_PREFIX}_3p.prof"

  # Ensure env vars are available to the Python block
  export FREQ5 FREQ3 F5P F3P

  python - <<'PY'
import os, sys
subs = ["A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"]

def need(name):
    v = os.environ.get(name)
    if v is None or not v:
        sys.stderr.write(f"[ERROR] Missing environment variable: {name}\n")
        sys.exit(2)
    return v

def read_vals(path):
    vals=[]
    try:
        with open(path,'r',encoding='utf-8',errors='replace') as fh:
            for line in fh:
                s=line.strip()
                if not s or s.lower().startswith(('pos','position','#')):
                    continue
                parts = s.split('\t')
                try:
                    vals.append(float(parts[-1].strip()))
                except Exception:
                    try:
                        vals.append(float(s.split()[-1]))
                    except Exception:
                        pass
    except FileNotFoundError:
        sys.stderr.write(f"[ERROR] Missing DP file: {path}\n")
        sys.exit(2)
    return vals

def write_prof(path, vals, which_col):
    with open(path,'w',encoding='ascii') as out:
        out.write("\t".join(subs) + "\n")
        for v in vals:
            row=[0.0]*12
            row[which_col]=v
            out.write("\t".join(f"{x:.10f}" for x in row) + "\n")

f5p = need('F5P'); f3p = need('F3P')
freq5 = need('FREQ5'); freq3 = need('FREQ3')
vals5 = read_vals(freq5)
vals3 = read_vals(freq3)

# Always require both profiles to be created
if not vals5:
    sys.stderr.write("[ERROR] No 5' C>T values read from DamageProfiler.\n")
    sys.exit(2)
if not vals3:
    sys.stderr.write("[ERROR] No 3' G>A values read from DamageProfiler.\n")
    sys.exit(2)

write_prof(f5p, vals5, 5)  # C>T (5' file)
write_prof(f3p, vals3, 6)  # G>A (3' file)
PY

  # Validate that both files exist and are non-empty before normalization
  for prof in "$F5P" "$F3P"; do
    if [[ -z "${prof:-}" ]]; then
      echo "[ERROR] Empty .prof path variable." >&2
      exit 1
    fi
    if [[ ! -s "$prof" ]]; then
      echo "[ERROR] .prof file not created or empty: $prof" >&2
      exit 1
    fi
  done

  # --- Normalize and validate .prof (strict tabs; 12 fields) ---
  normalize_prof() {
    local prof="${1:-}"
    [[ -z "$prof" ]] && { echo "[ERROR] normalize_prof: missing path argument." >&2; return 1; }
    local tmp="${prof}.tmp"
    sed -E '1s/^\xEF\xBB\xBF//; s/[[:space:]]+/\t/g; s/\t$//' "$prof" | tr -d '\r' > "$tmp"
    mv -f "$tmp" "$prof"
  }

  validate_prof() {
    local prof="${1:-}"
    [[ -z "$prof" ]] && { echo "[ERROR] validate_prof: missing path argument." >&2; return 1; }
    local hdr_nf
    hdr_nf=$(awk -F'\t' 'NR==1{print NF; exit}' "$prof")
    [[ "${hdr_nf:-0}" -ne 12 ]] && { echo "[ERROR] Header NF=${hdr_nf} (expected 12): $prof" >&2; return 1; }
    awk -F'\t' 'NR>1 && NF!=12{printf("[ERROR] %s: BAD line %d (NF=%d)\n", FILENAME, NR, NF); exit 1}' "$prof" || return 1
  }

  # Apply to both profiles
  for prof in "$F5P" "$F3P"; do
    normalize_prof "$prof"
    validate_prof "$prof"
  done

  echo "[MAP+DMG] Stage complete for $SAMPLE."
  ;;

###############################################################################
# ASSEMBLY — CarpeDeam ancient_assemble
###############################################################################
assemble_carpedeam)
  set -euo pipefail
  module load miniconda || true
  mkdir -p "$ASM_DIR" "$TMP_DIR"

  # Resolve sample and merged reads (prefers min20 if present)
  idx="${SLURM_ARRAY_TASK_ID:-0}"
  R1="${R1_FILES[$idx]}"
  SAMPLE="$(basename "$R1" _R1_001.fastq.gz)"

  MERGED_CAND=(
    "$MERGE_DIR/${SAMPLE}.merged.min20.fq.gz"
    "$MERGE_DIR/${SAMPLE}.merged.fq.gz"
  )
  MERGED=""
  for f in "${MERGED_CAND[@]}"; do
    [[ -s "$f" ]] && MERGED="$f" && break
  done
  if [[ -z "$MERGED" ]]; then
    echo "[ERROR] Merged reads not found for $SAMPLE" >&2
    echo "  Checked:"
    printf '    %s\n' "${MERGED_CAND[@]}"
    exit 1
  fi

  # Damage matrix prefix (CarpeDeam will append _5p.prof and _3p.prof)
  DMG_PREFIX="$DMG_DIR/${SAMPLE}_CarpeDeam"
  F5P="${DMG_PREFIX}_5p.prof"
  F3P="${DMG_PREFIX}_3p.prof"
  if [[ ! -s "$F5P" || ! -s "$F3P" ]]; then
    echo "[ERROR] Damage matrices missing for $SAMPLE. Run map_damage first." >&2
    echo "  Expected:"
    echo "    $F5P"
    echo "    $F3P"
    exit 1
  fi

  # Preflight: normalize headers/whitespace/CRs; re-validate strict 12-column tab separation.
  # (This mirrors the CarpeDeam README/FAQ suggestion to fix 'Profile not 12 fields' failures.)
  normalize_prof() {
    local prof="$1"
    local tmp="${prof}.tmp"
    # Normalize: strip BOM, collapse whitespace→TAB, strip trailing TAB, remove CRs.
    sed -E '1s/^\xEF\xBB\xBF//; s/[[:space:]]+/\t/g; s/\t$//' "$prof" | tr -d '\r' > "$tmp"
    mv -f "$tmp" "$prof"
  }
  validate_prof() {
    local prof="$1"
    local bad=0
    local hdr_nf
    hdr_nf=$(awk -F'\t' 'NR==1{print NF; exit}' "$prof")
    if [[ "${hdr_nf:-0}" -ne 12 ]]; then
      echo "[ERROR] Header has NF=${hdr_nf} (expected 12) in: $prof" >&2
      return 1
    fi
    awk -F'\t' 'NR>1 && NF!=12{printf("[ERROR] %s: BAD line %d (NF=%d)\n", FILENAME, NR, NF); exit 1}' "$prof" || true
  }

  echo "[ASM] Normalizing and validating damage matrices for $SAMPLE ..."
  normalize_prof "$F5P"
  normalize_prof "$F3P"
  validate_prof "$F5P"
  validate_prof "$F3P"

  # Print a concise summary (full paths + row counts)
  rows_5p=$(( $(wc -l < "$F5P") - 1 ))
  rows_3p=$(( $(wc -l < "$F3P") - 1 ))
  echo "[ASM] Damage matrices:"
  echo "  5p: $F5P  (rows excluding header: $rows_5p)"
  echo "  3p: $F3P  (rows excluding header: $rows_3p)"

  # Prepare output and tmp directories (full paths)
  ASM_OUT="$ASM_DIR/${SAMPLE}.carpedeam.fasta"
  TMP_SUB="$TMP_DIR/${SAMPLE}"
  mkdir -p "$TMP_SUB"

  # Activate environment and run CarpeDeam (capture detailed logs)
  echo "[ASM] Running CarpeDeam..."
  echo "  Command:"
  echo "    carpedeam ancient_assemble \\"
  echo "      $MERGED \\"
  echo "      $ASM_OUT \\"
  echo "      $TMP_SUB \\"
  echo "      --ancient-damage $DMG_PREFIX \\"
  echo "      --num-iter-reads-only 4 \\"
  echo "      --num-iterations 10 \\"
  echo "      --min-merge-seq-id 99"

  source activate carpedeam_env
  carpedeam ancient_assemble \
    "$MERGED" \
    "$ASM_OUT" \
    "$TMP_SUB" \
    --ancient-damage "$DMG_PREFIX" \
    --num-iter-reads-only 4 \
    --num-iterations 10 \
    --min-merge-seq-id 99 \
    1> "$ASM_DIR/${SAMPLE}.carpedeam.out" \
    2> "$ASM_DIR/${SAMPLE}.carpedeam.err"
  conda deactivate || true

  echo "[ASM] Done: $ASM_OUT"
  echo "[ASM] Logs:"
  echo "  stdout: $ASM_DIR/${SAMPLE}.carpedeam.out"
  echo "  stderr: $ASM_DIR/${SAMPLE}.carpedeam.err"
  ;;

###############################################################################
# Default
###############################################################################
*)
  echo "Usage: sbatch carpedeam_pipeline.sh <stage>"
  echo "Stages:"
  echo "  ref_prep_targeted   # Prepare reference FASTA and minimap index (run once)"
  echo "  merge_reads         # Merge paired reads (array job)"
  echo "  map_damage          # Map merged reads; compute damage; create CarpeDeam .prof (array job)"
  echo "  assemble_carpedeam  # Run CarpeDeam ancient_assemble (array job)"
  exit 1
  ;;
esac
