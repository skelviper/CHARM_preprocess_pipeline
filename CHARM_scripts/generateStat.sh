#!/usr/bin/env bash
set -euo pipefail

if (( $# != 6 )); then
    printf 'usage: %s WORK_DIR INPUT_CONTRACT STRUCTURE_ENABLED EXPERIMENT_TYPE RESOLUTIONS_CSV RNA_COUNT_BAM\n' "$0" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "$1" && pwd)"
INPUT_CONTRACT=$2
STRUCTURE_ENABLED=$3
EXPERIMENT_TYPE=$4
RESOLUTIONS_CSV=$5
RNA_COUNT_BAM=$6
HELPER="${SCRIPT_DIR}/generate_stat_contract.py"

case "$STRUCTURE_ENABLED" in True|False) ;; *) printf 'invalid structure flag: %s\n' "$STRUCTURE_ENABLED" >&2; exit 2 ;; esac
case "$EXPERIMENT_TYPE" in charm|hires) ;; *) printf 'invalid experiment type: %s\n' "$EXPERIMENT_TYPE" >&2; exit 2 ;; esac
[[ -f "$INPUT_CONTRACT" ]] || { printf 'missing frozen input contract: %s\n' "$INPUT_CONTRACT" >&2; exit 2; }

cd "$WORK_DIR"
temporary_dir=$(mktemp -d "${WORK_DIR}/.generate-stat.XXXXXX")
trap 'rm -rf -- "$temporary_dir"' EXIT HUP INT TERM
temporary_stat="${temporary_dir}/stat"
manifest="${temporary_dir}/samples.tsv"
mkdir -p "$temporary_stat"

python "$HELPER" manifest \
    --contract "$INPUT_CONTRACT" \
    --work-dir "$WORK_DIR" \
    --output "$manifest"

gzip_line_count() {
    gzip -cd -- "$1" | awk 'END { print NR + 0 }'
}

gzip_nonheader_count() {
    gzip -cd -- "$1" | awk '$0 !~ /^#/ && NF { count++ } END { print count + 0 }'
}

require_file() {
    [[ -f "$1" ]] || { printf 'missing required QC input: %s\n' "$1" >&2; return 1; }
}

python "$HELPER" raw-fastq-stat \
    --contract "$INPUT_CONTRACT" \
    --work-dir "$WORK_DIR" \
    --workers "${CHARM_STAT_THREADS:-1}" \
    --output "$temporary_stat/raw.fq.stat"

: > "$temporary_stat/rna.fq.stat"
: > "$temporary_stat/dna.fq.stat"
: > "$temporary_stat/pairs.c1.stat"
: > "$temporary_stat/pairs.c12.stat"
while IFS=$'\t' read -r sample safe_code; do
    [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
    rna_r1="processed/${sample}/RNA/${sample}.rna.clean.R1.fq.gz"
    dna_r1="processed/${sample}/DNA/${sample}.dna.clean.R1.fq.gz"
    clean1="result/cleaned_pairs/c1/${sample}.pairs.gz"
    clean12="result/cleaned_pairs/c12/${sample}.pairs.gz"
    require_file "$rna_r1"
    require_file "$dna_r1"
    require_file "$clean1"
    require_file "$clean12"
    printf './%s\t%s\n' "$rna_r1" "$(gzip_line_count "$rna_r1")" >> "$temporary_stat/rna.fq.stat"
    printf './%s\t%s\n' "$dna_r1" "$(gzip_line_count "$dna_r1")" >> "$temporary_stat/dna.fq.stat"
    printf './%s\t%s\n' "$clean1" "$(gzip_nonheader_count "$clean1")" >> "$temporary_stat/pairs.c1.stat"
    printf './%s\t%s\n' "$clean12" "$(gzip_nonheader_count "$clean12")" >> "$temporary_stat/pairs.c12.stat"
done < "$manifest"

: > "$temporary_stat/raw.pairs.stat"
: > "$temporary_stat/pairs.dedup.stat"
while IFS=$'\t' read -r sample safe_code; do
    [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
    raw_pairs="processed/${sample}/2d_info/raw.pairs.gz"
    dedup_pairs="processed/${sample}/2d_info/contacts.pairs.gz"
    require_file "$raw_pairs"
    require_file "$dedup_pairs"
    printf './%s\t%s\n' "$raw_pairs" "$(gzip_nonheader_count "$raw_pairs")" >> "$temporary_stat/raw.pairs.stat"
    printf './%s\t%s\n' "$dedup_pairs" "$(gzip_nonheader_count "$dedup_pairs")" >> "$temporary_stat/pairs.dedup.stat"
done < "$manifest"

python "$HELPER" cleaned-pair-stats \
    --contract "$INPUT_CONTRACT" \
    --work-dir "$WORK_DIR" \
    --pairs-output "$temporary_stat/pairs.c123.stat" \
    --inter-output "$temporary_stat/inter.pairs.c123.stat" \
    --workers "${CHARM_STAT_THREADS:-1}"

if [[ "$STRUCTURE_ENABLED" == True ]]; then
    IFS=',' read -r -a structure_resolutions <<< "$RESOLUTIONS_CSV"
    (( ${#structure_resolutions[@]} > 0 )) || { printf 'enabled structure QC has no resolutions\n' >&2; exit 2; }
    python "$HELPER" rmsd-stats \
        --contract "$INPUT_CONTRACT" \
        --work-dir "$WORK_DIR" \
        --resolutions "${structure_resolutions[@]}" \
        --output "$temporary_stat/rmsd.info"
fi

: > "$temporary_stat/yperx.stat"
while IFS=$'\t' read -r sample safe_code; do
    [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
    yperx="processed/${sample}/2d_info/${sample}.yperx.txt"
    require_file "$yperx"
    value=$(<"$yperx")
    [[ -n "$value" ]] || { printf 'empty required yperx input: %s\n' "$yperx" >&2; exit 1; }
    printf '%s\t%s\n' "$yperx" "$value" >> "$temporary_stat/yperx.stat"
done < "$manifest"

for split in atac ct; do
    if [[ "$EXPERIMENT_TYPE" == charm ]]; then
        : > "$temporary_stat/${split}.read.stat"
        : > "$temporary_stat/${split}.dedup_rate.stat"
        while IFS=$'\t' read -r sample safe_code; do
            [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
            bam="processed/${sample}/${split}/${sample}.sort.bam"
            require_file "$bam"
            records=$(samtools flagstat "$bam" | awk 'NR == 1 { print $1; found=1 } END { if (!found) exit 1 }')
            printf '%s,%s\n' "$sample" "$records" >> "$temporary_stat/${split}.read.stat"
        done < "$manifest"

        while IFS=$'\t' read -r sample safe_code; do
            [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
            fragment_log="processed/${split}_all/${sample}.${split}.frag.log"
            require_file "$fragment_log"
            if grep -q 'Combined duplication rate' "$fragment_log"; then
                rate=$(sed -n 's/Combined duplication rate is \(.*\)%/\1/p' "$fragment_log")
                [[ $(printf '%s\n' "$rate" | awk 'NF { count++ } END { print count + 0 }') -eq 1 ]]
            elif grep -q 'Duplication rate is' "$fragment_log"; then
                rate=$(sed -n 's/Duplication rate is \(.*\)%/\1/p' "$fragment_log" | \
                    awk 'NR==1{r1=$1} NR==2{r2=$1} END{if(NR==2) printf "%.2f", (1-(1-r1/100)*(1-r2/100))*100; else if(NR==1) printf "%.2f", r1; else exit 1}')
            else
                rate=NA
            fi
            printf '%s\t%s\n' "$sample" "$rate" >> "$temporary_stat/${split}.dedup_rate.stat"
        done < "$manifest"
    fi
done

require_file "$RNA_COUNT_BAM"
samtools view "$RNA_COUNT_BAM" | python "$HELPER" rna-alignment-stat \
    --contract "$INPUT_CONTRACT" \
    --work-dir "$WORK_DIR" \
    --output "$temporary_stat/rna.reads_per_cell.stat"
[[ -s "$temporary_stat/rna.reads_per_cell.stat" ]] || {
    printf 'empty RNA per-cell statistics: %s\n' "$RNA_COUNT_BAM" >&2
    exit 1
}

: > "$temporary_stat/rna.dna_contam.stat"
while IFS=$'\t' read -r sample safe_code; do
    [[ "$sample" != sample || "$safe_code" != safe_code ]] || continue
    rna_r2="processed/${sample}/RNA/${sample}.rna.clean.R2.fq.gz"
    require_file "$rna_r2"
    gzip -cd -- "$rna_r2" | awk -v sample="$sample" \
        'NR%4==2 {total++; if ($0 ~ /TTTTTTTT[ACG][ACGT]GATC/) gatc++} END {print sample"\t"total+0"\t"gatc+0}' \
        >> "$temporary_stat/rna.dna_contam.stat"
done < "$manifest"

mkdir -p qc/stat
base_outputs=(
    raw.pairs.stat raw.fq.stat rna.fq.stat dna.fq.stat pairs.dedup.stat
    pairs.c1.stat pairs.c12.stat pairs.c123.stat inter.pairs.c123.stat
    yperx.stat rna.reads_per_cell.stat rna.dna_contam.stat
)
if [[ "$EXPERIMENT_TYPE" == charm ]]; then
    base_outputs+=(atac.read.stat ct.read.stat atac.dedup_rate.stat ct.dedup_rate.stat)
fi
if [[ "$STRUCTURE_ENABLED" == True ]]; then
    base_outputs+=(rmsd.info)
fi
for output_name in "${base_outputs[@]}"; do
    require_file "$temporary_stat/$output_name"
    mv -f -- "$temporary_stat/$output_name" "qc/stat/$output_name"
done
