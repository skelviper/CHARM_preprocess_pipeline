#!/usr/bin/env bash

# Record-count and valid empty-output helpers shared by Snakemake rules.

charm_gzip_nonheader_records() {
    gzip -cd -- "$1" | awk '$0 !~ /^#/ && NF { count++ } END { print count + 0 }'
}

charm_plain_nonheader_records() {
    awk '$0 !~ /^#/ && NF { count++ } END { print count + 0 }' "$1"
}

charm_copy_gzip_headers() {
    local input_path=$1
    local output_path=$2
    local compressor=${3:-gzip}

    case "$compressor" in
        gzip)
            gzip -cd -- "$input_path" | awk '/^#/' | gzip > "$output_path"
            ;;
        bgzip)
            gzip -cd -- "$input_path" | awk '/^#/' | bgzip > "$output_path"
            ;;
        *)
            return 2
            ;;
    esac
}

charm_make_empty_seg_from_bam() {
    local bam_path=$1
    local output_path=$2

    samtools view -H "$bam_path" | awk '
        $1 == "@SQ" {
            chromosome = ""
            sequence_length = ""
            for (field = 2; field <= NF; field++) {
                if (substr($field, 1, 3) == "SN:") {
                    chromosome = substr($field, 4)
                } else if (substr($field, 1, 3) == "LN:") {
                    sequence_length = substr($field, 4)
                }
            }
            if (chromosome != "" && sequence_length != "") {
                print "#chromosome: " chromosome " " sequence_length
                chromosomes++
            }
        }
        END { exit chromosomes > 0 ? 0 : 1 }
    ' | gzip > "$output_path"
}

charm_make_empty_pairs_from_seg() {
    local seg_path=$1
    local output_path=$2
    local compressor=${3:-gzip}

    case "$compressor" in
        gzip|bgzip) ;;
        *)
            return 2
            ;;
    esac

    gzip -cd -- "$seg_path" | awk '
        BEGIN {
            print "## pairs format v1.0"
            print "#sorted: chr1-chr2-pos1-pos2"
            print "#shape: upper triangle"
        }
        /^#chromosome:/ {
            print
            chromosomes++
        }
        END {
            if (chromosomes == 0) {
                exit 1
            }
            print "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1"
        }
    ' | "$compressor" > "$output_path"
}
