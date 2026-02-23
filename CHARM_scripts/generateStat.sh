#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_FILE="${PIPELINE_DIR}/config.yaml"

WORK_DIR="$(awk -F': *' '/^work_dir:/ {print $2; exit}' "${CONFIG_FILE}" | tr -d '\r')"
if [ -z "${WORK_DIR}" ]; then
  WORK_DIR=".."
fi
if [[ "${WORK_DIR}" != /* ]]; then
  WORK_DIR="${PIPELINE_DIR}/${WORK_DIR}"
fi
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

cd "${WORK_DIR}"
mkdir ./stat

find ./ -name "raw.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#"| wc -l' | sort > ./stat/raw.pairs.stat
find -L ./Rawdata/ -name "*1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/raw.fq.stat
find ./processed -name "*.rna.clean.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/rna.fq.stat
find ./processed -name "*.dna.clean.R1.fq.gz" | parallel --tag 'zcat {} | wc -l' | sort > ./stat/dna.fq.stat
find ./processed -name "contacts.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.dedup.stat
find ./result/cleaned_pairs/c1 -name "*.pairs.gz" | parallel --tag 'zcat {} |grep -v "^#" | wc -l' | sort > ./stat/pairs.c1.stat
find ./result/cleaned_pairs/c12 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c12.stat
find ./result/cleaned_pairs/c123 -name "*.pairs.gz" | parallel --tag 'zcat {} | grep -v "^#" |wc -l' | sort > ./stat/pairs.c123.stat
for i in `find ./result/cleaned_pairs/c123 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c123.stat

for i in `find ./result/cleaned_pairs/c123 -name "*.pairs.gz" |sort`; do echo -n $i;echo -n -e "\t";zcat $i| grep -v "^#" | awk '$2!=$4 {print $0}' | wc -l; done > ./stat/inter.pairs.c123.stat

find result/ -name "*rms.info" | xargs -I {} grep --with-filename "top3 RMS RMSD" {} | sed -e "s/processed\///g" -e "s/\/3d_info\/50k.align.rms.info:\[M::__main__\] top3 RMS RMSD: /\t/g" | sort > ./stat/rmsd.info

find processed -name "*.yperx.txt" | parallel 'paste <(echo {}) <(cat {})'| sort > ./stat/yperx.stat

ls Rawdata | parallel "echo -n '{},'; samtools flagstat processed/{}/atac/{}.sort.bam | head -n 1 | sed -e 's/ /\t/g' | cut -f 1" | sort > stat/atac.read.stat
ls Rawdata | parallel "echo -n '{},'; samtools flagstat processed/{}/ct/{}.sort.bam | head -n 1 | sed -e 's/ /\t/g' | cut -f 1" | sort > stat/ct.read.stat

# Collect per-cell dedup rates for ATAC and CUT&Tag from frag.log files
for split in atac ct; do
    > stat/${split}.dedup_rate.stat
    for f in $(find processed/${split}_all -name "*.${split}.frag.log" | sort); do
        sample=$(basename "$f" ".${split}.frag.log")
        if grep -q "Combined duplication rate" "$f" 2>/dev/null; then
            rate=$(grep "Combined duplication rate" "$f" | sed 's/Combined duplication rate is \(.*\)%/\1/')
        elif grep -q "Duplication rate is" "$f" 2>/dev/null; then
            rate=$(grep "Duplication rate is" "$f" | sed 's/Duplication rate is \(.*\)%/\1/' | \
                awk 'NR==1{r1=$1} NR==2{r2=$1} END{if(NR==2) printf "%.2f", (1-(1-r1/100)*(1-r2/100))*100; else if(NR==1) printf "%.2f", r1}')
        else
            rate="NA"
        fi
        echo -e "${sample}\t${rate}" >> stat/${split}.dedup_rate.stat
    done
done

# Count RNA total and assigned reads per cell from featureCounts BAM
if [ -f processed/RNA_all/feature_gene_total/samsort.bam ]; then
    samtools view processed/RNA_all/feature_gene_total/samsort.bam | \
        awk -F'\t' '{
            split($1, a, "_");
            cell=a[2];
            total[cell]++;
            for(i=12;i<=NF;i++) {
                if($i ~ /^XS:Z:Assigned/) {
                    assigned[cell]++;
                    break
                }
            }
        } END {
            for (c in total) print c"\t"total[c]"\t"assigned[c]+0
        }' | sort > stat/rna.reads_per_cell.stat
fi

# Count RNA reads with GATC pattern (DNA contamination) per cell
for sample in $(ls Rawdata | sort); do
    if [ -f "processed/${sample}/RNA/${sample}.rna.clean.R2.fq.gz" ]; then
        zcat "processed/${sample}/RNA/${sample}.rna.clean.R2.fq.gz" | \
            awk -v sample="$sample" 'NR%4==2 {total++; if ($0 ~ /TTTTTTTT[ACG][ACGT]GATC/) gatc++} END {print sample"\t"total+0"\t"gatc+0}'
    fi
done | sort > stat/rna.dna_contam.stat
