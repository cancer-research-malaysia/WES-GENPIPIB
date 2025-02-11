#!/usr/bin/env bash
set -euo pipefail

# Function to display usage
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] MANIFEST.TSV

Options:
    -h, --help       Show this help message
    -d, --dry-run    Run in dry-run mode to test S3 file paths
    -o, --output     Specify data directory (default: script directory)
    -j, --jobs       Number of parallel jobs (default: 2)
    -r, --run-id     Specify previous run ID to rerun a pipeline (default: random UUID)

Argument:
    MANIFEST.TSV     File path to the data manifest file
EOF
    exit 0
}

# Initialize variables
DRY_RUN=false
OUTPUT_DIR=$(dirname "$0")
WORKING_DIR=$(pwd)
JOBS=2
RUN_ID=$(uuidgen | cut -d'-' -f1)
READ_DATATYPE="TUM"
S3_LOC="s3://crm.sequencing.raw.data.sharing/batch1/SLX"
S3_DEST="s3://crm.tumorstudy.analysis/suffian/WES.genotyping.outputs/WES-TUM-iter4"

# Parse command line arguments
while getopts "hdo:j:r:" opt; do
    case $opt in
        h) usage ;;
        d) DRY_RUN=true ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        r) RUN_ID="$OPTARG" ;;
        \?) usage ;;
    esac
done
shift $((OPTIND - 1))

# Check required argument
if [ $# -ne 1 ]; then
    echo "Error: Missing input manifest file." >&2
    usage
fi

MANIFEST_FILE=$1

# Directory setup function
# Initialize directory structure
init_directories() {
    local workdir=$1
    local dry_run=$2
    local run_id=$3

    if [ "$dry_run" = false ]; then
        # Create base directories
        mkdir -p "${workdir}/logs"
        # Create run-specific directories
        mkdir -p "${workdir}/flagfiles/${run_id}"

        # Check if this is a rerun by looking at log file
        if [ -f "${workdir}/logs/${run_id}-WES-pipeline.log" ]; then
            # Add clear separator for rerun
            echo -e "\n\n" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            echo "=====================================================" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            echo "RERUNNING PIPELINE WITH RUN ID: ${run_id}" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            echo "RERUN STARTED AT: $(date '+%Y-%m-%d %H:%M:%S')" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            echo "=====================================================" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            echo -e "\n" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
            log "INFO" "${workdir}" "${run_id}" "Restarting..."
        else
            # First time run - just create the file
            touch "${workdir}/logs/${run_id}-WES-pipeline.log"
        fi
    else
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create directory structure for run ${run_id}"
    fi
}

# Logging function
log() {
    local level=$1
    local workdir=$2
    local run_id=$3
    shift 3
    local color="\033[0m"
    case "$level" in
        "INFO")  color="\033[0;32m" ;;
        "WARN")  color="\033[0;33m" ;;
        "ERROR") color="\033[0;31m" ;;
    esac
    # Create the log message with timestamp
    local log_message="[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $*"
    # Output colored version to terminal
    echo -e "${color}${log_message}\033[0m" >&2
    # Output non-colored version to log file
    echo "${log_message}" >> "${workdir}/logs/${run_id}-WES-pipeline.log"
}

# Checkpoint management functions
create_checkpoint() {
    local stage=$1
    local identifier=$2
    local workdir=$3
    local dry_run=$4
    local run_id=$5

    if [ "$dry_run" = false ]; then
        touch "${workdir}/flagfiles/${run_id}/${run_id}.${identifier}.${stage}.success"
        log "INFO" "${workdir}" "${run_id}" "Created checkpoint for ${run_id}:${stage}:${identifier}"
    else
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create checkpoint for ${run_id}:${stage}:${identifier}"
    fi
}

check_checkpoint() {
    local stage=$1
    local identifier=$2
    local workdir=$3
    local dry_run=$4
    local run_id=$5

    if [ "$dry_run" = false ]; then
        [ -f "${workdir}/flagfiles/${run_id}/${run_id}.${identifier}.${stage}.success" ]
    else
        # In dry-run mode, always return false to ensure all steps would be executed
        return 1
    fi
}

mark_failure() {
    local stage=$1
    local identifier=$2
    local workdir=$3
    local dry_run=$4
    local message=$5
    local run_id=$6

    if [ "$dry_run" = false ]; then
        echo "$message" > "${workdir}/flagfiles/${run_id}/${run_id}.${identifier}.${stage}.failed"
        log "ERROR" "${workdir}" "${run_id}" "$message"
    else
        log "ERROR" "${workdir}" "${run_id}" "DRY-RUN: Would mark failure for ${run_id}:${stage}:${identifier} - $message"
    fi
    return 1
}

# helper function to initialize checkpoint directory
init_checkpoints() {
    local workdir=$1
    local dry_run=$2
    local run_id=$3

    if [ "$dry_run" = false ]; then
        log "INFO" "${workdir}" "${run_id}" "Initialized checkpoint directory"
    else
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create checkpoint directory at ${workdir}/flagfiles/${run_id}"
    fi
}

# S3 file mapping function
get_s3_files() {
    local slx_id=$1
    local tum_id=$2
    local file_prefix=$3
    local workdir=$4
    local run_id=$5
    local s3_loc=$6
    local s3_mapping=$7
    
    log "INFO" "${workdir}" "${run_id}" "Fetching file prefixes for SLX-${slx_id}..."
    aws s3 ls "${s3_loc}-${slx_id}/" | \
        grep "SLX-${file_prefix}\." | \
        awk '{print $NF}' | \
        grep 'fq.gz$' | \
        grep -v 0000 | \
        cut -d '.' -f 1,2,3,4 | \
        sort | uniq | \
        awk -v tid="$tum_id" -v slx="$slx_id" '{print slx ":" $0 ":" tid}' >> "$s3_mapping"
        # # debug lines
        # awk -v tid="$tum_id" -v slx="$slx_id" '{print slx ":" $0 ":" tid}' | \
        # tee -a "$s3_mapping"
}

# Processing functions
trim_reads() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    local dry_run=$4
    local run_id=$5
    local outdir=$6
    local workdir=$7
    local s3_loc=$8
    local s3_dest=$9
    
    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN mode enabled!"
        aws s3 cp --dryrun "${s3_loc}-${slx_id}/${prefix}.r_1.fq.gz" "${outdir}/"
        aws s3 cp --dryrun "${s3_loc}-${slx_id}/${prefix}.r_2.fq.gz" "${outdir}/"
        return 0
    fi

    if check_checkpoint "trim" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping trimming for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Trimming FASTQ files for SLX-${prefix} of sample ${tum_id}..."

    # create log directory
    mkdir -p "${workdir}/logs/1_trim_galore_out"

    # Download and process
    aws s3 cp "${s3_loc}-${slx_id}/${prefix}.r_1.fq.gz" "${outdir}/" && \
    aws s3 cp "${s3_loc}-${slx_id}/${prefix}.r_2.fq.gz" "${outdir}/" && \
    trim_galore "${outdir}/${prefix}.r_1.fq.gz" "${outdir}/${prefix}.r_2.fq.gz" \
        --paired --gzip -o "${outdir}" -a CTGTCTCTTATACACATCT \
        2>"${workdir}/logs/1_trim_galore_out/${run_id}--${prefix}_${tum_id}.trimgalore.log" && \
    aws s3 cp "${outdir}/${prefix}.r_1_val_1.fq.gz" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    aws s3 cp "${outdir}/${prefix}.r_2_val_2.fq.gz" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    aws s3 cp "${outdir}/${prefix}.r_1.fq.gz_trimming_report.txt" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    aws s3 cp "${outdir}/${prefix}.r_2.fq.gz_trimming_report.txt" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    create_checkpoint "trim" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "trim" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "TrimGalore failed for SLX-${slx_id} of sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "trim" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        rm "${outdir}/${prefix}.r_1.fq.gz" "${outdir}/${prefix}.r_2.fq.gz"
        rm "${outdir}/${prefix}.r_1.fq.gz_trimming_report.txt" "${outdir}/${prefix}.r_2.fq.gz_trimming_report.txt"
        log "INFO" "${workdir}" "${run_id}" "TrimGalore completed successfully for ${prefix}"
    fi
}

map_reads() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    local dry_run=$4
    local run_id=$5
    local outdir=$6
    local workdir=$7
    local s3_dest=$8
    
    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN enabled: Would map ${prefix}"
        return 0
    fi
    
    if ! check_checkpoint "trim" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with mapping - trimming not completed for ${prefix} of sample ${tum_id}"
        return 1
    fi

    if check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping mapping for ${prefix}_${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Mapping ${prefix} to reference genome..."

    # create log directory
    mkdir -p "${workdir}/logs/2_bwa_out"
    
    zcat "${outdir}/${prefix}.r_1_val_1.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" && \
    zcat "${outdir}/${prefix}.r_2_val_2.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" && \
    bwa mem -M -t 4 refs/GRCh38-109_bwa_db \
        "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" \
        "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" \
        2>"${workdir}/logs/2_bwa_out/${run_id}--${prefix}_${tum_id}.bwa.log" | \
    samtools view --threads 8 -b - | \
    samtools sort --threads 8 > "${outdir}/${prefix}_${tum_id}.sorted.bam" && \
    aws s3 cp "${outdir}/${prefix}_${tum_id}.sorted.bam" "${s3_dest}/${tum_id}/2_bwa_out/" && \
    create_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "BWA mapping failed for ${prefix}_${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        rm "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq"
        rm "${outdir}/${prefix}.r_1_val_1.fq.gz" "${outdir}/${prefix}.r_2_val_2.fq.gz"
        log "INFO" "${workdir}" "${run_id}" "BWA mapping completed successfully for ${prefix}_${tum_id}"
    fi
}

add_readgroups () {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    local dry_run=$4
    local run_id=$5
    local outdir=$6
    local workdir=$7
    local s3_dest=$8
    local read_datatype=$9

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would add read groups for file {$prefix} of sample ${tum_id}"
        return 0
    fi

    if ! check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with adding read groups - mapping not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "addrg" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping adding read groups for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Adding read groups for file {$prefix} of sample ${tum_id}..."

    # create log directory
    mkdir -p "${workdir}/logs/3_rg_added_bams"

    gatk AddOrReplaceReadGroups \
    -I "${outdir}/${prefix}_${tum_id}.sorted.bam" \
    -O "${outdir}/${prefix}_${tum_id}.sorted.RG-added.bam" \
    --RGID "SLX-${slx_id}" \
    --RGLB $(echo "${prefix}" | cut -d '.' -f 1,2) \
    --RGPL Illumina \
    --RGPU "SLX-${slx_id}.${tum_id}" \
    --RGSM "${tum_id}_${read_datatype}" 2>"${workdir}/logs/3_rg_added_bams/${run_id}--${prefix}_${tum_id}.RG_added.log" && \
    aws s3 cp "${outdir}/${prefix}_${tum_id}.sorted.RG-added.bam" "${s3_dest}/${tum_id}/3_rg_added_bams/" && \
    create_checkpoint "addrg" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "addrg" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "Adding read groups failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "addrg" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        rm "${outdir}/${prefix}_${tum_id}.sorted.bam"
        log "INFO" "${workdir}" "${run_id}" "Read groups added successfully for ${tum_id}"
    fi
    
}

list_bams() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6
    local s3_mapping_file="${workdir}/manifests/${run_id}--data-s3-mapping.txt"
    
    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create ${tum_id} bamlist"
        return 0
    fi

    # First check if we've already completed this step for this tumor ID
    if check_checkpoint "listbams" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping listing bams for ${tum_id} - already completed"
        return 0
    fi
    
    # Get all file prefixes for this tumor ID and verify their checkpoints
    # set ALLFILES_DONE boolean
    local ALLFILES_DONE=true
    while IFS=':' read -r slx_id prefix tumor_id; do
        if [ "${tumor_id}" = "${tum_id}" ]; then
            if ! check_checkpoint "addrg" "${prefix}_${tumor_id}" "${workdir}" "${dry_run}" "${run_id}"; then
                log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with bamlisting - mapping not completed for ${prefix}_${tumor_id}"
                ALLFILES_DONE=false
                break
            fi
        fi
    done < "${s3_mapping_file}"
    
    if [ "${ALLFILES_DONE}" = false ]; then
        return 1
    else
        log "INFO" "${workdir}" "${run_id}" "All fastq files corresponding to sample ${tum_id} have been processed."
    fi

    log "INFO" "${workdir}" "${run_id}" "Listing BAM files for sample ${tum_id}..."

    # create bamlist directory
    mkdir -p "${workdir}/logs/4_list_bams"
    
    local bam_list="${workdir}/logs/4_list_bams/${run_id}--${tum_id}_bams.list"
    
    aws s3 ls "${s3_dest}/${tum_id}/3_rg_added_bams/" | \
        grep "sorted.RG-added.bam$" | \
        awk '{print $NF}' | \
        while read -r l; do 
            echo "$l"
        done > "$bam_list" && \
    
    create_checkpoint "listbams" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "listbams" "${tum_id}" "${workdir}" "${dry_run}" "Listing bams failed for sample ${tum_id}" "${run_id}"
}

merge_bams() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    local bam_list="${workdir}/logs/4_list_bams/${run_id}--${tum_id}_bams.list"

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would merge BAMs for sample ${tum_id}"
        return 0
    fi
    
    if ! check_checkpoint "listbams" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with merging - listing bams not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping merging for ${tum_id} - already completed"
        return 0
    fi

    # check if bam list contains at least 1 line
    if [ ! -f "${bam_list}" ] || [ $(wc -l < "${bam_list}") -eq 0 ]; then
        log "ERROR" "${workdir}" "${run_id}" "No BAM files found for sample ${tum_id}"
        return 1
    fi

    log "INFO" "${workdir}" "${run_id}" "Merging BAM files for sample ${tum_id}..."

    # create log directory
    mkdir -p "${workdir}/logs/5_merged_bams"

    gatk MergeSamFiles --USE_THREADING true \
        --arguments_file "${workdir}/manifests/${tum_id}_bams.list" \
        -O "${outdir}/${tum_id}.sorted.RG-added.merged.bam" 2>"${workdir}/logs/5_merged_bams/${run-id}--${tum_id}.merge-bam.log" && \
    samtools index "${outdir}/${tum_id}.sorted.RG-added.merged.bam" "${outdir}/${tum_id}.sorted.RG-added.merged.bai" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.bam" "${s3_dest}/${tum_id}/4_merged_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.bai" "${s3_dest}/${tum_id}/4_merged_bams/" && \
    create_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "merge" "${tum_id}" "${workdir}" "${dry_run}" "Merging failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove associated RG-added bams
        rm "${outdir}/*${tum_id}.sorted.RG-added.bam"
        log "INFO" "${workdir}" "${run_id}" "Merging completed successfully for ${tum_id}"
    fi
}

# mark duplicates function
mark_dupes() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would mark duplicates for sample ${tum_id} BAM"
        return 0
    fi

    if ! check_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with marking duplicates - merging not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "markdupes" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping marking duplicates for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Marking duplicates for sample ${tum_id} BAM..."

    # create log directory
    mkdir -p "${workdir}/logs/6_dedupped_bams"

    gatk MarkDuplicates \
    -I "${outdir}/${tum_id}.sorted.RG-added.merged.bam" \
    -O "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.bam" \
    -M "${outdir}/${tum_id}.normal.MarkDup.metrics" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT 2>"${workdir}/logs/6_dedupped_bams/${run_id}--${tum_id}-markduplicates-bam.log" && \ 
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.bam" "${s3_dest}/${tum_id}/5_dedupped_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.bai" "${s3_dest}/${tum_id}/5_dedupped_bams/" && \ 
    aws s3 cp "${outdir}/${tum_id}.normal.MarkDup.metrics" "${s3_dest}/${tum_id}/5_dedupped_bams/" && \

    create_checkpoint "markdupes" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "markdupes" "${tum_id}" "${workdir}" "${dry_run}" "Marking duplicates failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "markdupes" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove associated merged bams
        rm "${outdir}/*${tum_id}.sorted.RG-added.merged.ba?"
        log "INFO" "${workdir}" "${run_id}" "Marking duplicates completed successfully for ${tum_id}"
    fi
}

# split read function
split_reads() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would split reads for sample ${tum_id} BAM"
        return 0
    fi

    if ! check_checkpoint "markdupes" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with splitting reads - marking duplicates not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "splitreads" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping splitting reads for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Splitting reads for sample ${tum_id} BAM..."

    # create directories
    mkdir -p "${workdir}/logs/7_split_reads"
    mkdir -p "${workdir}/${tum_id}_split-reads-tmp"

    gatk SplitNCigarReads \
    -I "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.bam" \
    -O "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.split.bam" \
    -R refs/GRCh38.109.fa \
    --tmp-dir "${workdir}/${tum_id}_split-reads-tmp" 2>"${workdir}/logs/7_split_reads/${run_id}--${tum_id}-splitNcigarReads.log"  && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.split.bam" "${s3_dest}/${tum_id}/6_split_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.split.bai" "${s3_dest}/${tum_id}/6_split_bams/" && \
    create_checkpoint "splitreads" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "splitreads" "${tum_id}" "${workdir}" "${dry_run}" "Splitting reads failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "splitreads" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove associated dedupped bams
        rm "${outdir}/*${tum_id}.sorted.RG-added.merged.dedup.ba?"
        rm -rf "${workdir}/${tum_id}_split-reads-tmp"
        log "INFO" "${workdir}" "${run_id}" "Splitting reads completed successfully for ${tum_id}"
    fi
}

# Base Quality Score Recalibration function - first pass (create recalibration table)
model_bqsr() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create BQSR model for sample ${tum_id} BAM"
        return 0
    fi

    if ! check_checkpoint "splitreads" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with BQSR model - splitting reads not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "modelbqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping BQSR modeling for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Creating BQSR model for sample ${tum_id} BAM..."

    # create directories
    mkdir -p "${workdir}/logs/8_model_BQSR"
    mkdir -p "${workdir}/${tum_id}_model-bqsr-tmp"

    gatk BaseRecalibrator \
    -I "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.split.bam" \
    -O "${outdir}/${tum_id}.recal_data.grp" \
    -R refs/GRCh38.109.fa \
    --tmp-dir "${workdir}/${tum_id}_model-bqsr-tmp" \
    --known-sites refs/Homo_sapiens_assembly38.known_indels.renamed.vcf \
    --known-sites refs/Homo_sapiens_assembly38.dbsnp138.renamed.vcf \
    --known-sites refs/Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf 2>"${workdir}/logs/8_model_BQSR/${run_id}--${tum_id}_model-bqsr.log" && \
    aws s3 cp "${outdir}/${tum_id}.recal_data.grp" "${s3_dest}/${tum_id}/7_bqsr_models/" && \
    create_checkpoint "modelbqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "modelbqsr" "${tum_id}" "${workdir}" "${dry_run}" "Creating recalibration model failed for sample ${tum_id}" "${run_id}"

    # clean up
    if check_checkpoint "modelbqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove tmp
        rm -rf "${workdir}/${tum_id}_model-bqsr-tmp"
        log "INFO" "${workdir}" "${run_id}" "Creating recalibration model completed successfully for ${tum_id}"
    fi
}

# Base Quality Score Recalibration function - second pass (adjust recalibration score)
apply_bqsr() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would apply BQSR for sample ${tum_id} BAM"
        return 0
    fi

    if ! check_checkpoint "modelbqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with BQSR application - model not created for ${tum_id}"
        return 1
    fi

    if check_checkpoint "applybqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping BQSR application for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Applying BQSR for sample ${tum_id} BAM..."

    # create directories
    mkdir -p "${workdir}/logs/9_apply_BQSR"
    mkdir -p "${workdir}/${tum_id}_apply-bqsr-tmp"

    gatk ApplyBQSR \
    -I "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.split.bam" \
    -O "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.recal.bam" \
    -R refs/GRCh38.109.fa \
    -bqsr "${outdir}/${tum_id}.recal_data.grp" \
    --tmp-dir "${workdir}/${tum_id}_apply-bqsr-tmp" 2>"${workdir}/logs/9_apply_BQSR/${run_id}--${tum_id}_apply-bqsr.log" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.recal.bam" "${s3_dest}/${tum_id}/8_recal_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.recal.bai" "${s3_dest}/${tum_id}/8_recal_bams/" && \

    create_checkpoint "applybqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "applybqsr" "${tum_id}" "${workdir}" "${dry_run}" "Applying recalibration failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "applybqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove associated recal bams
        rm "${outdir}/*${tum_id}.sorted.RG-added.merged.dedup.split.ba?" "${outdir}/*${tum_id}.recal_data.grp"
        rm -rf "${workdir}/${tum_id}_apply-bqsr-tmp"
        log "INFO" "${workdir}" "${run_id}" "Applying recalibration completed successfully for ${tum_id}"
    fi
}

# haplotype caller function
call_haps() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would call haplotypes for sample ${tum_id} BAM"
        return 0
    fi

    if ! check_checkpoint "applybqsr" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with haplotype calling - BQSR not applied for ${tum_id}"
        return 1
    fi

    if check_checkpoint "callhaps" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping haplotype calling for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Calling haplotypes for sample ${tum_id} BAM..."

    # create directories
    mkdir -p "${workdir}/logs/10_call_haps"
    mkdir -p "${workdir}/${tum_id}_call-haps-tmp"

    gatk HaplotypeCaller \
    -R refs/GRCh38.109.fa \
    -I "${outdir}/${tum_id}.sorted.RG-added.merged.dedup.recal.bam" \
    -O "${outdir}/${tum_id}.gatk4.all.germline.haps.vcf.gz" \
    --tmp-dir "${workdir}/${tum_id}_call-haps-tmp" \
    --dont-use-soft-clipped-bases 2>"${workdir}/logs/10_call_haps/${run_id}--${tum_id}.gatk4.all.germline.haps.vcf.gz.log" && \
    gzip -dk "${outdir}/${tum_id}.gatk4.all.germline.haps.vcf.gz" && \
    aws s3 cp "${outdir}/${tum_id}.gatk4.all.germline.haps.vcf.gz" "${s3_dest}/${tum_id}/9_call_germline_haps_vcfs/" && \
    aws s3 cp "${outdir}/${tum_id}.gatk4.all.germline.haps.vcf.gz.log" "${s3_dest}/${tum_id}/9_call_germline_haps_vcfs/" && \
    aws s3 cp "${outdir}/${tum_id}.gatk4.all.germline.haps.vcf" "${s3_dest}/${tum_id}/9_call_germline_haps_vcfs/" && \

    create_checkpoint "callhaps" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "callhaps" "${tum_id}" "${workdir}" "${dry_run}" "Haplotype calling failed for sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "callhaps" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        # remove associated recal bams
        rm "${outdir}/*${tum_id}.sorted.RG-added.merged.dedup.recal.ba?"
        rm "${outdir}/*${tum_id}.gatk4.all.germline.haps.vcf.gz" "${outdir}/*${tum_id}.gatk4.all.germline.haps.vcf"
        rm -rf "${workdir}/${tum_id}_call-haps-tmp"
        log "INFO" "${workdir}" "${run_id}" "Haplotype calling completed successfully for ${tum_id}"
    fi
}

# Main function
main() {
    local manifest=$1
    local jobs=$2
    local workdir=$3
    local outdir=$4
    local dry_run=$5
    local run_id=$6
    local s3_loc=$7
    local s3_dest=$8
    local read_datatype=$9

    local s3_mapping_file="${workdir}/manifests/${run_id}--data-s3-mapping.txt"
    local s3_tum_id_file="${workdir}/manifests/${run_id}--data-s3-tum-ids.txt"
    
    # Initialize checkpoint system
    init_checkpoints "${workdir}" "${dry_run}" "${run_id}"

    # Generate S3 mappings
    while IFS=$'\t' read -r -a line; do
        get_s3_files "${line[2]}" "${line[1]}" "${line[0]}" "${workdir}" "${run_id}" "${s3_loc}" "${s3_mapping_file}"
    done < <(tail -n +2 "${manifest}")
    
    if [ ! -s "${s3_mapping_file}" ]; then
        log "ERROR" "${workdir}" "${run_id}" "S3 data mapping file is empty. No data found at target bucket."
        exit 1
    fi
    
    log "INFO" "${workdir}" "${run_id}" "S3 prefix file created: ${s3_mapping_file}"
    awk -F':' '{print $3}' "$s3_mapping_file" | uniq > "$s3_tum_id_file"
    
    ###### MAIN PIPELINE ######
    # 1. TRIMMING
    log "INFO" "${workdir}" "${run_id}" "Starting trimming stage..."
    parallel --colsep ':' -j "$jobs" trim_reads {1} {2} {3} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_loc}" "${s3_dest}" :::: "$s3_mapping_file" && \
    log "INFO" "${workdir}" "${run_id}" "Trimming completed successfully!"

    # 2. MAPPING
    log "INFO" "${workdir}" "${run_id}" "Starting mapping stage..."
    parallel --colsep ':' -j "$jobs" map_reads {1} {2} {3} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_mapping_file" && \
    log "INFO" "${workdir}" "${run_id}" "Mapping completed successfully!"

    # 3. ADDING READGROUPS
    log "INFO" "${workdir}" "${run_id}" "Starting readgroup addition stage..."
    parallel --colsep ':' -j "$jobs" add_readgroups {1} {2} {3} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" "${read_datatype}" :::: "$s3_mapping_file" && \
    log "INFO" "${workdir}" "${run_id}" "Readgroup addition completed successfully!"

    # 4. LISTING BAMS
    log "INFO" "${workdir}" "${run_id}" "Listing bam-level tumor IDs..."
    parallel --colsep ':' -j "$jobs" list_bams {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "Listing bams completed successfully!"

    # 5. MERGING BAMS
    log "INFO" "${workdir}" "${run_id}" "Starting merging stage..."
    parallel --colsep ':' -j "$jobs" merge_bams {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "Merging completed successfully!"

    # 6. MARKING DUPLICATES
    log "INFO" "${workdir}" "${run_id}" "Starting marking duplicates stage..."
    parallel --colsep ':' -j "$jobs" mark_dupes {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "Marking duplicates completed successfully!"

    # 7. SPLITTING READS
    log "INFO" "${workdir}" "${run_id}" "Starting splitting reads stage..."
    parallel --colsep ':' -j "$jobs" split_reads {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "Splitting reads completed successfully!"

    # 8. BASE QUALITY SCORE RECALIBRATION (first pass) - model creation
    log "INFO" "${workdir}" "${run_id}" "Starting BQSR model creation stage..."
    parallel --colsep ':' -j "$jobs" model_bqsr {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "BQSR model creation completed successfully!"

    # 9. BASE QUALITY SCORE RECALIBRATION (second pass) - apply recalibration
    log "INFO" "${workdir}" "${run_id}" "Starting BQSR application stage..."
    parallel --colsep ':' -j "$jobs" apply_bqsr {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "BQSR application completed successfully!"

    # 10. HAPLOTYPE CALLING
    log "INFO" "${workdir}" "${run_id}" "Starting haplotype calling stage..."
    parallel --colsep ':' -j "$jobs" call_haps {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file" && \
    log "INFO" "${workdir}" "${run_id}" "Haplotype calling completed successfully!"


    log "INFO" "${workdir}" "${run_id}" "Pipeline has finished without errors at $(date +"%Y-%m-%d %H:%M:%S"). Please check logs for more details."
}

# Export functions for parallel execution
export -f \
        init_directories \
        log \
        create_checkpoint \
        check_checkpoint \
        mark_failure \
        init_checkpoints \
            trim_reads \
            map_reads \
            add_readgroups \
            list_bams \
            merge_bams \
            mark_dupes \
            split_reads \
            model_bqsr \
            apply_bqsr \
            call_haps

# Initialize directory structure before any logging occurs
init_directories "${WORKING_DIR}" "${DRY_RUN}" "${RUN_ID}"
log "INFO" "${WORKING_DIR}" "${RUN_ID}" "Initialized directory structure."
log "INFO" "${WORKING_DIR}" "${RUN_ID}" "Current RUN ID is: ${RUN_ID}"

# check for flagfile run-id directory to see if it was a previous run
if [ -d "${WORKING_DIR}/flagfiles/${RUN_ID}" ]; then
    log "INFO" "${WORKING_DIR}" "${RUN_ID}" "Found previous run with same RUN ID. This will be treated as a re-run/resume."
    # empty mapping file to avoid redundant appending
    : > "${WORKING_DIR}/manifests/${RUN_ID}--data-s3-mapping.txt"

fi

# run the main function
main "$MANIFEST_FILE" "$JOBS" "$WORKING_DIR" "$OUTPUT_DIR" "$DRY_RUN" "$RUN_ID" "$S3_LOC" "$S3_DEST" "$READ_DATATYPE"

