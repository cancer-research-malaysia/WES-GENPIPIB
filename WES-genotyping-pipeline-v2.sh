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
    -j, --jobs       Number of parallel jobs (default: 4)

Argument:
    MANIFEST.TSV     File path to the data manifest file
EOF
    exit 0
}

# Initialize variables
DRY_RUN=false
OUTPUT_DIR=$(dirname "$0")
WORKING_DIR=$(pwd)
JOBS=4

S3_LOC="s3://crm.sequencing.raw.data.sharing/batch1/SLX"
S3_DEST="s3://crm.tumorstudy.analysis/suffian/WES.genotyping.outputs/WES-TUM"

# Parse command line arguments
while getopts "hdo:j:" opt; do
    case $opt in
        h) usage ;;
        d) DRY_RUN=true ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
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

# Logging function
log() {
    local level=$1
    shift
    local color="\033[0m"
    case "$level" in
        "INFO")  color="\033[0;32m" ;;
        "WARN")  color="\033[0;33m" ;;
        "ERROR") color="\033[0;31m" ;;
    esac
    echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $*\033[0m" >&2
}

# Checkpoint management functions
create_checkpoint() {
    local stage=$1
    local identifier=$2
    touch "${WORKING_DIR}/flagfiles/${identifier}.${stage}.success"
}

check_checkpoint() {
    local stage=$1
    local identifier=$2
    [ -f "${WORKING_DIR}/flagfiles/${identifier}.${stage}.success" ]
}

mark_failure() {
    local stage=$1
    local identifier=$2
    local message=$3
    echo "$message" > "${WORKING_DIR}/flagfiles/${identifier}.${stage}.failed"
    log "ERROR" "$message"
    return 1
}

# S3 file mapping function
get_s3_files() {
    local slx_id=$1
    local tum_id=$2
    local file_prefix=$3
    local s3_loc=$4
    local s3_mapping=$5
    
    log "INFO" "Fetching file prefixes for SLX-${slx_id}..."
    aws s3 ls "${s3_loc}-${slx_id}/" | \
        grep "SLX-${file_prefix}\." | \
        awk '{print $NF}' | \
        grep 'fq.gz$' | \
        grep -v 0000 | \
        cut -d '.' -f 1,2,3,4 | \
        sort | uniq | \
        awk -v tid="$tum_id" -v slx="$slx_id" '{print slx ":" $0 ":" tid}' | \
        tee -a "$s3_mapping"
}

# Processing functions
process_trimming() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    
    if check_checkpoint "trim" "${prefix}"; then
        log "INFO" "Skipping trimming for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "Trimming FASTQ files for SLX-${slx_id} of sample ${tum_id}..."
    
    if [ "$DRY_RUN" = true ]; then
        aws s3 cp --dryrun "${S3_LOC}-${slx_id}/${prefix}.r_1.fq.gz" "${OUTPUT_DIR}/"
        aws s3 cp --dryrun "${S3_LOC}-${slx_id}/${prefix}.r_2.fq.gz" "${OUTPUT_DIR}/"
        return 0
    fi

    # Download and process
    aws s3 cp "${S3_LOC}-${slx_id}/${prefix}.r_1.fq.gz" "${OUTPUT_DIR}/" && \
    aws s3 cp "${S3_LOC}-${slx_id}/${prefix}.r_2.fq.gz" "${OUTPUT_DIR}/" && \
    trim_galore "${OUTPUT_DIR}/${prefix}.r_1.fq.gz" "${OUTPUT_DIR}/${prefix}.r_2.fq.gz" \
        --paired --gzip -o "${OUTPUT_DIR}" -a CTGTCTCTTATACACATCT \
        2>"${WORKING_DIR}/logs/${prefix}.trimgalore.log" && \
    aws s3 cp "${OUTPUT_DIR}/${prefix}.r_1_val_1.fq.gz" "${S3_DEST}/${tum_id}/1_trim_galore_out/" && \
    aws s3 cp "${OUTPUT_DIR}/${prefix}.r_2_val_2.fq.gz" "${S3_DEST}/${tum_id}/1_trim_galore_out/" && \
    create_checkpoint "trim" "${prefix}" || \
    mark_failure "trim" "${prefix}" "TrimGalore failed for SLX-${slx_id} of sample ${tum_id}"

    # Cleanup on success
    if check_checkpoint "trim" "${prefix}"; then
        rm "${OUTPUT_DIR}/${prefix}.r_1.fq.gz" "${OUTPUT_DIR}/${prefix}.r_2.fq.gz"
        log "INFO" "TrimGalore completed successfully for ${prefix}"
    fi
}

process_mapping() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    
    if ! check_checkpoint "trim" "${prefix}"; then
        log "ERROR" "Cannot proceed with mapping - trimming not completed for ${prefix}"
        return 1
    fi

    if check_checkpoint "map" "${prefix}"; then
        log "INFO" "Skipping mapping for ${prefix} - already completed"
        return 0
    fi

    if [ "$DRY_RUN" = true ]; then
        log "INFO" "DRY-RUN: Would map ${prefix}"
        return 0
    fi

    log "INFO" "Mapping ${prefix} to reference genome..."
    
    zcat "${OUTPUT_DIR}/${prefix}.r_1_val_1.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${OUTPUT_DIR}/${prefix}_${tum_id}.r_1_bwa_in.fq" && \
    zcat "${OUTPUT_DIR}/${prefix}.r_2_val_2.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${OUTPUT_DIR}/${prefix}_${tum_id}.r_2_bwa_in.fq" && \
    bwa mem -M -t 4 "${WORKING_DIR}/refs/GRCh38.109.fa" \
        "${OUTPUT_DIR}/${prefix}_${tum_id}.r_1_bwa_in.fq" \
        "${OUTPUT_DIR}/${prefix}_${tum_id}.r_2_bwa_in.fq" \
        2>"${WORKING_DIR}/logs/${prefix}.bwa.log" | \
    samtools view --threads 8 -b - | \
    samtools sort --threads 8 > "${OUTPUT_DIR}/${prefix}_${tum_id}.sorted.bam" && \
    aws s3 cp "${OUTPUT_DIR}/${prefix}_${tum_id}.sorted.bam" "${S3_DEST}/${tum_id}/2_bwa_out/" && \
    create_checkpoint "map" "${prefix}" || \
    mark_failure "map" "${prefix}" "BWA mapping failed for ${prefix}"

    # Cleanup on success
    if check_checkpoint "map" "${prefix}"; then
        rm "${OUTPUT_DIR}/${prefix}_${tum_id}.r_1_bwa_in.fq" "${OUTPUT_DIR}/${prefix}_${tum_id}.r_2_bwa_in.fq"
        log "INFO" "BWA mapping completed successfully for ${prefix}"
    fi
}

process_merging() {
    local slx_id=$1
    local tum_id=$2
    
    if check_checkpoint "merge" "${tum_id}"; then
        log "INFO" "Skipping merging for ${tum_id} - already completed"
        return 0
    fi

    if [ "$DRY_RUN" = true ]; then
        log "INFO" "DRY-RUN: Would merge BAMs for ${tum_id}"
        return 0
    fi

    log "INFO" "Merging BAM files for sample ${tum_id}..."
    
    local bam_list="${WORKING_DIR}/manifests/${tum_id}_bams.list"
    aws s3 ls "${S3_DEST}/${tum_id}/2_bwa_out/" | \
        grep "SLX-${slx_id}." | \
        grep ".sorted.bam$" | \
        awk '{print $NF}' | \
        while read -r l; do 
            echo "-I $l"
        done >> "$bam_list"

    gatk MergeSamFiles --USE_THREADING true \
        --arguments_file "$bam_list" \
        -O "${OUTPUT_DIR}/${tum_id}.merged.bam" && \
    samtools index "${OUTPUT_DIR}/${tum_id}.merged.bam" "${OUTPUT_DIR}/${tum_id}.merged.bai" && \
    aws s3 cp "${OUTPUT_DIR}/${tum_id}.merged.bam" "${S3_DEST}/${tum_id}/3_merged_bams/" && \
    aws s3 cp "${OUTPUT_DIR}/${tum_id}.merged.bai" "${S3_DEST}/${tum_id}/3_merged_bams/" && \
    create_checkpoint "merge" "${tum_id}" || \
    mark_failure "merge" "${tum_id}" "Merging failed for sample ${tum_id}"
}

main() {
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local s3_mapping_file="${WORKING_DIR}/manifests/data-s3-mapping--${timestamp}.txt"
    local s3_idpairs_file="${WORKING_DIR}/manifests/data-s3-IDpairs--${timestamp}.txt"
    
    # Create necessary directories
    mkdir -p "${WORKING_DIR}"/{logs,flagfiles}
    
    # Generate S3 mappings
    while IFS=$'\t' read -r -a line; do
        get_s3_files "${line[2]}" "${line[1]}" "${line[0]}" "$S3_LOC" "$s3_mapping_file"
    done < <(tail -n +2 "${MANIFEST_FILE}")
    
    if [ ! -s "${s3_mapping_file}" ]; then
        log "ERROR" "S3 data mapping file is empty. No data found at target bucket."
        exit 1
    fi
    
    log "INFO" "S3 prefix file created: ${s3_mapping_file}"
    awk -F':' '{print $1":"$3}' "$s3_mapping_file" | uniq > "$s3_idpairs_file"
    
    # Process each stage
    log "INFO" "Starting trimming stage..."
    parallel --colsep ':' -j "$JOBS" process_trimming {1} {2} {3} :::: "$s3_mapping_file"
    
    log "INFO" "Starting mapping stage..."
    parallel --colsep ':' -j "$JOBS" process_mapping {1} {2} {3} :::: "$s3_mapping_file"
    
    log "INFO" "Starting merging stage..."
    parallel --colsep ':' -j "$JOBS" process_merging {1} {2} :::: "$s3_idpairs_file"
    
    log "INFO" "Pipeline completed successfully"
}

# Export functions for parallel execution
export -f log create_checkpoint check_checkpoint mark_failure \
    process_trimming process_mapping process_merging

main
