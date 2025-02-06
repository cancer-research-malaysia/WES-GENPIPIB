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
    local workdir=$3
    local dry_run=$4

    if [ "$dry_run" = false ]; then
        touch "${workdir}/flagfiles/${identifier}.${stage}.success"
        log "INFO" "Created checkpoint for ${stage}:${identifier}"
    else
        log "INFO" "DRY-RUN: Would create checkpoint for ${stage}:${identifier}"
    fi
}

check_checkpoint() {
    local stage=$1
    local identifier=$2
    local workdir=$3
    local dry_run=$4

    if [ "$dry_run" = false ]; then
        [ -f "${workdir}/flagfiles/${identifier}.${stage}.success" ]
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

    if [ "$dry_run" = false ]; then
        echo "$message" > "${workdir}/flagfiles/${identifier}.${stage}.failed"
        log "ERROR" "$message"
    else
        log "ERROR" "DRY-RUN: Would mark failure for ${stage}:${identifier} - $message"
    fi
    return 1
}

# helper function to initialize checkpoint directory
init_checkpoints() {
    local workdir=$1
    local dry_run=$2

    if [ "$dry_run" = false ]; then
        mkdir -p "${workdir}/flagfiles"
        log "INFO" "Initialized checkpoint directory"
    else
        log "INFO" "DRY-RUN: Would create checkpoint directory at ${workdir}/flagfiles"
    fi
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
        awk -v tid="$tum_id" -v slx="$slx_id" '{print slx ":" $0 ":" tid}' >> "$s3_mapping"
        # # debug lines
        # awk -v tid="$tum_id" -v slx="$slx_id" '{print slx ":" $0 ":" tid}' | \
        # tee -a "$s3_mapping"
}

# Processing functions
process_trimming() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    local dry_run=$4
    local outdir=$5
    local workdir=$6
    local s3_loc=$7
    local s3_dest=$8
    
    if check_checkpoint "trim" "${prefix}" "${workdir}" "$dry_run"; then
        log "INFO" "Skipping trimming for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "Trimming FASTQ files for SLX-${slx_id} of sample ${tum_id}..."

    if [ "$dry_run" = true ]; then
        log "INFO" "DRY-RUN mode enabled!"
        aws s3 cp --dryrun "${s3_loc}-${slx_id}/${prefix}.r_1.fq.gz" "${outdir}/"
        aws s3 cp --dryrun "${s3_loc}-${slx_id}/${prefix}.r_2.fq.gz" "${outdir}/"
        return 0
    fi

    # Download and process
    aws s3 cp "${s3_loc}-${slx_id}/${prefix}.r_1.fq.gz" "${outdir}/" && \
    aws s3 cp "${s3_loc}-${slx_id}/${prefix}.r_2.fq.gz" "${outdir}/" && \
    trim_galore "${outdir}/${prefix}.r_1.fq.gz" "${outdir}/${prefix}.r_2.fq.gz" \
        --paired --gzip -o "${outdir}" -a CTGTCTCTTATACACATCT \
        2>"${workdir}/logs/${prefix}.trimgalore.log" && \
    aws s3 cp "${outdir}/${prefix}.r_1_val_1.fq.gz" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    aws s3 cp "${outdir}/${prefix}.r_2_val_2.fq.gz" "${s3_dest}/${tum_id}/1_trim_galore_out/" && \
    create_checkpoint "trim" "${prefix}" "${workdir}" "$dry_run" || \
    mark_failure "trim" "${prefix}" "${workdir}" "$dry_run" "TrimGalore failed for SLX-${slx_id} of sample ${tum_id}"

    # Cleanup on success
    if check_checkpoint "trim" "${prefix}" "${workdir}" "$dry_run"; then
        rm "${outdir}/${prefix}.r_1.fq.gz" "${outdir}/${prefix}.r_2.fq.gz"
        log "INFO" "TrimGalore completed successfully for ${prefix}"
    fi
}

process_mapping() {
    local slx_id=$1
    local prefix=$2
    local tum_id=$3
    local dry_run=$4
    local outdir=$5
    local workdir=$6
    local s3_dest=$7
    
    if [ "$dry_run" = true ]; then
        log "INFO" "DRY-RUN enabled: Would map ${prefix}"
        return 0
    fi
    
    if ! check_checkpoint "trim" "${prefix}" "${workdir}" "$dry_run"; then
        log "ERROR" "Cannot proceed with mapping - trimming not completed for ${prefix}"
        return 1
    fi

    if check_checkpoint "map" "${prefix}" "${workdir}" "$dry_run"; then
        log "INFO" "Skipping mapping for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "Mapping ${prefix} to reference genome..."
    
    zcat "${outdir}/${prefix}.r_1_val_1.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" && \
    zcat "${outdir}/${prefix}.r_2_val_2.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" && \
    bwa mem -M -t 4 "${workdir}/refs/GRCh38.109.fa" \
        "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" \
        "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" \
        2>"${workdir}/logs/${prefix}.bwa.log" | \
    samtools view --threads 8 -b - | \
    samtools sort --threads 8 > "${outdir}/${prefix}_${tum_id}.sorted.bam" && \
    aws s3 cp "${outdir}/${prefix}_${tum_id}.sorted.bam" "${s3_dest}/${tum_id}/2_bwa_out/" && \
    create_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "$dry_run" || \
    mark_failure "map" "${prefix}_${tum_id}" "${workdir}" "$dry_run" "BWA mapping failed for ${prefix}_${tum_id}"

    # Cleanup on success
    if check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "$dry_run"; then
        rm "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq"
        log "INFO" "BWA mapping completed successfully for ${prefix}_${tum_id}"
    fi
}

process_bamlisting() {
    local tum_id=$1
    local dry_run=$2
    local outdir=$3
    local workdir=$4
    local s3_dest=$5
    
    if [ "$dry_run" = true ]; then
        log "INFO" "DRY-RUN: Would create ${tum_id} bamlist"
        return 0
    fi
    
    # if ! check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "$dry_run"; then
    #     log "ERROR" "Cannot proceed with bamlisting - mapping not completed for ${prefix}_${tum_id}"
    #     return 1
    # fi

    if check_checkpoint "bamlist" "${tum_id}" "${workdir}" "$dry_run"; then
        log "INFO" "Skipping bamlisting for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "Listing BAM files for sample ${tum_id}..."
    
    local bam_list="${workdir}/manifests/${tum_id}_bams.list"
    
    aws s3 ls "${s3_dest}/${tum_id}/2_bwa_out/" | \
        grep ".sorted.bam$" | \
        awk '{print $NF}' | \
        while read -r l; do 
            echo "-I $l"
        done >> "$bam_list"
    
    create_checkpoint "bamlist" "${tum_id}" "${workdir}" "$dry_run" || \
    mark_failure "bamlist" "${tum_id}" "${workdir}" "$dry_run" "Bamlisting failed for sample ${tum_id}"
}

process_merging() {
    local tum_id=$1
    local dry_run=$2
    local outdir=$3
    local workdir=$4
    local s3_dest=$5

    if [ "$dry_run" = true ]; then
        log "INFO" "DRY-RUN: Would merge BAMs for sample ${tum_id}"
        return 0
    fi
    
    if ! check_checkpoint "bamlist" "${tum_id}" "${workdir}" "$dry_run"; then
        log "ERROR" "Cannot proceed with merging - bamlisting not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "merge" "${tum_id}" "${workdir}" "$dry_run"; then
        log "INFO" "Skipping merging for ${tum_id} - already completed"
        return 0
    fi

    # check if bam list contains at least 1 line
    if [ ! -f "${workdir}/manifests/${tum_id}_bams.list" ] || [ $(wc -l < "${workdir}/manifests/${tum_id}_bams.list") -eq 0 ]; then
        log "ERROR" "No BAM files found for sample ${tum_id}"
        return 1
    fi

    log "INFO" "Merging BAM files for sample ${tum_id}..."

    gatk MergeSamFiles --USE_THREADING true \
        --arguments_file "${workdir}/manifests/${tum_id}_bams.list" \
        -O "${outdir}/${tum_id}.merged.bam" && \
    samtools index "${outdir}/${tum_id}.merged.bam" "${outdir}/${tum_id}.merged.bai" && \
    aws s3 cp "${outdir}/${tum_id}.merged.bam" "${s3_dest}/${tum_id}/3_merged_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.merged.bai" "${s3_dest}/${tum_id}/3_merged_bams/" && \
    create_checkpoint "merge" "${tum_id}" "${workdir}" "$dry_run" || \
    mark_failure "merge" "${tum_id}" "${workdir}" "$dry_run" "Merging failed for sample ${tum_id}"
}

main() {
    local manifest=$1
    local jobs=$2
    local workdir=$3
    local outdir=$4
    local dry_run=$5
    local s3_loc=$6
    local s3_dest=$7

    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local s3_mapping_file="${workdir}/manifests/data-s3-mapping--${timestamp}.txt"
    local s3_tum_id_file="${workdir}/manifests/data-s3-tum-ids--${timestamp}.txt"
    
    # Create necessary directories
    if [ "$dry_run" = false ]; then
        mkdir -p "${workdir}"/{logs,flagfiles}
    else
        log "INFO" "DRY-RUN: Would create directories: ${workdir}/{logs,flagfiles}"
    fi
    
    # Initialize checkpoint system
    init_checkpoints "${workdir}" "${dry_run}"

    # Generate S3 mappings
    while IFS=$'\t' read -r -a line; do
        get_s3_files "${line[2]}" "${line[1]}" "${line[0]}" "${s3_loc}" "${s3_mapping_file}"
    done < <(tail -n +2 "${manifest}")
    
    if [ ! -s "${s3_mapping_file}" ]; then
        log "ERROR" "S3 data mapping file is empty. No data found at target bucket."
        exit 1
    fi
    
    log "INFO" "S3 prefix file created: ${s3_mapping_file}"
    awk -F':' '{print $3}' "$s3_mapping_file" | uniq > "$s3_tum_id_file"
    
    # Process each stage
    log "INFO" "Starting trimming stage..."
    parallel --colsep ':' -j "$jobs" process_trimming {1} {2} {3} "${dry_run}" "${outdir}" "${workdir}" "${s3_loc}" "${s3_dest}" :::: "$s3_mapping_file"
    
    log "INFO" "Starting mapping stage..."
    parallel --colsep ':' -j "$jobs" process_mapping {1} {2} {3} "${dry_run}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_mapping_file"
    
    log "INFO" "Starting bamlisting stage..."
    parallel --colsep ':' -j "$jobs" process_bamlisting {1} "${dry_run}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file"

    log "INFO" "Starting merging stage..."
    parallel --colsep ':' -j "$jobs" process_merging {1} "${dry_run}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file"
    
    log "INFO" "Pipeline completed successfully!"
}

# Export functions for parallel execution
export -f log create_checkpoint check_checkpoint mark_failure init_checkpoints \
    process_trimming process_mapping process_bamlisting process_merging

main "$MANIFEST_FILE" "$JOBS" "$WORKING_DIR" "$OUTPUT_DIR" "$DRY_RUN" "$S3_LOC" "$S3_DEST"

