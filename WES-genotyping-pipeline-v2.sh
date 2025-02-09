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
JOBS=3

S3_LOC="s3://crm.sequencing.raw.data.sharing/batch1/SLX"
S3_DEST="s3://crm.tumorstudy.analysis/suffian/WES.genotyping.outputs/WES-TUM"

# generate a ID for the run
RUN_ID=$(head -c 8 /dev/urandom | xxd -p)

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
    echo "${log_message}" >> "${workdir}/logs/${run_id}-pipeline-$(date '+%Y-%m-%d').log"
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
        mkdir -p "${workdir}/flagfiles/${run_id}"
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
    local run_id=$5
    local outdir=$6
    local workdir=$7
    local s3_loc=$8
    local s3_dest=$9
    
    if check_checkpoint "trim" "${prefix}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping trimming for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Trimming FASTQ files for SLX-${slx_id} of sample ${tum_id}..."

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN mode enabled!"
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
    create_checkpoint "trim" "${prefix}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "trim" "${prefix}" "${workdir}" "${dry_run}" "TrimGalore failed for SLX-${slx_id} of sample ${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "trim" "${prefix}" "${workdir}" "${dry_run}" "${run_id}"; then
        rm "${outdir}/${prefix}.r_1.fq.gz" "${outdir}/${prefix}.r_2.fq.gz"
        log "INFO" "${workdir}" "${run_id}" "TrimGalore completed successfully for ${prefix}"
    fi
}

process_mapping() {
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
    
    if ! check_checkpoint "trim" "${prefix}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with mapping - trimming not completed for ${prefix}"
        return 1
    fi

    if check_checkpoint "map" "${prefix}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping mapping for ${prefix} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Mapping ${prefix} to reference genome..."
    
    zcat "${outdir}/${prefix}.r_1_val_1.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" && \
    zcat "${outdir}/${prefix}.r_2_val_2.fq.gz" | \
        awk '(NR%2==0){$0=substr($0,1,75)}{print}' > "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" && \
    bwa mem -M -t 4 refs/GRCh38-109_bwa_db \
        "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" \
        "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq" \
        2>"${workdir}/logs/${prefix}.bwa.log" | \
    samtools view --threads 8 -b - | \
    samtools sort --threads 8 > "${outdir}/${prefix}_${tum_id}.sorted.bam" && \
    aws s3 cp "${outdir}/${prefix}_${tum_id}.sorted.bam" "${s3_dest}/${tum_id}/2_bwa_out/" && \
    create_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "BWA mapping failed for ${prefix}_${tum_id}" "${run_id}"

    # Cleanup on success
    if check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        rm "${outdir}/${prefix}_${tum_id}.r_1_bwa_in.fq" "${outdir}/${prefix}_${tum_id}.r_2_bwa_in.fq"
        log "INFO" "${workdir}" "${run_id}" "BWA mapping completed successfully for ${prefix}_${tum_id}"
    fi
}

process_bamlisting() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6
    
    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create ${tum_id} bamlist"
        return 0
    fi
    
    # if ! check_checkpoint "map" "${prefix}_${tum_id}" "${workdir}" "$dry_run"; then
    #     log "ERROR" "Cannot proceed with bamlisting - mapping not completed for ${prefix}_${tum_id}"
    #     return 1
    # fi

    if check_checkpoint "bamlist" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping bamlisting for ${tum_id} - already completed"
        return 0
    fi

    log "INFO" "${workdir}" "${run_id}" "Listing BAM files for sample ${tum_id}..."
    
    local bam_list="${workdir}/manifests/${tum_id}_bams.list"
    
    aws s3 ls "${s3_dest}/${tum_id}/2_bwa_out/" | \
        grep ".sorted.bam$" | \
        awk '{print $NF}' | \
        while read -r l; do 
            echo "-I $l"
        done >> "$bam_list"
    
    create_checkpoint "bamlist" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "bamlist" "${tum_id}" "${workdir}" "${dry_run}" "Bamlisting failed for sample ${tum_id}" "${run_id}"
}

process_merging() {
    local tum_id=$1
    local dry_run=$2
    local run_id=$3
    local outdir=$4
    local workdir=$5
    local s3_dest=$6

    if [ "$dry_run" = true ]; then
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would merge BAMs for sample ${tum_id}"
        return 0
    fi
    
    if ! check_checkpoint "bamlist" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "ERROR" "${workdir}" "${run_id}" "Cannot proceed with merging - bamlisting not completed for ${tum_id}"
        return 1
    fi

    if check_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}"; then
        log "INFO" "${workdir}" "${run_id}" "Skipping merging for ${tum_id} - already completed"
        return 0
    fi

    # check if bam list contains at least 1 line
    if [ ! -f "${workdir}/manifests/${tum_id}_bams.list" ] || [ $(wc -l < "${workdir}/manifests/${tum_id}_bams.list") -eq 0 ]; then
        log "ERROR" "${workdir}" "${run_id}" "No BAM files found for sample ${tum_id}"
        return 1
    fi

    log "INFO" "${workdir}" "${run_id}" "Merging BAM files for sample ${tum_id}..."

    gatk MergeSamFiles --USE_THREADING true \
        --arguments_file "${workdir}/manifests/${tum_id}_bams.list" \
        -O "${outdir}/${tum_id}.merged.bam" && \
    samtools index "${outdir}/${tum_id}.merged.bam" "${outdir}/${tum_id}.merged.bai" && \
    aws s3 cp "${outdir}/${tum_id}.merged.bam" "${s3_dest}/${tum_id}/3_merged_bams/" && \
    aws s3 cp "${outdir}/${tum_id}.merged.bai" "${s3_dest}/${tum_id}/3_merged_bams/" && \
    create_checkpoint "merge" "${tum_id}" "${workdir}" "${dry_run}" "${run_id}" || \
    mark_failure "merge" "${tum_id}" "${workdir}" "${dry_run}" "Merging failed for sample ${tum_id}" "${run_id}"
}

main() {
    local manifest=$1
    local jobs=$2
    local workdir=$3
    local outdir=$4
    local dry_run=$5
    local run_id=$6
    local s3_loc=$7
    local s3_dest=$8

    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local s3_mapping_file="${workdir}/manifests/data-s3-mapping--${timestamp}.txt"
    local s3_tum_id_file="${workdir}/manifests/data-s3-tum-ids--${timestamp}.txt"

    # print the run id
    log "INFO" "${workdir}" "${run_id}" "The current run ID: ${run_id}"
    
    # Create necessary directories
    if [ "$dry_run" = false ]; then
        mkdir -p "${workdir}"/{logs,flagfiles}
        log "INFO" "${workdir}" "${run_id}" "Created directories: ${workdir}/{logs,flagfiles}"
    else
        log "INFO" "${workdir}" "${run_id}" "DRY-RUN: Would create directories: ${workdir}/{logs,flagfiles}"
    fi
    
    # Initialize checkpoint system
    init_checkpoints "${workdir}" "${dry_run}" "${run_id}"

    # Generate S3 mappings
    while IFS=$'\t' read -r -a line; do
        get_s3_files "${line[2]}" "${line[1]}" "${line[0]}" "${s3_loc}" "${s3_mapping_file}"
    done < <(tail -n +2 "${manifest}")
    
    if [ ! -s "${s3_mapping_file}" ]; then
        log "ERROR" "${workdir}" "${run_id}" "S3 data mapping file is empty. No data found at target bucket."
        exit 1
    fi
    
    log "INFO" "${workdir}" "${run_id}" "S3 prefix file created: ${s3_mapping_file}"
    awk -F':' '{print $3}' "$s3_mapping_file" | uniq > "$s3_tum_id_file"
    
    # Process each stage
    log "INFO" "${workdir}" "${run_id}" "Starting trimming stage..."
    parallel --colsep ':' -j "$jobs" process_trimming {1} {2} {3} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_loc}" "${s3_dest}" :::: "$s3_mapping_file"
    
    log "INFO" "${workdir}" "${run_id}" "Starting mapping stage..."
    parallel --colsep ':' -j "$jobs" process_mapping {1} {2} {3} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_mapping_file"
    
    log "INFO" "${workdir}" "${run_id}" "Starting bamlisting stage..."
    parallel --colsep ':' -j "$jobs" process_bamlisting {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file"

    log "INFO" "${workdir}" "${run_id}" "Starting merging stage..."
    parallel --colsep ':' -j "$jobs" process_merging {1} "${dry_run}" "${run_id}" "${outdir}" "${workdir}" "${s3_dest}" :::: "$s3_tum_id_file"
    
    log "INFO" "${workdir}" "${run_id}" "Pipeline completed successfully!"
}

# Export functions for parallel execution
export -f log create_checkpoint check_checkpoint mark_failure init_checkpoints \
    process_trimming process_mapping process_bamlisting process_merging

main "$MANIFEST_FILE" "$JOBS" "$WORKING_DIR" "$OUTPUT_DIR" "$DRY_RUN" "$RUN_ID" "$S3_LOC" "$S3_DEST"

