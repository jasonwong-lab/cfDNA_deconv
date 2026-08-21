import os
import subprocess
import logging
from pathlib import Path
from cfDNA_deconv.utils import marker_bed, resolve_marker_dir, run_command, validate_inputs

logger = logging.getLogger(__name__)

def deconvolute_short_read(
    marker_dir: str,
    ref_dir: str,
    bam_file: str,
    output_dir: str,
    mode: str = "PE",
    aligner: str = "bismark",
    cores: int = 4,
    genome: str = "hg19",
    marker_set: str = "U250_36celltype"
):
    """
    Deconvolute cfDNA composition from short-read sequencing data (BAM file)
    based on methylation markers.
    
    Args:
        marker_dir: Path to folder with .bed markers and ReferenceList.txt
        ref_dir: Path to folder with the selected genome FASTA
        bam_file: Path to sorted/indexed BAM file
        output_dir: Output directory for results
        mode: Sequencing mode (PE/SE), default=PE
        aligner: Aligner used (bsmap/bismark/segemehl/gem), default=bismark
        cores: Number of threads, default=4
        genome: Reference genome build (hg19/hg38)
        marker_set: Marker set name (for example, U250_36celltype)
    """
    # Validate inputs
    validate_inputs([marker_dir, ref_dir, bam_file], ["dir", "dir", "file"])
    marker_dir = resolve_marker_dir(marker_dir, genome, marker_set)
    ref_dir = Path(ref_dir)
    bam_file = Path(bam_file)
    output_dir = Path(output_dir)
    
    # Preprocessing
    sample_name = bam_file.stem
    logger.info(f"Starting short-read deconvolution for sample: {sample_name}")
    
    # Create output directories
    output_subdirs = [
        output_dir / "Region_Bam",
        output_dir / "Region_Bed",
        output_dir / "ReadsRatio",
        output_dir / "RPKM"
    ]
    for subdir in output_subdirs:
        subdir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Extract tissue-specific reads (samtools view)
    reference_list = marker_dir / "ReferenceList.txt"
    if not reference_list.exists():
        raise FileNotFoundError(f"ReferenceList.txt not found in {marker_dir}")
    
    with open(reference_list, "r") as f:
        for line in f:
            cell_type = line.strip()
            if not cell_type:
                continue
            bed_file = marker_bed(marker_dir, cell_type, genome, marker_set)
            if not bed_file.exists():
                logger.warning(f"Marker bed file missing: {bed_file}")
                continue
            
            output_bam = output_subdirs[0] / f"{sample_name}_{cell_type}.bam"
            cmd = [
                "samtools", "view", str(bam_file),
                "-M", "-L", str(bed_file),
                "-b", "-o", str(output_bam),
                "-@", str(cores)
            ]
            run_command(cmd, f"Extracting reads for {cell_type}")
    
    # Step 2: Calculate read-level methylation (RLM)
    region_bam_dir = output_subdirs[0]
    region_bed_dir = output_subdirs[1]
    bam_files = list(region_bam_dir.glob(f"{sample_name}*.bam"))
    for bam in bam_files:
        bed_output = region_bed_dir / f"{bam.stem}.bed"
        cmd = [
            "RLM", "-b", str(bam),
            "-r", str(ref_dir / f"{genome}.fa"),
            "-m", mode, "-s", "all",
            "-a", aligner, "-o", str(bed_output)
        ]
        run_command(cmd, f"Running RLM for {bam.name}")
    
    # Step 3: Calculate ReadsRatio
    reads_ratio_file = output_subdirs[2] / f"{sample_name}.txt"
    with open(reference_list, "r") as f, open(reads_ratio_file, "w") as out_f:
        for line in f:
            cell_type = line.strip()
            if not cell_type:
                continue
            
            # Find corresponding BED file
            bed_files = list(region_bed_dir.glob(f"{sample_name}_{cell_type}*.bed"))
            if not bed_files:
                logger.warning(f"No BED files found for {cell_type}")
                out_f.write(f"{cell_type}\t0.0000\n")
                continue
            
            # Calculate count/total with awk (wrap in bash for pipeline)
            count_cmd = f"""sed '1d' {" ".join(str(b) for b in bed_files)} | awk '{{if($6>=3 &&$NF<=0.25) print $0}}' | wc -l"""
            total_cmd = f"""sed '1d' {" ".join(str(b) for b in bed_files)} | awk '{{if($6>=3) print $0}}' | wc -l"""
            
            count = int(subprocess.check_output(count_cmd, shell=True).strip())
            total = int(subprocess.check_output(total_cmd, shell=True).strip())
            
            ratio = round(count / total, 4) if total > 0 else 0.0
            out_f.write(f"{cell_type}\t{ratio}\n")
    
    # Step 4: Calculate RPKM
    rpkm_file = output_subdirs[3] / f"{sample_name}.txt"
    with open(reference_list, "r") as f, open(rpkm_file, "w") as out_f:
        for line in f:
            cell_type = line.strip()
            if not cell_type:
                continue
            
            bed_files = list(region_bed_dir.glob(f"{sample_name}*{cell_type}*.bed"))
            marker_bed_file = marker_bed(marker_dir, cell_type, genome, marker_set)
            
            if not bed_files or not marker_bed_file.exists():
                logger.warning(f"Missing files for RPKM calculation: {cell_type}")
                out_f.write(f"{cell_type}\t0.000000\n")
                continue
            
            # Calculate region length
            region_len = 0
            with open(marker_bed_file, "r") as marker_handle:
                for marker_line in marker_handle:
                    fields = marker_line.rstrip("\n").split("\t")
                    if len(fields) >= 3:
                        region_len += int(fields[2]) - int(fields[1])
            
            # Calculate count/total
            count_cmd = f"""sed '1d' {" ".join(str(b) for b in bed_files)} | awk '{{if($6>=3 && $NF<=0.25) print $0}}' | wc -l"""
            total_cmd = f"""sed '1d' {" ".join(str(b) for b in bed_files)} | awk '{{if($6>=3) print $0}}' | wc -l"""
            
            count = int(subprocess.check_output(count_cmd, shell=True).strip())
            total = int(subprocess.check_output(total_cmd, shell=True).strip())
            
            # Calculate RPKM
            if total == 0 or region_len == 0:
                rpkm = 0.0
            else:
                per_k = region_len / 1000
                per_m = total / 1000000
                denom = per_k * per_m
                rpkm = round(count / denom, 6) if denom > 0 else 0.0
            
            out_f.write(f"{cell_type}\t{rpkm}\n")
    
    logger.info(f"Short-read deconvolution complete. Results in {output_dir}")
