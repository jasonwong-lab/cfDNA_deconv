import os
import subprocess
import logging
from pathlib import Path
from cfDNA_deconv.utils import run_command, validate_inputs

logger = logging.getLogger(__name__)

def deconvolute_nanopore(
    marker_dir: str,
    ref_dir: str,
    mbtools_path: str,
    bam_file: str,
    output_dir: str,
    cores: int = 4
):
    """
    Deconvolute cfDNA composition from Nanopore sequencing data (BAM file)
    based on methylation markers.
    
    Args:
        marker_dir: Path to folder with .bed markers and ReferenceList.txt
        ref_dir: Path to folder with hg19.fa reference
        mbtools_path: Path to mbtools executable
        bam_file: Path to sorted/indexed BAM file (with methylation info)
        output_dir: Output directory for results
        cores: Number of threads, default=4
    """
    # Validate inputs
    validate_inputs([marker_dir, ref_dir, bam_file, mbtools_path], ["dir", "dir", "file", "file"])
    marker_dir = Path(marker_dir)
    ref_dir = Path(ref_dir)
    mbtools_path = Path(mbtools_path)
    bam_file = Path(bam_file)
    output_dir = Path(output_dir)
    
    # Preprocessing
    sample_name = bam_file.stem
    logger.info(f"Starting Nanopore deconvolution for sample: {sample_name}")
    
    # Create output directories
    output_subdirs = [
        output_dir / "Region_Methylome",
        output_dir / "ReadsRatio"
    ]
    for subdir in output_subdirs:
        subdir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Calculate read-level methylation with mbtools
    reference_list = marker_dir / "ReferenceList.txt"
    if not reference_list.exists():
        raise FileNotFoundError(f"ReferenceList.txt not found in {marker_dir}")
    
    with open(reference_list, "r") as f:
        for line in f:
            cell_type = line.strip()
            if not cell_type:
                continue
            bed_file = marker_dir / f"U250_hg19_{cell_type}.bed"
            if not bed_file.exists():
                logger.warning(f"Marker bed file missing: {bed_file}")
                continue
            
            output_tsv = output_subdirs[0] / f"{sample_name}_{cell_type}_methylome.tsv"
            cmd = [
                str(mbtools_path), "read-region-frequency",
                "-t", "0.1",
                "-r", str(bed_file),
                "-g", str(ref_dir / "hg19.fa"),
                str(bam_file)
            ]
            # Redirect output to TSV
            with open(output_tsv, "w") as out_f:
                run_command(cmd, f"Running mbtools for {cell_type}", stdout=out_f)
    
    # Step 2: Calculate ReadsRatio
    reads_ratio_file = output_subdirs[1] / f"{sample_name}.txt"
    with open(reference_list, "r") as f, open(reads_ratio_file, "w") as out_f:
        for line in f:
            cell_type = line.strip()
            if not cell_type:
                continue
            
            # Find corresponding TSV file
            tsv_files = list(output_subdirs[0].glob(f"{sample_name}_{cell_type}*_methylome.tsv"))
            if not tsv_files:
                logger.warning(f"No methylome TSV files found for {cell_type}")
                out_f.write(f"{cell_type}\t0.0000\n")
                continue
            
            # Calculate count/total with awk
            count_cmd = f"""sed '1d' {" ".join(str(t) for t in tsv_files)} | awk '{{if($8>=3 && $10<=0.25) print $0}}' | wc -l"""
            total_cmd = f"""sed '1d' {" ".join(str(t) for t in tsv_files)} | awk '{{if($8>=3) print $0}}' | wc -l"""
            
            count = int(subprocess.check_output(count_cmd, shell=True).strip())
            total = int(subprocess.check_output(total_cmd, shell=True).strip())
            
            ratio = round(count / total, 4) if total > 0 else 0.0
            out_f.write(f"{cell_type}\t{ratio}\n")
    
    logger.info(f"Nanopore deconvolution complete. Results in {output_dir}")
