import click
import logging
from cfDNA_deconv.short_read import deconvolute_short_read
from cfDNA_deconv.nanopore import deconvolute_nanopore

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

@click.group()
def cli():
    """cfDNA_deconv: Deconvolute cfDNA composition from methylation data"""
    pass

@cli.command(name="short-read")
@click.option("--marker-dir", "-m", required=True, help="Path to marker directory (contains .bed files and ReferenceList.txt)")
@click.option("--ref-dir", "-r", required=True, help="Path to reference directory (contains hg19.fa)")
@click.option("--bam-file", "-b", required=True, help="Path to sorted/indexed BAM file")
@click.option("--output-dir", "-o", required=True, help="Path to output directory (will be created if missing)")
@click.option("--mode", "-M", default="PE", type=click.Choice(["PE", "SE"]), help="Sequencing mode (PE/SE)")
@click.option("--aligner", "-a", default="bismark", type=click.Choice(["bsmap", "bismark", "segemehl", "gem"]), help="Aligner used for short-read data")
@click.option("--cores", "-c", default=4, type=int, help="Number of threads to use")
def short_read_cli(marker_dir, ref_dir, bam_file, output_dir, mode, aligner, cores):
    """Deconvolute cfDNA composition from short-read sequencing data"""
    try:
        deconvolute_short_read(marker_dir, ref_dir, bam_file, output_dir, mode, aligner, cores)
        click.echo(f"Success! Results saved to {output_dir}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()

@cli.command(name="nanopore")
@click.option("--marker-dir", "-m", required=True, help="Path to marker directory (contains .bed files and ReferenceList.txt)")
@click.option("--ref-dir", "-r", required=True, help="Path to reference directory (contains hg19.fa)")
@click.option("--mbtools-path", "-p", required=True, help="Path to mbtools executable")
@click.option("--bam-file", "-b", required=True, help="Path to sorted/indexed BAM file (with methylation info)")
@click.option("--output-dir", "-o", required=True, help="Path to output directory (will be created if missing)")
@click.option("--cores", "-c", default=4, type=int, help="Number of threads to use")
def nanopore_cli(marker_dir, ref_dir, mbtools_path, bam_file, output_dir, cores):
    """Deconvolute cfDNA composition from Nanopore sequencing data"""
    try:
        deconvolute_nanopore(marker_dir, ref_dir, mbtools_path, bam_file, output_dir, cores)
        click.echo(f"Success! Results saved to {output_dir}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    cli()
