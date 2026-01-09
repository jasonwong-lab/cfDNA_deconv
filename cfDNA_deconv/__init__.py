"""
cfDNA_deconv: A Python package for deconvoluting cfDNA composition from methylation data
Supports short-read and Nanopore sequencing data.
"""
__version__ = "0.1.0"
__author__ = "Qiuyu Jing"

from .cli import cli
from .short_read import deconvolute_short_read
from .nanopore import deconvolute_nanopore

__all__ = ["cli", "deconvolute_short_read", "deconvolute_nanopore"]
