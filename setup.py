from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="cfDNA_deconv",
    version="1.0.0",
    author="Qiuyu Jing",
    author_email="qiuyu.jing@coinno.hk",
    description="Deconvolute cfDNA composition from methylation data (short-read/Nanopore)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jasonwong-lab/cfDNA_deconv",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.8,<3.12",
    install_requires=[
        "click>=8.0,<9.0",
        "pandas>=1.4,<2.0",
        "numpy>=1.21,<1.25",
    ],
    entry_points={
        "console_scripts": [
            "cfDNA_deconv = cfDNA_deconv.cli:cli",
        ],
    },
)