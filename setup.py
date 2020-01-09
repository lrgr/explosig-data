import os
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="explosig-data",
    version="0.0.3",
    author="Leiserson Research Group",
    description="Process mutation data into standard formats originally developed for the ExploSig family of tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lrgr/explosig-data",
    packages=setuptools.find_packages(),
    package_data={
        'explosig_data': [
            os.path.join('snakefiles', 'genes', 'human.smk'),
            os.path.join('snakefiles', 'genomes', 'human.smk')
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'requests>=2.22.0',
        'pandas>=0.25.1',
        'numpy>=1.17.0',
        'snakemake>=5.3',
        'biopython>=1.75',
        'twobitreader>=3.1',
        'tqdm>=4.39.0'
    ],
)
