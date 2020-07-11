from setuptools import setup

import sys
import os

import cigarco

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

setup(
    name="CIGARCO",
    version=cigarco.version,
    author="Sergey Aganezov",
    author_email="sergeyaganezovjr@gmail.com",
    description="A tool for computationally efficient transformation of coordinates between query and target aligned sequences",
    license="MIT",
    keywords="Bioinformatics, CIGAR, alignments, coordinate transformation",
    url="https://github.com/aganezov/rck",
    packages=["", "cigarco"],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "cigarco = cigarco.app:main",
        ]
    },
    install_requires=["pytest", "hypothesis"]
)
