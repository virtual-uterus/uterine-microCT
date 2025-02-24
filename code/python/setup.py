from setuptools import setup, find_packages
from pathlib import Path


# Function to read requirements.txt
def read_requirements():
    requirements_file = Path(__file__).parent / "requirements.txt"
    with requirements_file.open(encoding="utf-8") as f:
        return f.read().splitlines()


setup(
    name="thickness",
    version="0.1.0",
    description="A Python package to compute the thickness of a uterine"
    "micro-computed tomography dataset",
    long_description_content_type="text/markdown",  # Use Markdown for README
    author="Mathias Roesler",
    author_email="mroe734@aucklanduni.ac.nz",
    url="https://github.com/virtual-uterus/uterine-microCT",
    license="Apache License 2.0",  # License type
    packages=find_packages(),
    python_requires=">=3.8",  # Minimum Python version required
    install_requires=read_requirements(),
    zip_safe=True,
)
