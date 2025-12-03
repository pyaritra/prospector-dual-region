from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="prospector-dual-region",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Simultaneous SED fitting of two spatially resolved regions with Prospector",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/YOUR_USERNAME/prospector-dual-region",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.18",
        "scipy>=1.5",
        "astro-prospector>=1.0",
        "python-fsps>=0.4",
        "dynesty>=1.1",
        "sedpy>=0.2",
        "h5py>=3.0",
        "astropy>=4.0",
        "matplotlib>=3.3",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "corner",
            "jupyter",
        ],
    },
)
