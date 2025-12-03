# Quick Start Guide

## 1. Personalize the Files

Before pushing to GitHub, update these files with your information:

**setup.py** (lines 8-10):
```python
author="Your Name",
author_email="your.email@example.com",
url="https://github.com/YOUR_USERNAME/prospector-dual-region",
```

**prospector_dual_region/__init__.py** (line 8):
```python
__author__ = "Your Name"
```

**LICENSE** (line 3):
```
Copyright (c) 2024 [Your Name]
```

## 2. Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `prospector-dual-region`
3. Description: "Simultaneous SED fitting of spatially resolved regions with Prospector"
4. Choose: **Public**
5. Do NOT initialize with README
6. Click "Create repository"

## 3. Push Your Code

```bash
cd prospector-dual-region

git init
git add .
git commit -m "Initial commit: Dual-region SED fitting for Prospector"
git remote add origin https://github.com/YOUR_USERNAME/prospector-dual-region.git
git branch -M main
git push -u origin main
```

## 4. Configure Repository

On GitHub:
- Click Settings → About (gear icon)
- Add topics: `prospector`, `sed-fitting`, `astronomy`, `stellar-populations`
- Add website/description

## 5. Make First Release

1. Go to Releases → "Create a new release"
2. Tag: `v0.1.0`
3. Title: `v0.1.0 - Initial Release`
4. Description:
   ```
   Initial release of dual-region SED fitting for Prospector.
   
   Features:
   - Simultaneous fitting of two regions
   - Blended photometry constraint
   - Full Bayesian inference with dynesty
   ```
5. Publish release

## 6. Share with Community

Post in [Prospector Discussions](https://github.com/bd-j/prospector/discussions):

**Title:** Extension: Dual-Region SED Fitting with Blended Photometry

**Body:**
```
Hi everyone,

I've developed an extension for Prospector that enables simultaneous 
SED fitting of two spatially resolved regions with blended photometry.

Use cases:
- Lensed galaxies with multiple images
- Star-forming clumps in high-z galaxies
- Merging systems

Repository: https://github.com/YOUR_USERNAME/prospector-dual-region

Feedback welcome!
```

## That's It!

Your repository is ready to share. See `docs/README_DUAL_REGION.md` for detailed usage instructions.

## Optional Enhancements

- Add DOI via Zenodo: https://zenodo.org/
- Create example Jupyter notebook
- Add unit tests
- Set up GitHub Actions for CI
