# TempEdge: Artefact Detection and Filtering in Landsat Surface Temperature Data

**Author**: Gulam Mohiuddin (2025)  
**License**: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)  
**Repository Purpose**: Methodological transparency and reproducibility for the future users.

---

## 🔍 Description

**TempEdge** is a methodological framework designed to detect and filter artefacts in Landsat Collection 2 Level-2 (LC2L2) Surface Temperature (ST) data, particularly focused on tropical urban environments. The approach ensures high-quality LST observations by combining:

- Cloud masking via CFMask
- Monthly LST plausibility checks
- Threshold derivation using TempEdge logic
- Artefact filtering with statistical comparisons

This repository is not intended for end-users or application interfaces, but instead for scientific reviewers and researchers assessing the reproducibility and robustness of the method.

---
<details open>
<summary><strong>📁 Repository Structure</strong></summary>

tempedge/
├── LICENSE ← License file (CC BY-NC 4.0)
├── .gitignore ← Ignored files
├── README.md ← This file
├── requirements.txt ← Required Python packages
└── modules/ ← Main codebase
├── apply_tempedge/ ← Final application of TempEdge thresholds
│ └── apply_tempedge.py
├── dev_analysis_tempedge/ ← Development and statistical evaluation
│ └── dev_analysis_tempedge.py
└── cross_site_validation/ ← Cross-site testing across tropical cities
└── cross_site_validation.py

```bash
</details>

---

## ⚙️ Requirements

Install required libraries using:

```bash
pip install -r requirements.txt

You must also authenticate and initialize Google Earth Engine for Python.

🧪 Scripts Overview
dev_analysis_tempedge.py
Full statistical workflow

TempEdge threshold derivation

Artefact comparison with alternative thresholding methods

Includes visualizations from manuscript

apply_tempedge.py
Loads GEE data

Applies TempEdge filtering to retain only plausible LST values

cross_site_validation.py
Evaluates TempEdge in multiple tropical cities (e.g., Lagos, Mérida, Kuala Lumpur)

Generates comparison figure for raw vs. filtered LST values

📌 Notes
Default ROI: Phnom Penh, Cambodia (editable)

EPSG codes are defined per location

Ensure your GEE account has access to LC2L2 datasets

🔒 License
This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license.

You may share and adapt the material for non-commercial purposes with appropriate credit.

📫 Contact
For queries, please contact: Gulam Mohiuddin at: Gulam.Mohiuddin@hnee.de
