# TempEdge: Artefact Detection and Filtering in Analysis-ready Landsat Surface Temperature Data

**Author of the code**: ***[Gulam Mohiuddin](https://www.linkedin.com/in/mohigeo33/)*** (2025)  
**License**: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)  
**Repository Purpose**: Methodological transparency and reproducibility for future users.

## 📑 1. Citation for the code

If you use TempEdge in your work, please cite the Zenodo archive:

> Mohiuddin, G. (2025). *TempEdge: Artefact Detection and Filtering in Landsat Analysis-Ready Surface Temperature Data* (v1.0.0) [Software].[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15543544.svg)](https://doi.org/10.5281/zenodo.15543544)

---

## 🔍 2. Description

**TempEdge** is a methodological framework designed to detect and filter artefacts in Landsat Collection 2 Level-2 (LC2L2) Surface Temperature (ST) data, particularly focused on tropical urban environments. The approach ensures high-quality LST observations by combining:

- Cloud masking via CFMask
- Monthly LST plausibility checks
- Threshold derivation using TempEdge logic
- Applying TempEdge Threshold to filter artefacts

## 📘 3. How to use the repository?
This repository has three modules:

***apply_tempedge:*** 
If you are interested in applying TempEdge in your study to detect and filter out artefacts, you should use this module

***dev_analysis_tempedge:*** 
If you are interested in seeing how TempEdge was developed and how statistical analysis is performed to reinforce the method's effectiveness

***cross_site_validation:*** 
If you are interested in seeing how the cross-site validation was performed for the TempEdge

---
## 📁 4. Repository Structure
 ```bash
tempedge/
├── LICENSE ← License file (CC BY-NC 4.0)
├── .gitignore ← Ignored files
├── README.md ← This file
├── requirements.txt ← Required Python packages
└── modules/ ← Main codebase
├── apply_tempedge/ ← ready code for the future TempEdge users
│ └── apply_tempedge.py
├── dev_analysis_tempedge/ ← development and statistical evaluation of the method
│ └── dev_analysis_tempedge.py
└── cross_site_validation/ ← cross-site testing across tropical cities
└── cross_site_validation.py
```
---

## ⚙️ 5. Requirements

Install required libraries using:

```bash
pip install -r requirements.txt

You must also authenticate and initialise Google Earth Engine for Python.
```
## 🧪 6. Scripts overview
***apply_tempedge.py***

Full implementation code for future TempEdge users
Loads Landsat image collections
Applies TempEdge filtering to retain only plausible LST values
Return the collections with filtered artefacts and further use

***dev_analysis_tempedge.py***

Full development and statistical workflow
TempEdge threshold derivation
Artefact comparison with alternative thresholding methods
Includes visualisations from the manuscript

***cross_site_validation.py***

Evaluates TempEdge in multiple tropical cities (e.g., Lagos, Mérida, Kuala Lumpur)
Generates comparison figure for raw vs. filtered LST values used in the manuscript

## 📊 7. Data overview
Data used for statistical analysis in the manuscript:

***dfallmonths.csv***

Data used for cross-site validations used for the manuscript

***dfallmonths_kuala.csv***

***dfallmonths_lagos.csv***

***dfallmonths_merida***

***dfallmonths_prudente.csv***

## 📌 8. Notes
Default ROI: Phnom Penh, Cambodia (editable)

EPSG codes are defined per location

## 🔒 9. License
This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license.
You may share and adapt the material for non-commercial purposes with appropriate credit as follows:

Mohiuddin, G. (2025). TempEdge: Artefact Detection and Filtering in Landsat Analysis-Ready Surface Temperature Data (v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.15543544

## 📫 10. Contact
For queries, please contact: Gulam Mohiuddin at: Gulam.Mohiuddin@hnee.de
