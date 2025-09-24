# Changelog
All notable changes to **Peptalyzer™** will be documented in this file.  
This project follows [Semantic Versioning](https://semver.org/).

---

## [1.1.0] - 2025-09-21
### Added
- **Amino Acid Composition**
  - Percentage column added to the Amino Acid Counts table (now always visible, even before calculation)
- **Descriptor Links**
  - GRAVY Score, Kyte–Doolittle hydropathy plot, and Hopp–Woods hydrophilicity plot now link to dedicated explanatory articles on peptidechemistry.org
- **User Guidance**
  - Added explanatory note *(i)* to clarify the influence of the "Select pH" slider on results
- **Plots**
  - Net Charge vs pH plot now indicates when no charged amino acids are present, improving interpretation

### Improved
- **Hydropathy Profile (Kyte–Doolittle)**
  - Refinement of calculation and smoothing approach for more consistent representation of hydropathy profiles  
  *(previous version was functional but less aligned with reference standards)*

### Fixed
- Tooltip hover text alignment issues
- Minor favicon path and stylesheet reference corrections
- Corrected SEO keyword handling in meta tags
- Corrected Hopp–Woods hydrophilicity scale dictionary: **Thr (T)** value updated from `+0.4` ➜ `-0.4`, now consistent with the original Hopp & Woods (1981) publication and ExPASy reference implementation.
- Updated default window size in `smooth_hopp_woods_hydropathy` from **5 ➜ 6**.

---

## [1.0.0] - 2025-07-01
### Initial Release
- Core peptide property calculations:
  - Molecular weight (monoisotopic & average)
  - Net charge & isoelectric point (pI)
  - Hydropathy index (Kyte–Doolittle)
  - Boman Index, Aliphatic Index
  - Molecular formula & amino acid composition
- Results export to CSV
- Responsive frontend layout with real-time calculation