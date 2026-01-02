# Changelog
All notable changes to **Peptalyzer™** will be documented in this file.  
This project follows [Semantic Versioning](https://semver.org/).

---

## [1.2.0] - 2026-01-02
### Added
- **Multi-Scale pI Comparison**
  - Integrated support for four distinct pKa scales: **IPC 2.0** (default), **Bjellqvist**, **EMBOSS**, and **Lehninger**.
  - Enabled side-by-side comparison of results to facilitate alignment with varied laboratory and synthesis standards.
- **Educational Integration**
  - Added direct links from the dashboard results (pI and Net Charge) to the comprehensive guide on **Peptide pI Calculation** for better user context and theoretical transparency.

### Improved
- **Isoelectric Point (pI) Accuracy**
  - Enhanced the bisection algorithm to support highly acidic modified peptides by extending the search range down to **pH 0.0**.
  - Decoupled side-chain ionization from terminal ionization logic to improve accuracy for modified sequences.
- **Charge Calculation Logic**
  - Refined terminal handling for **Acetylation (Ac-)** and **Amidation (-NH2)** to ensure correct charge suppression in synthetic peptides.

### Fixed
- **The "0.000 pI" Bug**
  - Resolved a calculation error where N-terminally modified acidic sequences (e.g., Ac-HE) incorrectly defaulted to a pI of 0.000.
- **Data Synchronization**
  - Fixed a state-management issue that prevented comparison scales from updating in real-time on the dashboard.
- **UI/UX Refinements**
  - Corrected table layout breaks and ensured responsive rendering of nested results on mobile and desktop views.

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