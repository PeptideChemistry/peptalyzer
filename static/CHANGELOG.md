# Changelog
All notable changes to **Peptalyzer™** will be documented in this file.  
This project follows [Semantic Versioning](https://semver.org/).

---

## [1.4.0] - 2026-01-05
### Added
- **Extinction Coefficient (ε205)**: Added calculation for extinction coefficient at 205 nm to complement the existing 280 nm metric.
- **Amino Acid Distribution Visualization**: Enhanced the Amino Acid Counts table with a new "Nature & Distribution" column featuring color-coded visual bars representing residue abundance and chemical nature.

### Changed
- **UI Enhancements**: Added tooltips and direct resource links for Extinction Coefficients.
- **Color Updates**: Adjusted the color for "Others" (Glycine/Proline) in the distribution visualization for better contrast.
- **Reporting**: Updated the PDF report to include the new ε205 metric and the enhanced amino acid distribution table.

## [1.3.0] - 2026-01-05
### Added
- **Neutral pH Aromaticity**: Implemented a new sub-metric for Aromaticity that accounts for Histidine (His) contribution at neutral pH.

### Changed
- **UI/UX Refactoring**: Reorganized the results table to better categorize properties.
  - Moved "Aromaticity" to the "Molecular Characteristics" section.
  - Renamed the "Indices" section to "Bio-Predictive Metrics".
- **Resource Linking**: Added direct links to educational resources for Aromaticity, Aliphatic Index, and Instability Index.
- **Reporting**: Updated the PDF report template to match the new dashboard layout and include the Histidine-inclusive aromaticity metric.

## [1.2.2] - 2026-01-02
### Fixed
- **Deployment**: Bumped static asset versions to force browser cache refresh and ensure all users receive the latest JavaScript logic for report generation.

## [1.2.1] - 2026-01-02
### Fixed
- **Critical Chemistry Logic**: Corrected `calculate_net_charge` to include Cysteine (C) as an ionizable acidic residue. Previous versions treated it as neutral in specific calculation paths.
- **Concurrency/State Management**: Replaced global state with a unique ID system (`CALCULATION_CACHE`) to prevent users from overwriting each other's reports.
- **Input Validation**: Fixed a bug in `input_validator.py` where a premature return caused most validation checks (pH range, Cysteine counts) to be skipped.
- **PDF Report**: Fixed missing "Selected pH" value and corrected the Amino Acid Composition table to properly display counts and percentages.
- **Frontend Performance**: Fixed an issue where sorting event listeners were duplicated on every calculation, and optimized disulfide bond constraint logic.

### Changed
- **Code Architecture**: Refactored `app.py` to separate business logic from request handling.
- **Aliphatic Index**: Updated formula to strictly match the Ikai (1980) standard.
- **Hydropathy Plots**: Standardized smoothing logic for Hopp-Woods plots to handle edge cases consistently with Kyte-Doolittle.

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