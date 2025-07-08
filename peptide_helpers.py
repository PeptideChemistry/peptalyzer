from collections import Counter
import numpy as np
import re

# DICTIONARIES

# 🔢 Monoisotopic mass (residue only, water already removed, in Da)
MONOISOTOPIC_MASS = {
    'A': 71.03711, 'R': 156.10110, 'N': 114.04292, 'D': 115.02693, 'C': 103.009180,
    'E': 129.04258, 'Q': 128.05857, 'G': 57.02146, 'H': 137.05890, 'I': 113.08405,
    'L': 113.08405, 'K': 128.09495, 'M': 131.04048, 'F': 147.06840, 'P': 97.05276,
    'S': 87.03202, 'T': 101.04767, 'W': 186.07930, 'Y': 163.06332, 'V': 99.06840
}

# 🔢 Average mass (residue only, water already removed, in Da)
AVERAGE_MASS = {
    'A': 71.07884,  'R': 156.18764, 'N': 114.10392, 'D': 115.08864, 'C': 103.14484,
    'E': 129.11552, 'Q': 128.13080, 'G': 57.05196,  'H': 137.14120, 'I': 113.15947,
    'L': 113.15947, 'K': 128.17416, 'M': 131.19860, 'F': 147.17660, 'P': 97.11672,
    'S': 87.07824,  'T': 101.10512, 'W': 186.21327, 'Y': 163.17600, 'V': 99.13259
}

# 🧬 Valid standard amino acids
VALID_AMINO_ACIDS = set(MONOISOTOPIC_MASS.keys())

# ⚛️ Residue-specific ionizable pKa values for net charge/pI
AA_PROPERTIES = {
    'D': {'pKa': 3.86,  'acidic': True},
    'E': {'pKa': 4.25,  'acidic': True},
    'C': {'pKa': 8.33,  'acidic': True},
    'Y': {'pKa': 10.07, 'acidic': True},
    'H': {'pKa': 6.00,  'basic': True},
    'K': {'pKa': 10.53, 'basic': True},
    'R': {'pKa': 12.48, 'basic': True}
}

# ⚛️ N-terminal modifications and their charge behavior
N_TERM_MODIFICATIONS = {
    "H": {
        "monoisotopic_mass_shift": 1.00782,
        "average_mass_shift": 1.00794,
        "label": "H- (free amine)",
        "ionizable_groups": [
            {
                "pKa": 8.6,
                "acidic": False,
                "charge": 1
            }
        ],
        "atomic_composition": {'H': 1}  # Free amine group (NH2 is in amino acid formulas)
    },
    "Ac": {
        "monoisotopic_mass_shift": 43.01838,
        "average_mass_shift": 43.04522,
        "label": "Ac- (acetylated)",
        "ionizable_groups": [],
        "atomic_composition": {'C': 2, 'H': 3, 'O': 1}  # CH3-CO group
    }
}

# ⚛️ C-terminal modifications and their charge behavior
C_TERM_MODIFICATIONS = {
    "OH": {
        "monoisotopic_mass_shift": 17.00273,
        "average_mass_shift": 17.00734,
        "label": "-OH (carboxylic acid)",
        "ionizable_groups": [
            {
                "pKa": 3.6,
                "acidic": True,
                "charge": 1
            }
        ],
        "atomic_composition": {'O': 1, 'H': 1}
    },
    "NH2": {
        "monoisotopic_mass_shift": 16.01872,
        "average_mass_shift": 16.02262,
        "label": "-NH2 (amide)",
        "ionizable_groups": [],
        "atomic_composition": {'N': 1, 'H': 2}

    }
}

# 📊 Secondary structure scoring (Chou-Fasman)
CHOU_FASMAN = {
    'A': {'helix': 1.42, 'sheet': 0.83, 'coil': 0.66},
    'R': {'helix': 0.98, 'sheet': 0.93, 'coil': 0.95},
    'N': {'helix': 0.67, 'sheet': 0.89, 'coil': 1.56},
    'D': {'helix': 1.01, 'sheet': 0.54, 'coil': 1.46},
    'C': {'helix': 0.70, 'sheet': 1.19, 'coil': 1.19},
    'Q': {'helix': 1.11, 'sheet': 1.10, 'coil': 0.98},
    'E': {'helix': 1.51, 'sheet': 0.37, 'coil': 1.00},
    'G': {'helix': 0.57, 'sheet': 0.75, 'coil': 1.64},
    'H': {'helix': 1.00, 'sheet': 0.87, 'coil': 0.95},
    'I': {'helix': 1.08, 'sheet': 1.60, 'coil': 0.47},
    'L': {'helix': 1.21, 'sheet': 1.30, 'coil': 0.59},
    'K': {'helix': 1.16, 'sheet': 0.74, 'coil': 1.01},
    'M': {'helix': 1.45, 'sheet': 1.05, 'coil': 0.60},
    'F': {'helix': 1.13, 'sheet': 1.38, 'coil': 0.60},
    'P': {'helix': 0.57, 'sheet': 0.55, 'coil': 1.52},
    'S': {'helix': 0.77, 'sheet': 0.75, 'coil': 1.43},
    'T': {'helix': 0.83, 'sheet': 1.19, 'coil': 0.96},
    'W': {'helix': 1.08, 'sheet': 1.37, 'coil': 0.96},
    'Y': {'helix': 0.69, 'sheet': 1.47, 'coil': 1.14},
    'V': {'helix': 1.06, 'sheet': 1.70, 'coil': 0.50}
}

# 💧 Hydropathy scales - Kyte-Doolittle (Primary Scale for GRAVY and Plot)
HYDROPATHY = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5,  'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5,  'L': 3.8,  'K': -3.9,
    'M': 1.9,  'F': 2.8,  'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# 💧 Hydropathy scales - Hopp-Woods (Secondary Scale for Alternative Plot)
HYDROPATHY_HOPP_WOODS = {
    'A': -0.5, 'R': 3.0,  'N': 0.2,  'D': 3.0,  'C': -1.0,
    'E': 3.0,  'Q': 0.2,  'G': 0.0,  'H': -0.5, 'I': -1.8,
    'L': -1.8, 'K': 3.0,  'M': -1.3, 'F': -2.5, 'P': 0.0,
    'S': 0.3,  'T': 0.4,  'W': -3.4, 'Y': -2.3, 'V': -1.5
}

# 🔌 Binding affinity estimates (Boman Index, kcal/mol, from Kyte-Doolittle)
AA_BINDING_VALUES = {
    'A': 0.17, 'R': 0.81, 'N': 0.42, 'D': 0.19,
    'C': 0.24, 'Q': 0.58, 'E': 0.21, 'G': 0.01,
    'H': 0.70, 'I': 0.73, 'L': 0.53, 'K': 0.99,
    'M': 0.44, 'F': 0.61, 'P': -0.07, 'S': 0.13,
    'T': 0.14, 'W': 1.02, 'Y': 0.45, 'V': 0.54
}

# 🔥 Aliphatic index from hydrophobic side chains (Ikai, 1980)
ALIPHATIC_WEIGHTS = {
    'A': 0.0, 'V': 1.0, 'I': 2.9, 'L': 3.9
}

# ⚙️ Thresholds for hydropathy interpretation (used in GRAVY score)
GRAVY_THRESHOLDS = {
    "hydrophilic_limit": -0.5,
    "hydrophobic_limit": 0.5,
    "charged_residues_min": 5,
    "charged_residues_max": 2
}

# 🎨 Amino acid color codes (PepCalc-inspired)
RESIDUE_COLORS = {
    'D': '#ffb6c1',  # Acidic - Light Pink
    'E': '#ffb6c1',  # Acidic - Light Pink
    'R': '#1e90ff',  # Basic - Dodger Blue
    'K': '#1e90ff',  # Basic - Dodger Blue
    'H': '#1e90ff',  # Basic - Dodger Blue
    'F': '#90ee90',  # Aromatic - Light Green
    'W': '#90ee90',  # Aromatic - Light Green
    'Y': '#90ee90',  # Aromatic - Light Green
    'A': '#d3d3d3',  # Aliphatic - Light Gray
    'I': '#d3d3d3',  # Aliphatic - Light Gray
    'L': '#d3d3d3',  # Aliphatic - Light Gray
    'V': '#d3d3d3',  # Aliphatic - Light Gray
    'M': '#d3d3d3',  # Aliphatic - Light Gray
    'S': '#228B22',  # Polar - Forest Green
    'T': '#228B22',  # Polar - Forest Green
    'N': '#228B22',  # Polar - Forest Green
    'Q': '#228B22',  # Polar - Forest Green
    'C': '#ffd700',  # Cysteine - Gold
    'G': '#333333',  # Others - Neutral Gray
    'P': '#333333'   # Others - Neutral Gray
}

# ================================
# 🔧 FUNCTION DEFINITIONS
# ================================

# 🧬 Checks for invalid amino acids in input sequence
def validate_sequence(sequence):
    sequence = sequence.upper()
    invalid = [aa for aa in sequence if aa not in VALID_AMINO_ACIDS]
    if invalid:
        raise ValueError(f"Invalid residues in sequence: {invalid}")
    return sequence

# ⚡ Computes net charge at specified pH with terminal mods
def calculate_net_charge(sequence, pH, n_term="H", c_term="OH"):
    sequence = validate_sequence(sequence)
    aa_counts = Counter(sequence)
    charge = 0.0
    charge += apply_terminal_modification_charge(pH, N_TERM_MODIFICATIONS[n_term], is_n_term=True)
    charge += apply_terminal_modification_charge(pH, C_TERM_MODIFICATIONS[c_term], is_n_term=False)
    for aa, props in AA_PROPERTIES.items():
        count = aa_counts.get(aa, 0)
        if props.get("basic"):
            charge += count * (10 ** props["pKa"]) / (10 ** props["pKa"] + 10 ** pH)
        if props.get("acidic"):
            charge -= count * (10 ** pH) / (10 ** props["pKa"] + 10 ** pH)
    return charge

# 🧮 Computes average molecular weight with terminal mods
def calculate_average_mass(sequence, n_term="H", c_term="OH"):
    sequence = validate_sequence(sequence)
    base_mass = sum(AVERAGE_MASS.get(aa, 0.0) for aa in sequence)
    n_shift = N_TERM_MODIFICATIONS[n_term]["average_mass_shift"]
    c_shift = C_TERM_MODIFICATIONS[c_term]["average_mass_shift"]
    return round(base_mass + n_shift + c_shift, 4)

# 🧮 Computes monoisotopic mass with terminal mods
def calculate_monoisotopic_mass(sequence, n_term="H", c_term="OH"):
    base_mass = sum(MONOISOTOPIC_MASS.get(aa, 0) for aa in sequence)
    n_shift = N_TERM_MODIFICATIONS[n_term]["monoisotopic_mass_shift"]
    c_shift = C_TERM_MODIFICATIONS[c_term]["monoisotopic_mass_shift"]
    return round(base_mass + n_shift + c_shift, 4)

# 📈 Returns smoothed hydropathy profile using sliding window
def smooth_hydropathy(sequence, window_size=5):
    sequence = validate_sequence(sequence)
    if window_size < 1:
        window_size = 1

    values = np.array([HYDROPATHY.get(aa, 0) for aa in sequence])

    if window_size > len(sequence):
        # When the window is larger than the sequence, use the average value
        avg_value = np.mean(values) if len(values) > 0 else 0
        smoothed = np.array([avg_value] * len(sequence))
    else:
        kernel = np.ones(window_size) / window_size
        smoothed = np.convolve(values, kernel, mode='same')

    return [round(v, 2) for v in smoothed]

# 📈 Returns smoothed Hopp-Woods hydropathy profile using sliding window
def smooth_hopp_woods_hydropathy(sequence, window_size=5):
    sequence = validate_sequence(sequence)
    if window_size < 1:
        window_size = 1

    values = np.array([HYDROPATHY_HOPP_WOODS.get(aa, 0) for aa in sequence])

    if window_size > len(sequence):
        # When the window is larger than the sequence, use the average value
        avg_value = np.mean(values) if len(values) > 0 else 0
        smoothed = np.array([avg_value] * len(sequence))
    else:
        kernel = np.ones(window_size) / window_size
        smoothed = np.convolve(values, kernel, mode='same')

    return [round(v, 2) for v in smoothed]

# 🌊 Labels peptide as polar/nonpolar based on hydropathy and charge
def classify_peptide_polarity(sequence):
    sequence = validate_sequence(sequence)
    hydrophobicity = sum(HYDROPATHY.get(aa, 0) for aa in sequence)
    charged_residues = sum(1 for aa in sequence if aa in {'D', 'E', 'R', 'K', 'H'})
    length = len(sequence)
    charge_fraction = charged_residues / length if length > 0 else 0

    if hydrophobicity < 0 and charge_fraction >= 0.20:
        return "polar"
    elif hydrophobicity > 20 and charge_fraction <= 0.05:
        return "nonpolar"
    else:
        return "intermediate"

#  🧠 Predicts structural bias (secondary structure) using Chou-Fasman propensities
def estimate_secondary_structure(sequence):
    sequence = validate_sequence(sequence)
    counts = {'helix': 0.0, 'sheet': 0.0, 'coil': 0.0}
    for aa in sequence:
        if aa in CHOU_FASMAN:
            for s in counts:
                counts[s] += CHOU_FASMAN[aa][s]
    total = sum(counts.values())
    if total == 0:
        return {s: 0.0 for s in counts}
    return {s: round(100 * counts[s] / total, 2) for s in counts}

# 💥Predicts protein-binding propensity (Boman Index)
def calculate_boman_index(sequence):
    sequence = validate_sequence(sequence)
    if not sequence:
        return "0.000 kcal/mol (N/A)"
    boman = sum(AA_BINDING_VALUES.get(aa, 0) for aa in sequence) / len(sequence)
    rounded = round(boman, 3)
    if rounded < 1.0:
        rating = "Low"
    elif rounded < 2.5:
        rating = "Moderate"
    else:
        rating = "High"
    return f"{rounded:.3f} kcal/mol ({rating})"

# 🔥 Calculates aliphatic index using side-chain weights
def calculate_aliphatic_index(sequence):
    sequence = validate_sequence(sequence)
    if not sequence:
        return {"value": 0.0, "rating": "N/A"}

    length = len(sequence)
    counts = Counter(sequence)

    X_V = counts.get('V', 0) / length
    X_I = counts.get('I', 0) / length
    X_L = counts.get('L', 0) / length
    X_M = counts.get('M', 0) / length

    aliphatic_index = round(100 * (X_V * 1.0 + X_I * 2.9 + X_L * 3.9 + X_M * 1.9), 2)

    if aliphatic_index < 50:
        rating = "Low"
    elif aliphatic_index < 80:
        rating = "Moderate"
    elif aliphatic_index <= 100:
        rating = "High"
    else:
        rating = "Very High"

    return {"value": aliphatic_index, "rating": rating}

# 💧 Calculates GRAVY score and hydropathy annotation
def calculate_gravy_score(sequence):
    sequence = validate_sequence(sequence)
    if not sequence:
        return {"value": 0.0, "annotation": "N/A"}
    gravy = sum(HYDROPATHY.get(aa, 0) for aa in sequence) / len(sequence)
    rounded = round(gravy, 3)
    charged_residues = sum(1 for aa in sequence if aa in {'D', 'E', 'R', 'K', 'H'})

    # Use config values
    t = GRAVY_THRESHOLDS
    if rounded < t["hydrophilic_limit"] and charged_residues > t["charged_residues_min"]:
        annotation = "Strongly hydrophilic and charged"
    elif rounded > t["hydrophobic_limit"] and charged_residues < t["charged_residues_max"]:
        annotation = "Strongly hydrophobic and neutral"
    else:
        annotation = "Mixed hydropathy/charge behavior"
    return {"value": rounded, "annotation": annotation}


# 🧪 Empirical atomic formulas for peptide composition
def calculate_molecular_formula(sequence, n_term="H", c_term="OH"):
    sequence = validate_sequence(sequence)
    amino_acid_formulas = {
        'A': {'C':3,'H':5,'N':1,'O':1}, 'R': {'C':6,'H':12,'N':4,'O':1}, 'N': {'C':4,'H':6,'N':2,'O':2},
        'D': {'C':4,'H':5,'N':1,'O':3}, 'C': {'C':3,'H':5,'N':1,'O':1,'S':1}, 'E': {'C':5,'H':7,'N':1,'O':3},
        'Q': {'C':5,'H':8,'N':2,'O':2}, 'G': {'C':2,'H':3,'N':1,'O':1}, 'H': {'C':6,'H':7,'N':3,'O':1},
        'I': {'C':6,'H':11,'N':1,'O':1}, 'L': {'C':6,'H':11,'N':1,'O':1}, 'K': {'C':6,'H':12,'N':2,'O':1},
        'M': {'C':5,'H':9,'N':1,'O':1,'S':1}, 'F': {'C':9,'H':9,'N':1,'O':1}, 'P': {'C':5,'H':7,'N':1,'O':1},
        'S': {'C':3,'H':5,'N':1,'O':2}, 'T': {'C':4,'H':7,'N':1,'O':2}, 'W': {'C':11,'H':10,'N':2,'O':1},
        'Y': {'C':9,'H':9,'N':1,'O':2}, 'V': {'C':5,'H':9,'N':1,'O':1}
    }

    formula_counts = Counter()

    # Add amino acid atoms
    for aa in sequence:
        for atom, count in amino_acid_formulas.get(aa, {}).items():
            formula_counts[atom] += count

    # Get terminal atom compositions from the modification dictionaries
    n_term_atoms = N_TERM_MODIFICATIONS.get(n_term, {}).get('atomic_composition', {})
    c_term_atoms = C_TERM_MODIFICATIONS.get(c_term, {}).get('atomic_composition', {})

    # Add terminal atoms
    for atom, count in n_term_atoms.items():
        formula_counts[atom] += count
    for atom, count in c_term_atoms.items():
        formula_counts[atom] += count

    formatted = ''.join(f"{atom}{formula_counts[atom]}" for atom in sorted(formula_counts))
    return {
        "raw": dict(formula_counts),
        "formatted": formatted
    }

# ⚛️ Computes terminal charge contributions at pH (from mod + removes native charge if needed)
def apply_terminal_modification_charge(pH, mod, is_n_term):
    """
    Calculate the charge contribution of an N- or C-terminal modification at a given pH.
    Args:
        pH (float): pH at which the calculation is made.
        mod (dict): Modification dictionary.
        is_n_term (bool): True if N-terminal, False if C-terminal.
    Returns:
        float: Charge contribution at the given pH.
    """
    total = 0.0
    for group in mod.get("ionizable_groups", []):
        pKa = group["pKa"]
        acidic = group["acidic"]
        charge = group.get("charge", 1)
        if acidic:
            # Acidic group (e.g., COOH)
            fraction = 10 ** pH / (10 ** pKa + 10 ** pH)
            total -= charge * fraction
        else:
            # Basic group (e.g., NH2)
            fraction = 10 ** pKa / (10 ** pKa + 10 ** pH)
            total += charge * fraction
    return total

# ⚖️ Estimates isoelectric point via bisection method
def calculate_isoelectric_point(sequence, n_term="H", c_term="OH"):
    """
    Estimate the isoelectric point (pI) of a peptide sequence using bisection.
    Args:
        sequence (str): Peptide sequence (1-letter code).
        n_term (str): N-terminal modification key.
        c_term (str): C-terminal modification key.
    Returns:
        float: Estimated isoelectric point (pI).
    """
    sequence = validate_sequence(sequence)
    low = 0.0
    high = 14.0
    tolerance = 0.001
    for _ in range(100):
        mid = (low + high) / 2
        charge = calculate_net_charge(sequence, pH=mid, n_term=n_term, c_term=c_term)
        if abs(charge) < tolerance:
            return round(mid, 3)
        if charge > 0:
            low = mid
        else:
            high = mid
    return round(mid, 3)  # Fallback if tolerance not reached

# 🧮 Counts the number of each amino acid in the sequence
def count_amino_acids(sequence):
    sequence = validate_sequence(sequence)
    counts = Counter(sequence)
    # Keep only amino acids with counts > 0
    filtered_counts = {aa: count for aa, count in counts.items() if count > 0}
    return filtered_counts

# ✅ String formatting utility 🔤
# Adds HTML subscripts to numerical characters in a string.
# Example: 'C6H12O6' ➜ 'C<sub>6</sub>H<sub>12</sub>O<sub>6</sub>'
def format_subscripts(text):
    """
    Adds HTML subscript tags to numbers in a string.

    :param text: The input string to format.
    :return: The formatted string with subscripts.
    """
    return re.sub(r'(\d+)', r'<sub>\1</sub>', text)

# ⚡ Charge distribution calculator ⚡
def calculate_charge_distribution(sequence):
    """
    Calculates the number of acidic and basic residues in a peptide sequence.

    :param sequence: Peptide sequence (1-letter code)
    :return: Dictionary with counts of acidic and basic residues
    """
    acidic_residues = {'D', 'E'}
    basic_residues = {'R', 'K', 'H'}
    return {
        "acidic_count": sum(sequence.count(aa) for aa in acidic_residues),
        "basic_count": sum(sequence.count(aa) for aa in basic_residues)
    }

def format_molecular_formula(atom_counts):
    """
    Formats the molecular formula dictionary into a standard string like C6H12O6.

    :param atom_counts: Dictionary of atoms and their counts.
    :return: Molecular formula as a string.
    """
    formula = ''
    for atom in sorted(atom_counts.keys()):
        count = atom_counts[atom]
        if count == 1:
            formula += f"{atom}"
        elif count > 1:
            formula += f"{atom}{count}"
    return formula

def calculate_extinction_coefficient(sequence, total_disulfide_bonds, peptide_units=1):
    """
    Calculates the extinction coefficient based on user-defined disulfide bonds and aromatic residues.
    Extinction coefficient is not scaled by the number of peptide units.
    """
    sequence = validate_sequence(sequence)
    trp_count = sequence.count('W')
    tyr_count = sequence.count('Y')

    # Extinction coefficient for the entire system (not per unit)
    extinction_total = (trp_count * 5500) + (tyr_count * 1490) + (total_disulfide_bonds * 125)

    return {
        "extinction_coefficient": extinction_total,
        "unit": "M⁻¹·cm⁻¹"
    }