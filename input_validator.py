# input_validator.py

from peptide_helpers import N_TERM_MODIFICATIONS, C_TERM_MODIFICATIONS, validate_sequence

def validate_input(data):
    """
    Validates the user input for the peptide calculator.
    Returns a tuple: (is_valid: bool, message: str)
    """
    # ===== Validate peptide sequence =====
    sequence = data.get("sequence", "").strip().upper()
    if not sequence:
        return False, "Peptide sequence is required."

    try:
        validate_sequence(sequence)
    except ValueError as ve:
        return False, str(ve)

    # ===== Extract disulfide bond inputs =====
    try:
        intra_disulfide_bonds = int(data.get("intra_disulfide_bonds", 0))
        inter_disulfide_bonds = int(data.get("inter_disulfide_bonds", 0))
    except Exception:
        return False, "Disulfide bond counts must be integers."

    if intra_disulfide_bonds < 0 or inter_disulfide_bonds < 0:
        return False, "Disulfide bond counts cannot be negative."

    # ===== Count cysteines in sequence =====
    cys_count = sequence.count("C")

    # ===== Disallow disulfide bonds if no cysteines =====
    if cys_count == 0 and (intra_disulfide_bonds > 0 or inter_disulfide_bonds > 0):
        return False, "Disulfide bonds cannot be present when there are no cysteine residues in the sequence."

    # ===== Validate cysteine usage =====
    cysteines_used = (intra_disulfide_bonds * 2) + (inter_disulfide_bonds * 1)

    if cysteines_used > cys_count:
        return False, f"Total cysteines required exceed available cysteines. You have {cys_count} Cys but require {cysteines_used}."

    # ===== Validate pH values =====
    try:
        pH_min = float(data.get("pH_min", 0.0))
        pH_max = float(data.get("pH_max", 14.0))
        pH_step = float(data.get("pH_step", 0.1))
    except Exception:
        return False, "pH values must be valid numbers."

    if not (0.0 <= pH_min <= 14.0 and 0.0 <= pH_max <= 14.0):
        return False, "pH min and max must be between 0 and 14."

    if pH_min >= pH_max:
        return False, "pH min must be smaller than pH max."

    if pH_step <= 0:
        return False, "pH step must be greater than 0."

    # ===== Validate terminal modifications =====
    n_term = data.get("n_term", "H")
    c_term = data.get("c_term", "OH")
    if n_term not in N_TERM_MODIFICATIONS:
        return False, f"Unknown N-terminal modification: {n_term}"
    if c_term not in C_TERM_MODIFICATIONS:
        return False, f"Unknown C-terminal modification: {c_term}"

    # ===== Validate hydropathy window =====
    window_size = int(data.get("hydropathy_window", 5))
    if window_size <= 0:
        return False, "Hydropathy window must be greater than 0."

    return True, "Input is valid."