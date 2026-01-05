from flask import Flask, Blueprint, request, jsonify, render_template
from flask_cors import CORS
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from datetime import datetime, timezone
from peptide_helpers import (
    validate_sequence,
    calculate_net_charge,
    calculate_average_mass,
    calculate_monoisotopic_mass,
    smooth_hydropathy,
    smooth_hopp_woods_hydropathy,
    classify_peptide_polarity,
    estimate_secondary_structure,
    calculate_boman_index,
    calculate_aliphatic_index,
    calculate_gravy_score,
    calculate_molecular_formula,
    calculate_isoelectric_point,
    count_amino_acids,
    format_subscripts,
    calculate_charge_distribution,
    format_molecular_formula,
    calculate_extinction_coefficient,
    calculate_epsilon_205,
    calculate_all_pis,
    N_TERM_MODIFICATIONS,
    C_TERM_MODIFICATIONS,
    HYDROPATHY_HOPP_WOODS,
    RESIDUE_COLORS,
    RESIDUE_CATEGORIES,
)

from plot_helpers import (
    generate_hydropathy_plot,
    generate_hopp_woods_plot,
    generate_net_charge_plot,
    generate_charge_distribution_pie,
)

from input_validator import validate_input

import uuid

app = Flask(__name__)

@app.template_filter('datetimeformat')
def datetimeformat(value):
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return dt.strftime("%B %d, %Y at %H:%M")
peptalyzer = Blueprint(
    'peptalyzer',
    __name__,
    url_prefix='/peptalyzer/app',
    static_folder='static',
    template_folder='templates'
)
CORS(app)

# Global in-memory storage for reports.
# WARNING: In a production environment with many users, this should be replaced by a database or Redis cache to prevent memory leaks.
CALCULATION_CACHE = {}

@peptalyzer.route("/")
def index():
    """Render the main app interface."""
    return render_template("index.html", current_year=datetime.now().year)

def process_peptide_calculation(data: dict) -> dict:
    """
    Core business logic to calculate peptide properties.
    Separated from the HTTP route for better structure and testability.

    Args:
        data (dict): A dictionary containing user inputs (sequence, pH, disulfide bonds, etc.).

    Returns:
        dict: A dictionary containing all calculated chemical properties, plot data, and metadata.
    """
    sequence = data.get("sequence", "").strip().upper()
    pH = float(data.get("pH", 7.0))
    pH_min = float(data.get("pH_min", 0.0))
    pH_max = float(data.get("pH_max", 14.0))
    pH_step = float(data.get("pH_step", 0.1))
    debug = bool(data.get("debug", False))

    # Get separate disulfide bond inputs
    intra_disulfide_bonds = int(data.get("intra_disulfide_bonds", 0))  # Intramolecular
    inter_disulfide_bonds = int(data.get("inter_disulfide_bonds", 0))  # Intermolecular

    total_disulfide_bonds = intra_disulfide_bonds + inter_disulfide_bonds
    peptide_units = inter_disulfide_bonds + 1  # Each intermolecular bond links two peptides

    if inter_disulfide_bonds == 0:
        peptide_state = "Monomer"
    else:
        peptide_state = f"Oligomer ({peptide_units} units)"

    n_term = data.get("n_term", "H")    # Default = free amine
    c_term = data.get("c_term", "OH")   # Default = carboxylic acid
    net_charge_at_pH = round(calculate_net_charge(sequence, pH, n_term, c_term, "IPC2_peptide"), 2)

    # Step 1: Calculate base monoisotopic mass for one peptide unit
    base_mono_mw = calculate_monoisotopic_mass(sequence, n_term, c_term)

    # Step 2: Calculate hydrogen loss PER MONOMER
    # Logic: 2 Hydrogens lost per intramolecular bond, 1 Hydrogen lost per intermolecular bond connection.
    hydrogen_loss_per_monomer = (2 * intra_disulfide_bonds) + (1 * inter_disulfide_bonds)

    # Step 3: Apply hydrogen loss to the monomer mass, then scale to oligomer
    # 1.0078 is the monoisotopic mass of a single Hydrogen atom.
    mono_mw = (base_mono_mw - (hydrogen_loss_per_monomer * 1.0078)) * peptide_units

    charge_distribution = calculate_charge_distribution(sequence)

    analysis = ProteinAnalysis(sequence)
    peptide_polarity = classify_peptide_polarity(sequence)

    # Step 4: Average mass correction
    # 1.00794 is the average atomic mass of a single Hydrogen atom.
    base_avg_mw = calculate_average_mass(sequence, n_term, c_term)
    avg_mw = (base_avg_mw - (hydrogen_loss_per_monomer * 1.00794)) * peptide_units

    # Prepare hydropathy data
    window_size = int(data.get("hydropathy_window", 9))
    y = smooth_hydropathy(sequence, window_size)
    y_hw = [HYDROPATHY_HOPP_WOODS.get(aa, 0) for aa in sequence]
    x = list(sequence)

    colors = [RESIDUE_COLORS.get(aa, '#333333') for aa in sequence]
    hydropathy_base64 = generate_hydropathy_plot(x, y, colors)

    # --- Hopp-Woods Hydropathy Plot ---
    colors_hw = [RESIDUE_COLORS.get(aa, '#333333') for aa in sequence]
    hopp_woods_base64 = generate_hopp_woods_plot(x, y_hw, colors_hw)

    # Net charge calculation
    ph_values = [round(pH_min + i * pH_step, 2) for i in range(int((pH_max - pH_min) / pH_step) + 1)]

    # Calculate the three specific pI results
    pi_results = calculate_all_pis(sequence, n_term, c_term)
    pI_ipc2 = pi_results.get("IPC2_peptide", 0.0)
    pI_bjellqvist = pi_results.get("Bjellqvist", 0.0)
    pI_emboss = pi_results.get("EMBOSS", 0.0)
    pI_lehninger = pi_results.get("Lehninger", 0.0)

    charges_dict = {
        'ipc2': [calculate_net_charge(sequence, ph, n_term, c_term, "IPC2_peptide") for ph in ph_values],
        'bjellqvist': [calculate_net_charge(sequence, ph, n_term, c_term, "Bjellqvist") for ph in ph_values],
        'emboss': [calculate_net_charge(sequence, ph, n_term, c_term, "EMBOSS") for ph in ph_values],
        'lehninger': [calculate_net_charge(sequence, ph, n_term, c_term, "Lehninger") for ph in ph_values]
    }

    netcharge_base64 = generate_net_charge_plot(ph_values, charges_dict, pI_ipc2)

    # Pie Chart
    acidic = charge_distribution['acidic_count']
    basic = charge_distribution['basic_count']
    charge_pie_base64 = generate_charge_distribution_pie(acidic, basic)

    # Extinction coefficient
    extinction_data = calculate_extinction_coefficient(sequence, total_disulfide_bonds, peptide_units)
    epsilon_205 = calculate_epsilon_205(sequence)

    aromaticity = analysis.aromaticity()

    # Calculate Neutral pH Aromaticity (His-Included)
    his_count = sequence.count('H')
    std_aromatic_count = sequence.count('F') + sequence.count('W') + sequence.count('Y')
    aromaticity_his = (std_aromatic_count + his_count) / len(sequence) if len(sequence) > 0 else 0.0

    instability_index = analysis.instability_index()
    gravy_result = calculate_gravy_score(sequence)
    secondary_structure = estimate_secondary_structure(sequence)
    boman = calculate_boman_index(sequence)
    aliphatic_index = calculate_aliphatic_index(sequence)

    # Molecular formula
    molecular_formula = calculate_molecular_formula(sequence, n_term, c_term)
    if 'H' in molecular_formula['raw']:
        molecular_formula['raw']['H'] -= hydrogen_loss_per_monomer
        if molecular_formula['raw']['H'] < 0:
            molecular_formula['raw']['H'] = 0

    if peptide_units > 1:
        for atom in molecular_formula['raw']:
            molecular_formula['raw'][atom] *= peptide_units

    molecular_formula['formatted'] = format_molecular_formula(molecular_formula['raw'])

    total_charge_residues = charge_distribution["acidic_count"] + charge_distribution["basic_count"]
    
    # Enrich amino acid counts with metadata for the frontend table
    raw_counts = count_amino_acids(sequence)
    amino_acid_counts = {}
    for aa, count in raw_counts.items():
        amino_acid_counts[aa] = {
            "count": count,
            "category": RESIDUE_CATEGORIES.get(aa, "Unknown"),
            "color": RESIDUE_COLORS.get(aa, "#333333")
        }

    n_term_full = format_subscripts(N_TERM_MODIFICATIONS.get(n_term, {}).get("label", n_term))
    c_term_full = format_subscripts(C_TERM_MODIFICATIONS.get(c_term, {}).get("label", c_term))

    return {
        "sequence": sequence,
        "length": len(sequence),
        "terminal_modifications": f"{n_term_full} / {c_term_full}",
        "monoisotopic_mw": round(mono_mw, 4),
        "average_mw": round(avg_mw, 4),
        "pI": pI_ipc2,
        "pI_ipc2": pI_ipc2,
        "pI_bjellqvist": pI_bjellqvist,
        "pI_emboss": pI_emboss,
        "pI_lehninger": pI_lehninger,
        "amino_acid_counts": amino_acid_counts,
        "isoelectric_point": pI_ipc2,
        "pi_scales": {
            "bjellqvist": pI_bjellqvist,
            "emboss": pI_emboss,
            "lehninger": pI_lehninger
        },
        "charge_at_pH": net_charge_at_pH,
        "selected_pH": pH,
        "n_term": n_term,
        "c_term": c_term,
        "pka_set": "standard (Biopython)",
        "hydropathy_plot": hydropathy_base64,
        "hopp_woods_plot": hopp_woods_base64,
        "net_charge_plot": netcharge_base64,
        "charge_distribution_pie": charge_pie_base64,
        "extinction_coefficient": {
                "adjusted": extinction_data["extinction_coefficient"],
                "unit": extinction_data["unit"]
            },
        "extinction_coefficient_205": {"value": epsilon_205, "unit": "M⁻¹·cm⁻¹"},
        "aromaticity": {
            "value": round(aromaticity, 4),
            "value_his": round(aromaticity_his, 4),
            "has_his": his_count > 0
        },
        "instability_index": {"value": round(instability_index, 2)},
        "gravy_value": gravy_result["value"],
        "intra_disulfide_bonds": intra_disulfide_bonds,
        "inter_disulfide_bonds": inter_disulfide_bonds,
        "peptide_state": peptide_state,
        "gravy_classification": gravy_result["annotation"],
        "hydropathy_window": window_size,
        "secondary_structure": secondary_structure,
        "boman_index": boman,
        "aliphatic_index": aliphatic_index,
        "charge_distribution": charge_distribution,
        "molecular_formula": molecular_formula,
        "total_charge_residues": total_charge_residues,
        "peptide_polarity": peptide_polarity,
        "sequence_valid": True,
        "current_year": datetime.now().year,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "debug_info": {
            "ph_values": ph_values,
            "net_charges": charges_dict['ipc2'],
            "raw_extinction_coefficient": extinction_data,
            "net_charge_at_input_pH": net_charge_at_pH,
            "peptide_polarity": peptide_polarity,
            "hydrogen_loss_per_monomer": hydrogen_loss_per_monomer,
            "peptide_units": peptide_units,
            "formula_raw": dict(molecular_formula["raw"])
        } if debug else None
    }

@peptalyzer.route("/calculate", methods=["POST"])
def calculate_properties():
    """Receive peptide input via JSON, calculate all relevant properties, and return them."""
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        return jsonify({"error": "Request body must be valid JSON with Content-Type: application/json"}), 400
    debug = bool(data.get("debug", False))
    pH_min = float(data.get("pH_min", 0.0))
    pH_max = float(data.get("pH_max", 14.0))
    pH_step = float(data.get("pH_step", 0.1))
    if pH_min >= pH_max or pH_step <= 0:
        return jsonify({"error": "Invalid pH range or step size."}), 400

    is_valid, message = validate_input(data)
    if not is_valid:
        return jsonify({"error": message}), 400

    pH = float(data.get("pH", 7.0))
    if not (0.0 <= pH <= 14.0):
        return jsonify({"error": "pH must be between 0 and 14"}), 400

    try:
        results = process_peptide_calculation(data)
        
        # Generate a unique ID for this calculation
        report_id = str(uuid.uuid4())
        CALCULATION_CACHE[report_id] = results
        
        # Return the ID to the frontend so it can request the specific report
        results["report_id"] = report_id
        
        return jsonify(results)

    except Exception as e:
        import traceback
        import logging
        app.logger.error(traceback.format_exc())
        return jsonify({"error": "An internal error occurred.", "details": str(e)}), 500

@peptalyzer.route("/report/preview")
def preview_report():
    """Render the most recent peptide results as a printable report."""
    report_id = request.args.get("id")
    
    if not report_id or report_id not in CALCULATION_CACHE:
        return "Report not found or expired. Please recalculate.", 404
        
    return render_template("report_template.html", data=CALCULATION_CACHE[report_id])

app.register_blueprint(peptalyzer)

if __name__ == "__main__":
    app.run(debug=True)