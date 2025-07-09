from flask import Flask, Blueprint, request, jsonify, render_template, send_file
from flask_cors import CORS
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from datetime import datetime, timezone
from weasyprint import HTML
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
    N_TERM_MODIFICATIONS,
    C_TERM_MODIFICATIONS,
    HYDROPATHY_HOPP_WOODS,
    RESIDUE_COLORS,
)

from plot_helpers import (
    generate_hydropathy_plot,
    generate_hopp_woods_plot,
    generate_net_charge_plot,
    generate_charge_distribution_pie,
)

from input_validator import validate_input

import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for macOS compatibility
import matplotlib.pyplot as plt
import io
import base64

app = Flask(__name__)
peptideiq = Blueprint(
    'peptideiq',
    __name__,
    url_prefix='/peptideiq',
    static_folder='static',
    template_folder='templates'
)
CORS(app)
latest_results = {}

@peptideiq.route("/")
def index():
    return render_template("index.html")

@peptideiq.route("/calculate", methods=["POST"])
def calculate_properties():
    data = request.get_json()
    debug = bool(data.get("debug", False))
    pH_min = float(data.get("pH_min", 0.0))
    pH_max = float(data.get("pH_max", 14.0))
    pH_step = float(data.get("pH_step", 0.1))
    if pH_min >= pH_max or pH_step <= 0:
        return jsonify({"error": "Invalid pH range or step size."}), 400

    is_valid, message = validate_input(data)
    if not is_valid:
        return jsonify({"error": message}), 400

    sequence = data.get("sequence", "").strip().upper()
    # Get separate disulfide bond inputs
    intra_disulfide_bonds = int(data.get("intra_disulfide_bonds", 0))  # Intramolecular
    inter_disulfide_bonds = int(data.get("inter_disulfide_bonds", 0))  # Intermolecular

    total_disulfide_bonds = intra_disulfide_bonds + inter_disulfide_bonds  # Total SS bonds
    peptide_units = inter_disulfide_bonds + 1  # Each intermolecular bond links two peptides
    # Calculate hydrogen loss per monomer (intra: 2H per bond, inter: 1H per bond)
    hydrogen_loss_per_monomer = (2 * intra_disulfide_bonds) + 1  # Always 1 hydrogen loss per monomer for intermolecular SS

    if inter_disulfide_bonds == 0:
        peptide_state = "Monomer"
    else:
        peptide_state = f"Oligomer ({peptide_units} units)"

    pH = float(data.get("pH", 7.0))
    n_term = data.get("n_term", "H")    # Default = free amine
    c_term = data.get("c_term", "OH")   # Default = carboxylic acid
    net_charge_at_pH = round(calculate_net_charge(sequence, pH, n_term, c_term), 2)

    # Step 1: Calculate base monoisotopic mass for one peptide unit
    base_mono_mw = calculate_monoisotopic_mass(sequence, n_term, c_term)

    # Step 2: Calculate hydrogen loss PER MONOMER (before multiplication)
    hydrogen_loss_per_monomer = (2 * intra_disulfide_bonds) + (1 * inter_disulfide_bonds)

    # Step 3: Apply hydrogen loss to the monomer mass, then scale to oligomer
    mono_mw = (base_mono_mw - (hydrogen_loss_per_monomer * 1.0078)) * peptide_units

    charge_distribution = calculate_charge_distribution(sequence)

    if not (0.0 <= pH <= 14.0):
        return jsonify({"error": "pH must be between 0 and 14"}), 400

    try:
        analysis = ProteinAnalysis(sequence)

        peptide_polarity = classify_peptide_polarity(sequence)

        # Average MW from Biopython
        # Calculate average mass and adjust for disulfide bonds (-2.0156 Da per bond)
        # Step 4: Average mass correction using the same logic
        base_avg_mw = calculate_average_mass(sequence, n_term, c_term)
        avg_mw = (base_avg_mw - (hydrogen_loss_per_monomer * 1.00794)) * peptide_units

        # Prepare hydropathy data
        window_size = int(data.get("hydropathy_window", 5))  # default to 5
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
        net_charges = [round(calculate_net_charge(sequence, ph, n_term, c_term), 2) for ph in ph_values]
        pI = calculate_isoelectric_point(sequence, n_term, c_term)

        # Net charge plot
        netcharge_base64 = generate_net_charge_plot(ph_values, net_charges, pI)

        # Pie Chart or No Data Box
        acidic = charge_distribution['acidic_count']
        basic = charge_distribution['basic_count']

        charge_pie_base64 = generate_charge_distribution_pie(acidic, basic)

        pI_value = pI
        # Calculate monomer extinction coefficient (excluding cysteines)
        extinction_data = calculate_extinction_coefficient(sequence, total_disulfide_bonds, peptide_units)

        # Correct: extinction coefficient includes base monomer + disulfide bonds per oligomer
        aromaticity = analysis.aromaticity()
        instability_index = analysis.instability_index()
        gravy_result = calculate_gravy_score(sequence)
        gravy_value = gravy_result["value"]
        gravy_classification = gravy_result["annotation"]
        secondary_structure = estimate_secondary_structure(sequence)
        boman = calculate_boman_index(sequence)
        aliphatic_index = calculate_aliphatic_index(sequence)
        # Step 5: Calculate molecular formula for one peptide unit
        molecular_formula = calculate_molecular_formula(sequence, n_term, c_term)

        # Step 7: Remove hydrogens from the monomer BEFORE multiplication
        if 'H' in molecular_formula['raw']:
            molecular_formula['raw']['H'] -= hydrogen_loss_per_monomer
            if molecular_formula['raw']['H'] < 0:
                molecular_formula['raw']['H'] = 0  # Prevent negative hydrogen count

        # Step 8: Now scale the ENTIRE formula to the number of peptide units
        if peptide_units > 1:
            for atom in molecular_formula['raw']:
                molecular_formula['raw'][atom] *= peptide_units

        # Step 9: Reformat molecular formula for display
        molecular_formula['formatted'] = format_molecular_formula(molecular_formula['raw'])

        # Reformat the molecular formula for display after adjustments
        molecular_formula['formatted'] = format_molecular_formula(molecular_formula['raw'])

        total_charge_residues = charge_distribution["acidic_count"] + charge_distribution["basic_count"]

        debug_info = {}

        if debug:
            debug_info = {
                "ph_values": ph_values,
                "net_charges": net_charges,
                "raw_extinction_coefficient": extinction_data,
                "net_charge_at_input_pH": net_charge_at_pH,
                "peptide_polarity": peptide_polarity
            }
            debug_info["hydrogen_loss_per_monomer"] = hydrogen_loss_per_monomer
            debug_info["peptide_units"] = peptide_units
            debug_info["formula_raw"] = dict(molecular_formula["raw"])

        n_term_full_raw = N_TERM_MODIFICATIONS.get(n_term, {}).get("label", n_term)
        c_term_full_raw = C_TERM_MODIFICATIONS.get(c_term, {}).get("label", c_term)

        n_term_full = format_subscripts(n_term_full_raw)
        c_term_full = format_subscripts(c_term_full_raw)

        amino_acid_counts = count_amino_acids(sequence)

        results = {
            "sequence": sequence,
            "length": len(sequence),
            "terminal_modifications": f"{n_term_full} / {c_term_full}",
            "monoisotopic_mw": round(mono_mw, 4),
            "average_mw": round(avg_mw, 4),
            "amino_acid_counts": amino_acid_counts,
            "isoelectric_point": pI,
            "charge_at_pH": net_charge_at_pH,
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
            "aromaticity": {"value": round(aromaticity, 4)},
            "instability_index": {"value": round(instability_index, 2)},
            "gravy_value": gravy_value,
            "intra_disulfide_bonds": intra_disulfide_bonds,
            "inter_disulfide_bonds": inter_disulfide_bonds,
            "peptide_state": peptide_state,
            "gravy_classification": gravy_classification,
            "secondary_structure": secondary_structure,
            "boman_index": boman,
            "aliphatic_index": aliphatic_index,
            "charge_distribution": charge_distribution,
            "molecular_formula": molecular_formula,
            "total_charge_residues": total_charge_residues,
            "peptide_polarity": peptide_polarity,
            "sequence_valid": True,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "debug_info": debug_info if debug else None
        }
        
        global latest_results
        latest_results = results
        return jsonify(results)

    except Exception as e:
        import traceback
        import logging
        logging.basicConfig(level=logging.ERROR)
        app.logger.error(traceback.format_exc())
        return jsonify({"error": "An internal error occurred.", "details": str(e)}), 500

@peptideiq.route("/export_pdf")
def export_pdf():
    global latest_results

    if not latest_results:
        return "No calculation data available to export.", 400

    # Get the base URL for static assets
    base_url = request.host_url  # Example: 'http://127.0.0.1:5000/'

    # Render the HTML with absolute URLs
    rendered_html = render_template("report_template.html", data=latest_results)

    try:
        # Generate PDF using WeasyPrint
        pdf_io = io.BytesIO()
        HTML(string=rendered_html, base_url=request.host_url).write_pdf(pdf_io)
        pdf_io.seek(0)

        return send_file(
            pdf_io,
            as_attachment=True,
            download_name="peptide_report.pdf",
            mimetype='application/pdf'
        )

    except Exception as e:
        return f"Error generating PDF: {e}", 500

app.register_blueprint(peptideiq)

from flask import redirect

if __name__ == "__main__":
    app.run(debug=True)