document.addEventListener("DOMContentLoaded", () => {
  // ====== DOM Elements ======
  const form = document.getElementById("calcForm");
  const phInput = document.getElementById("phInput");
  const phValueDisplay = document.getElementById("phValue");
  const calcButton = document.getElementById("calcButton");
  const exportCsvBtn = document.getElementById("exportCsvBtn");

  // ====== Disable Export Buttons ======
  function disableExportButtons() {
    exportCsvBtn.disabled = true;
    exportPdfBtn.disabled = true;
  }

  // ====== New DOM Elements for Disulfide Bonds ======
  const intraInput = document.getElementById("intra_disulfide_bonds");
  const interInput = document.getElementById("inter_disulfide_bonds");
  const sequenceInput = document.getElementById("sequence");

  // ====== Real-Time pH Slider Display ======
  phInput.addEventListener("input", () => {
    phValueDisplay.textContent = phInput.value;
    document.getElementById("displayed-pH").textContent = phInput.value;
    phInput.setAttribute("aria-valuetext", phInput.value);
    disableExportButtons();
  });

  // ====== Real-Time Cysteine Validation ======
  sequenceInput.addEventListener("input", () => {
    disableExportButtons();

    const sequence = sequenceInput.value.toUpperCase();
    const cysCount = (sequence.match(/C/g) || []).length;

    if (cysCount === 0) {
      intraInput.value = 0;
      interInput.value = 0;
      intraInput.disabled = true;
      interInput.disabled = true;
    } else {
      intraInput.disabled = false;
      interInput.disabled = false;

      const maxIntraBonds = Math.floor(cysCount / 2);
      if (parseInt(intraInput.value, 10) > maxIntraBonds) {
        intraInput.value = maxIntraBonds;
      }
      intraInput.max = maxIntraBonds;

      const remainingCys = cysCount - (parseInt(intraInput.value, 10) * 2);
      const maxInterBonds = remainingCys;

      if (parseInt(interInput.value, 10) > maxInterBonds) {
        interInput.value = maxInterBonds;
      }
      interInput.max = maxInterBonds;
    }
  });

// ====== Additional Export Disable Triggers ======
document.getElementById("n_term").addEventListener("change", disableExportButtons);
document.getElementById("c_term").addEventListener("change", disableExportButtons);
intraInput.addEventListener("input", disableExportButtons);
interInput.addEventListener("input", disableExportButtons);

  // ====== Graph Reset Utility ======
  function resetGraphs() {
    const graphs = [
      { alt: "Graph showing net charge of the peptide as a function of pH" },
      { alt: "Bar graph showing Kyte-Doolittle hydropathy profile of the peptide" },
      { alt: "Bar graph showing Hopp-Woods hydrophilicity profile of the peptide" },
      { alt: "Pie chart showing distribution of acidic and basic residues in the peptide" }
    ];

    graphs.forEach(graph => {
      const img = document.querySelector(`img[alt='${graph.alt}']`);
      const placeholder = img.nextElementSibling;
      img.style.display = "none";
      img.src = "";
      placeholder.style.display = "flex";
    });
  }

  // ====== Update Graph Helper ======
  function updateGraph(imgSelector, plotData, customMessage = "Graph will appear here") {
    const imgElement = document.querySelector(imgSelector);
    const placeholder = imgElement.nextElementSibling;

    if (plotData) {
      imgElement.style.display = "block";
      imgElement.src = "data:image/png;base64," + plotData;
      placeholder.style.display = "none";
      placeholder.textContent = "Graph will appear here";
    } else {
      imgElement.style.display = "none";
      placeholder.style.display = "flex";
      placeholder.textContent = customMessage;
    }
  }

  // ====== Update Results Table ======
  function updateResults(data) {
    document.getElementById("result-length").textContent = data.length;
    document.getElementById("result-intra-disulfide-bonds").textContent = data.intra_disulfide_bonds;
    document.getElementById("result-inter-disulfide-bonds").textContent = data.inter_disulfide_bonds;
    document.getElementById("result-peptide-state").textContent = data.peptide_state;
    document.getElementById("result-termini").innerHTML = data.terminal_modifications;
    document.getElementById("result-pi").textContent = data.isoelectric_point.toFixed(2);
    document.getElementById("result-charge").textContent = data.charge_at_pH.toFixed(2);
    document.getElementById("result-gravy").textContent = data.gravy_value.toFixed(2);
    document.getElementById("result-gravy-classification").textContent = data.gravy_classification;
    document.getElementById("result-polarity").textContent = data.peptide_polarity;
    document.getElementById("result-mono-mw").textContent = data.monoisotopic_mw.toFixed(4) + " Da";
    document.getElementById("result-avg-mw").textContent = data.average_mw.toFixed(4) + " Da";
    document.getElementById("result-boman").textContent = data.boman_index;
    document.getElementById("result-aliphatic").textContent = data.aliphatic_index.value + "% (" + data.aliphatic_index.rating + ")";
    document.getElementById("result-aromaticity").textContent = data.aromaticity.value;
    document.getElementById("result-instability").textContent = data.instability_index.value;
    document.getElementById("result-extinction").textContent = data.extinction_coefficient.adjusted + ' ' + data.extinction_coefficient.unit;

    const formulaWithSubscripts = data.molecular_formula.formatted.replace(/(\d+)/g, "<sub>$1</sub>");
    document.getElementById("result-formula").innerHTML = formulaWithSubscripts;

    document.getElementById("result-helix").textContent = data.secondary_structure.helix;
    document.getElementById("result-sheet").textContent = data.secondary_structure.sheet;
    document.getElementById("result-coil").textContent = data.secondary_structure.coil;

    document.getElementById("result-acidic-count").textContent = data.charge_distribution.acidic_count;
    document.getElementById("result-basic-count").textContent = data.charge_distribution.basic_count;
    document.getElementById("result-total-charged").textContent = data.total_charge_residues;

    updateGraph("img[alt='Pie chart showing distribution of acidic and basic residues in the peptide']", data.charge_distribution_pie === "NO_DATA" ? null : data.charge_distribution_pie);
    updateGraph("img[alt='Graph showing net charge of the peptide as a function of pH']", data.net_charge_plot);
    updateGraph("img[alt='Bar graph showing Kyte-Doolittle hydropathy profile of the peptide']", data.hydropathy_plot);
    updateGraph("img[alt='Bar graph showing Hopp-Woods hydrophilicity profile of the peptide']", data.hopp_woods_plot);

    updateAminoAcidCountsTable(data.amino_acid_counts);

    const now = new Date();
    document.getElementById('timestamp').textContent = now.toLocaleString();
  }

  function updateAminoAcidCountsTable(counts) {
    const tableBody = document.querySelector("#aaCountsTable tbody");
    tableBody.innerHTML = "";

    if (counts && Object.keys(counts).length > 0) {
      for (const [aminoAcid, count] of Object.entries(counts)) {
        const row = document.createElement("tr");
        const aaCell = document.createElement("td");
        aaCell.textContent = aminoAcid;

        const countCell = document.createElement("td");
        countCell.textContent = count;

        row.appendChild(aaCell);
        row.appendChild(countCell);
        tableBody.appendChild(row);
      }
    } else {
      const row = document.createElement("tr");
      const cell = document.createElement("td");
      cell.colSpan = 2;
      cell.textContent = "No amino acids found.";
      row.appendChild(cell);
      tableBody.appendChild(row);
    }
  }

  exportCsvBtn.addEventListener("click", () => {
    const rows = document.querySelectorAll("#resultsTable tr");
    let csvContent = "";

    rows.forEach(row => {
      const cols = row.querySelectorAll("th, td");
      const rowData = Array.from(cols).map(col => col.innerText).join(",");
      csvContent += rowData + "\r\n";
    });

    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);

    const link = document.createElement("a");
    link.setAttribute("href", url);
    link.setAttribute("download", "peptide_results.csv");
    link.style.display = "none";
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  });

  exportPdfBtn.addEventListener("click", () => {
    window.location.href = "/peptideiq/export_pdf";
  });

  form.addEventListener("submit", async (e) => {
    e.preventDefault();

    calcButton.disabled = true;
    calcButton.textContent = "Calculating...";

    resetGraphs();
    disableExportButtons();

    const sequence = document.getElementById("sequence").value;
    const n_term = document.getElementById("n_term").value;
    const c_term = document.getElementById("c_term").value;
    const selected_pH = parseFloat(phInput.value);
    const intra_disulfide_bonds = parseInt(document.getElementById("intra_disulfide_bonds").value, 10);
    const inter_disulfide_bonds = parseInt(document.getElementById("inter_disulfide_bonds").value, 10);

    try {
      const response = await fetch("/peptideiq/calculate", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence,
          n_term,
          c_term,
          pH_min: 0.0,
          pH_max: 14.0,
          pH_step: 0.1,
          pH: selected_pH,
          intra_disulfide_bonds,
          inter_disulfide_bonds
        })
      });

      if (!response.ok) {
        const errorData = await response.json();
        alert(`Error: ${errorData.error}`);
        return;
      }

      const data = await response.json();
      updateResults(data);

      // Enable export buttons
      exportCsvBtn.disabled = false;
      exportPdfBtn.disabled = false;

    } catch (error) {
      console.error("Network error:", error);
      alert("A network error occurred. Please try again.");
    } finally {
      calcButton.disabled = false;
      calcButton.textContent = "Calculate";
    }
  });
});
