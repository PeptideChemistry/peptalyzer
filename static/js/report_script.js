document.addEventListener("DOMContentLoaded", () => {
  // ====== DOM Elements ======
  const exportCsvBtn = document.getElementById("exportCsvBtn");
  const exportPdfBtn = document.getElementById("exportPdfBtn");

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
      placeholder.textContent = "Graph will appear here"; // Reset to default when graph is valid
    } else {
      imgElement.style.display = "none";
      placeholder.style.display = "flex";
      placeholder.textContent = customMessage; // Show the specific message
    }
  }

  // ====== Update Results Table ======
  function updateResults(data) {
    if (data.length !== undefined) {
      document.getElementById("result-length").textContent = data.length;
    }

    if (data.intra_disulfide_bonds !== undefined) {
      document.getElementById("result-intra-disulfide-bonds").textContent = data.intra_disulfide_bonds;
    }

    if (data.inter_disulfide_bonds !== undefined) {
      document.getElementById("result-inter-disulfide-bonds").textContent = data.inter_disulfide_bonds;
    }

    if (data.peptide_state !== undefined) {
      document.getElementById("result-peptide-state").textContent = data.peptide_state;
    }

    if (data.terminal_modifications !== undefined) {
      document.getElementById("result-termini").innerHTML = data.terminal_modifications;
    }

    if (data.isoelectric_point !== undefined) {
      document.getElementById("result-pi").textContent = data.isoelectric_point.toFixed(2);
    }

    if (data.charge_at_pH !== undefined) {
      document.getElementById("result-charge").textContent = data.charge_at_pH.toFixed(2);
    }

    if (data.gravy_value !== undefined) {
      document.getElementById("result-gravy").textContent = data.gravy_value.toFixed(2);
      document.getElementById("result-gravy-classification").textContent = data.gravy_classification;
    }

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

    // Pie Chart
    if (data.charge_distribution_pie === "NO_DATA") {
      updateGraph(
        "img[alt='Pie chart showing distribution of acidic and basic residues in the peptide']",
        null,
        "No charged residues in the sequence. Pie chart not applicable."
      );
    } else {
      updateGraph(
        "img[alt='Pie chart showing distribution of acidic and basic residues in the peptide']",
        data.charge_distribution_pie
      );
    }

    // Net Charge vs pH Graph
    updateGraph("img[alt='Graph showing net charge of the peptide as a function of pH']", data.net_charge_plot);

    // Kyte-Doolittle Graph
    updateGraph("img[alt='Bar graph showing Kyte-Doolittle hydropathy profile of the peptide']", data.hydropathy_plot);

    // Hopp-Woods Graph
    updateGraph("img[alt='Bar graph showing Hopp-Woods hydrophilicity profile of the peptide']", data.hopp_woods_plot);

    // Amino Acid Counts Table
    updateAminoAcidCountsTable(data.amino_acid_counts);

    // Update timestamp
    const now = new Date();
    const formatted = now.toLocaleString();
    document.getElementById('timestamp').textContent = formatted;
  }

  // ====== Update Amino Acid Counts Table ======
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

  // ====== CSV Export Handler ======
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

  // ====== PDF Export Handler ======
  exportPdfBtn.addEventListener("click", () => {
    window.location.href = "/export_pdf"; // Calls your Flask backend to generate the PDF
  });

});