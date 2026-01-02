document.addEventListener("DOMContentLoaded", () => {
  // ====== DOM Elements ======
  const form = document.getElementById("calcForm");
  const phInput = document.getElementById("phInput");
  const phValueDisplay = document.getElementById("phValue");
  const calcButton = document.getElementById("calcButton");
  const exportCsvBtn = document.getElementById("exportCsvBtn");
  const exportPdfBtn = document.getElementById("exportPdfBtn");

  // ====== Disable Export Buttons ======
  function disableExportButtons() {
    exportCsvBtn.disabled = true;
    exportPdfBtn.disabled = true;
  }

  // ====== Disulfide Bond Constraint Logic ======
  function updateDisulfideConstraints() {
    const sequence = sequenceInput.value.toUpperCase();
    const cysCount = (sequence.match(/C/g) || []).length;
    const currentIntra = parseInt(intraInput.value, 10) || 0;

    if (cysCount === 0) {
      intraInput.value = 0;
      interInput.value = 0;
      intraInput.disabled = true;
      interInput.disabled = true;
    } else {
      intraInput.disabled = false;
      interInput.disabled = false;

      // 1. Calculate max Intra bonds allowed (2 Cys per bond)
      const maxIntraBonds = Math.floor(cysCount / 2);
      
      // Clamp Intra value if it exceeds new max
      if (currentIntra > maxIntraBonds) {
        intraInput.value = maxIntraBonds;
      }
      intraInput.max = maxIntraBonds;

      // 2. Calculate remaining Cys for Inter bonds
      // (Use the potentially clamped value from above)
      const usedCysForIntra = (parseInt(intraInput.value, 10) || 0) * 2;
      const remainingCys = cysCount - usedCysForIntra;
      const maxInterBonds = remainingCys;

      if (parseInt(interInput.value, 10) > maxInterBonds) {
        interInput.value = maxInterBonds;
      }
      interInput.max = maxInterBonds;
    }
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
    updateDisulfideConstraints();
  });

// ====== Additional Export Disable Triggers ======
document.getElementById("n_term").addEventListener("change", disableExportButtons);
document.getElementById("c_term").addEventListener("change", disableExportButtons);
intraInput.addEventListener("input", () => {
  disableExportButtons();
  updateDisulfideConstraints();
});
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
    document.getElementById("result-pi").textContent = data.isoelectric_point.toFixed(3);
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
    document.getElementById('res-pI-ipc2').textContent = data.pI_ipc2.toFixed(3);
    document.getElementById('res-pI-bjellqvist').textContent = data.pI_bjellqvist.toFixed(3);
    document.getElementById('res-pI-emboss').textContent = data.pI_emboss.toFixed(3);
    document.getElementById('res-pI-lehninger').textContent = data.pI_lehninger.toFixed(3);

    const formulaWithSubscripts = data.molecular_formula.formatted.replace(/(\d+)/g, "<sub>$1</sub>");
    document.getElementById("result-formula").innerHTML = formulaWithSubscripts;

    document.getElementById("result-helix").textContent = data.secondary_structure.helix;
    document.getElementById("result-sheet").textContent = data.secondary_structure.sheet;
    document.getElementById("result-coil").textContent = data.secondary_structure.coil;

    document.getElementById("result-acidic-count").textContent = data.charge_distribution.acidic_count;
    document.getElementById("result-basic-count").textContent = data.charge_distribution.basic_count;
    document.getElementById("result-total-charged").textContent = data.total_charge_residues;

    const cdNoData = data.charge_distribution_pie === "NO_DATA";
    updateGraph(
      "img[alt='Pie chart showing distribution of acidic and basic residues in the peptide']",
      cdNoData ? null : data.charge_distribution_pie,
      cdNoData ? "No charged residues detected" : "Graph will appear here"
    );

    updateGraph("img[alt='Graph showing net charge of the peptide as a function of pH']", data.net_charge_plot);
    // KD: show ProtScale-style notice if sequence is too short for the KD window
    const usedWindow = data.hydropathy_window || 9;
    const minLen = 2 * usedWindow;
    const tooShortKD = data.length < minLen;

    const kdNotice =
      `Sorry. Your sequence should be at least ${minLen} residues long (twice the window size of ${usedWindow}). ` +
      `Hydropathy profiles are unreliable below this length. For shorter peptides, use the GRAVY ` +
      `(Grand Average of Hydropathicity) value instead.`;

    updateGraph(
      "img[alt='Bar graph showing Kyte-Doolittle hydropathy profile of the peptide']",
      tooShortKD ? null : data.hydropathy_plot,
      tooShortKD ? kdNotice : "Graph will appear here"
    );

    // HW: unchanged (no smoothing/window used in app.py)
    updateGraph("img[alt='Bar graph showing Hopp-Woods hydrophilicity profile of the peptide']", data.hopp_woods_plot);

    updateAminoAcidCountsTable(data.amino_acid_counts);

    const now = new Date();
    document.getElementById('timestamp').textContent = now.toLocaleString();
  }

  function updateAminoAcidCountsTable(counts) {
    const table = document.getElementById("aaCountsTable");
    const thead = table.querySelector("thead");
    const tbody = table.querySelector("tbody");
    tbody.innerHTML = "";

    if (!counts || Object.keys(counts).length === 0) {
      const row = document.createElement("tr");
      const cell = document.createElement("td");
      cell.colSpan = 3; // now 3 columns
      cell.textContent = "No amino acids found.";
      row.appendChild(cell);
      tbody.appendChild(row);
      return;
    }

    // total length for percent calculation (one decimal)
    const total = Object.values(counts).reduce((a, b) => a + b, 0);

    // default alphabetical order (A→Z)
    const entries = Object.entries(counts).sort((a, b) => a[0].localeCompare(b[0]));

    for (const [aa, count] of entries) {
      const pct = ((count / total) * 100).toFixed(1);

      const tr = document.createElement("tr");
      // store raw values for sorting
      tr.innerHTML = `
        <td data-type="str">${aa}</td>
        <td data-type="num" data-value="${count}">${count}</td>
        <td data-type="num" data-value="${pct}">${pct}</td>
      `;
      tbody.appendChild(tr);
    }
  }
function attachAATableSortHandlers(table) {
  const thead = table.querySelector("thead");
  if (!thead) return;
  const ths = Array.from(thead.querySelectorAll("th"));

  ths.forEach((th, index) => {
    // make headers clickable for sorting
    th.style.cursor = "pointer";
    th.setAttribute("title", "Click to sort");
    th.addEventListener("click", () => {
      // toggle sort dir
      const dir = th.dataset.sortDir === "asc" ? "desc" : "asc";
      ths.forEach(h => delete h.dataset.sortDir); // clear others
      th.dataset.sortDir = dir;

      const tbody = table.querySelector("tbody");
      const rows = Array.from(tbody.querySelectorAll("tr"));
      const type = index === 0 ? "str" : "num"; // AA column is text; others numeric

      rows.sort((r1, r2) => {
        const c1 = r1.children[index];
        const c2 = r2.children[index];
        const v1 = type === "num" ? parseFloat(c1.dataset.value || c1.textContent) : c1.textContent;
        const v2 = type === "num" ? parseFloat(c2.dataset.value || c2.textContent) : c2.textContent;
        if (type === "num") {
          return dir === "asc" ? v1 - v2 : v2 - v1;
        } else {
          return dir === "asc" ? v1.localeCompare(v2) : v2.localeCompare(v1);
        }
      });

      // reattach in new order
      rows.forEach(r => tbody.appendChild(r));
    }, { passive: true });
  });
}

  // Attach sort handlers once on load (since headers are static in HTML)
  const aaTable = document.getElementById("aaCountsTable");
  if (aaTable) attachAATableSortHandlers(aaTable);

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
    const reportId = exportPdfBtn.dataset.reportId;
    if (reportId) {
      window.open(`/peptalyzer/app/report/preview?id=${reportId}`, "_blank");
    }
  });

  form.addEventListener("submit", async (e) => {
    e.preventDefault();

    calcButton.disabled = true;
    calcButton.textContent = "Calculating...";

    const piIds = ['result-pi', 'res-pI-bjellqvist', 'res-pI-emboss', 'res-pI-lehninger'];
    piIds.forEach(id => {
      const el = document.getElementById(id);
      if (el) el.textContent = "..."; 
    });

    resetGraphs();
    disableExportButtons();

    const sequence = document.getElementById("sequence").value;
    const n_term = document.getElementById("n_term").value;
    const c_term = document.getElementById("c_term").value;
    const selected_pH = parseFloat(phInput.value);
    const intra_disulfide_bonds = parseInt(document.getElementById("intra_disulfide_bonds").value, 10);
    const inter_disulfide_bonds = parseInt(document.getElementById("inter_disulfide_bonds").value, 10);

    try {
      const response = await fetch("/peptalyzer/app/calculate", {
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
        let msg = "Server error.";
        try {
          const ct = response.headers.get("content-type") || "";
          if (ct.includes("application/json")) {
            const err = await response.json();
            msg = err.details ? `${err.error}: ${err.details}` : (err.error || msg);
          } else {
            const text = await response.text();
            msg = text || msg;
          }
        } catch (_) { /* ignore parse errors */ }
        alert(`Error: ${msg}`);
        return;
      }

      const data = await response.json();

      // Store the report ID on the button for the click handler to use
      if (data.report_id) exportPdfBtn.dataset.reportId = data.report_id;

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