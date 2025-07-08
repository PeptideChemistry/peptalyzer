import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for macOS compatibility
import matplotlib.pyplot as plt
import io
import base64

# 🖼️ Utility to save plot to base64
def save_plot_to_base64():
    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png', dpi=300, transparent=True)  # Increased DPI for crispness
    plt.close()
    img_bytes.seek(0)
    return base64.b64encode(img_bytes.read()).decode('utf-8')

# 📊 Hydropathy Plot (Kyte-Doolittle)
def generate_hydropathy_plot(x, y, colors, ylabel="Hydropathy Index (Kyte-Doolittle)"):
    """
    Generates a Kyte-Doolittle hydropathy plot and returns base64 string.
    """
    plt.figure(figsize=(7.5, 3), facecolor='none')
    plt.gca().patch.set_alpha(0)
    plt.bar(range(len(x)), y, color=colors, edgecolor='#333333', linewidth=0.8)
    plt.xticks(range(len(x)), x, fontsize=10)
    plt.yticks(fontsize=10)
    plt.axhline(0, color='#333333', linewidth=0.8)
    plt.ylabel(ylabel, fontsize=11)
    plt.grid(True, linestyle=':', linewidth=0.5, color='#cccccc')  # Softer gridlines

    for spine in plt.gca().spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('#333333')

    plt.tight_layout()
    return save_plot_to_base64()

# 📊 Hydropathy Plot (Hopp-Woods)
def generate_hopp_woods_plot(x, y_hw, colors, ylabel="Hydropathy Index (Hopp-Woods)"):
    """
    Generates a Hopp-Woods hydropathy plot and returns base64 string.
    """
    plt.figure(figsize=(7.5, 3), facecolor='none')
    plt.gca().patch.set_alpha(0)
    plt.bar(range(len(x)), y_hw, color=colors, edgecolor='#333333', linewidth=0.8)
    plt.xticks(range(len(x)), x, fontsize=10)
    plt.yticks(fontsize=10)
    plt.axhline(0, color='#333333', linewidth=0.8)
    plt.ylabel(ylabel, fontsize=11)
    plt.grid(True, linestyle=':', linewidth=0.5, color='#cccccc')  # Softer gridlines

    for spine in plt.gca().spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('#333333')

    plt.tight_layout()
    return save_plot_to_base64()

# 📈 Net Charge vs pH Plot
def generate_net_charge_plot(ph_values, net_charges, pI):
    """
    Generates a Net Charge vs pH plot and returns base64 string.
    """
    plt.figure(figsize=(7.5, 4), facecolor='none')  # Slightly taller plot for balance
    plt.gca().patch.set_alpha(0)
    plt.plot(ph_values, net_charges, color='#1e90ff', linewidth=2)  # Updated to site blue
    plt.axhline(0, color='#333333', linewidth=0.8)
    plt.axvline(pI, color='#333333', linestyle='--', linewidth=0.8)
    plt.text(pI + 0.2, 0.5, f"pI ≈ {pI}", fontsize=9, color='#333333')
    plt.xlabel('pH', fontsize=11)
    plt.ylabel('Net Charge', fontsize=11)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(True, linestyle=':', linewidth=0.5, color='#cccccc')  # Softer gridlines

    for spine in plt.gca().spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('#333333')

    plt.tight_layout()
    return save_plot_to_base64()

# 🥧 Charge Distribution Pie Chart (Flat, No Branding)
def generate_charge_distribution_pie(acidic_count, basic_count):
    import matplotlib.pyplot as plt
    import numpy as np

    if acidic_count == 0 and basic_count == 0:
        return "NO_DATA"

    labels = ['Acidic (D, E)', 'Basic (R, K, H)']
    sizes = [acidic_count, basic_count]
    colors = ['#ff7f7f', '#1e90ff']  # Soft coral and site blue

    plt.figure(figsize=(5, 5), facecolor='none')
    plt.gca().patch.set_alpha(0)

    wedges, texts, autotexts = plt.pie(
        sizes,
        colors=colors,
        startangle=90,
        autopct='%1.0f%%',
        textprops={'color': 'white', 'weight': 'bold', 'fontsize': 12},
        wedgeprops={'linewidth': 1.5, 'edgecolor': '#333333'}
    )

    plt.legend(
        wedges,
        labels,
        title="Charge Type",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=2,
        frameon=False,
        fontsize=13,
        title_fontsize=14
    )

    for spine in plt.gca().spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('#333333')

    plt.tight_layout()
    return save_plot_to_base64()
