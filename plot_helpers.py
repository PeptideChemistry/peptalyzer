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
# 📈 Multi-Scale Net Charge vs pH Plot
def generate_net_charge_plot(ph_values, charges_dict, pI_ipc2):
    """
    Generates a multi-line Net Charge vs pH plot using brand colors.
    charges_dict should look like: {'ipc2': [...], 'bjellqvist': [...], 'emboss': [...]}
    """
    plt.figure(figsize=(7.5, 4), facecolor='none')
    plt.gca().patch.set_alpha(0)

    # 1. Plot each scale with brand colors
    # Primary: IPC 2.0 (Red)
    plt.plot(ph_values, charges_dict['ipc2'], color='#EB5247', linewidth=2.5, label='IPC 2.0')
    
    # Secondary: Bjellqvist (Blue)
    plt.plot(ph_values, charges_dict['bjellqvist'], color='#1E90FF', linewidth=1.5, label='Bjellqvist', alpha=0.8)
    
    # Accent/Element: EMBOSS (Yellow)
    plt.plot(ph_values, charges_dict['emboss'], color='#F1C40F', linewidth=1.5, label='EMBOSS', alpha=0.8)

    # Accent/Element: Lehninger (Green)
    plt.plot(ph_values, charges_dict['lehninger'], color='#2ECC71', linewidth=1.5, label='Lehninger', alpha=0.8)

    # 2. Add Reference Lines
    plt.axhline(0, color='#333333', linewidth=1.0) # Zero line
    plt.axvline(pI_ipc2, color='#EB5247', linestyle='--', linewidth=0.8, alpha=0.6) # pI marker
    
    # 3. Text and Labels
    plt.text(pI_ipc2 + 0.2, 0.5, f"pI (IPC2) ≈ {pI_ipc2}", fontsize=9, color='#EB5247', fontweight='bold')
    plt.xlabel('pH', fontsize=11, color='#333333')
    plt.ylabel('Net Charge', fontsize=11, color='#333333')
    
    # 4. Legend styling
    plt.legend(loc='upper right', fontsize=9, frameon=True, facecolor='#F9F9F9', edgecolor='#333333')
    
    plt.grid(True, linestyle=':', linewidth=0.5, color='#cccccc')

    # Spines/Border
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
