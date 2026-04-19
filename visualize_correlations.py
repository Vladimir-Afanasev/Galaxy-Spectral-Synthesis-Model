import matplotlib.pyplot as plt
import re
try:
    with open('galaxy_analysis_report.txt', 'r', encoding='utf-8') as f:
        content = f.read()
except FileNotFoundError:
    print("Ошибка: Файл galaxy_analysis_report.txt не найден!")
    content = ""

def parse_report(text):
    galaxies = []
    blocks = re.split(r'-{30,}', text)
    for block in blocks:
        block = block.strip()
        if not block: continue
        
        name_match = re.search(r'^(.*?)\s\(', block, re.MULTILINE)
        if not name_match: continue
        name = name_match.group(1).strip()
        
        # Разделяем на секции, чтобы не перепутать B-V модели и реальности
        sections = re.split(r'РЕАЛЬНОСТЬ:', block)
        if len(sections) < 2: continue
        model_sec, real_sec = sections[0], sections[1]
        
        m_mod = re.search(r'M_bol:\s*([-+]?\d*\.?\d+)', model_sec)
        bv_mod = re.search(r'B-V:\s*([-+]?\d*\.?\d+)', model_sec)
        m_real = re.search(r'M_bol:\s*([-+]?\d*\.?\d+)(?:\.\.([-+]?\d*\.?\d+))?', real_sec)
        bv_real = re.search(r'B-V:\s*([-+]?\d*\.?\d+)', real_sec)
        
        if m_mod and bv_mod and m_real and bv_real:
            my = float(m_mod.group(1))
            mx = float(bv_mod.group(1))
            ry = (float(m_real.group(1)) + float(m_real.group(2)))/2 if m_real.group(2) else float(m_real.group(1))
            rx = float(bv_real.group(1))
            galaxies.append({"name": name, "model": (mx, my), "real": (rx, ry)})
    return galaxies

data = parse_report(content)
fig, ax1 = plt.subplots(figsize=(14, 9))
M_SUN_BOL = 4.74

for i, g in enumerate(data, 1):
    mx, my = g["model"]; rx, ry = g["real"]
    ax1.plot([mx, rx], [my, ry], color='gray', linestyle=':', alpha=0.3, zorder=1)
    
    # Модель (номер СВЕРХУ)
    ax1.scatter(mx, my, color='crimson', s=120, edgecolors='white', zorder=4)
    ax1.text(mx, my - 0.3, str(i), fontsize=10, ha='center', va='bottom', fontweight='bold', color='crimson')
    
    # Реальность (номер СНИЗУ)
    ax1.scatter(rx, ry, color='dodgerblue', marker='x', s=100, zorder=3)
    ax1.text(rx, ry + 0.3, str(i), fontsize=9, ha='center', va='top', color='dodgerblue')

ax1.invert_yaxis()
ax1.set_xlabel(r'Color Index $(B-V)$', fontsize=12)
ax1.set_ylabel(r'Absolute Bolometric Magnitude $M_{bol}$', fontsize=12)
ax1.grid(alpha=0.3, linestyle='--')
ax1.set_title('Galaxy Correlation: Analytical Model vs Observations', fontsize=15, pad=35)

# Доп. ось Y: Светимость
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()); ax2.invert_yaxis()
ax2.set_ylabel(r'Luminosity $L/L_{\odot}$', fontsize=12, labelpad=15)
m_ticks = ax1.get_yticks()
ax2.set_yticks(m_ticks)
ax2.set_yticklabels([f"{10**(0.4*(M_SUN_BOL - m)):.1e}" for m in m_ticks])

# Доп. ось X: Температура (T = 4600 / (0.92 + BV))
ax3 = ax1.twiny()
ax3.set_xlim(ax1.get_xlim())
bv_ticks = ax1.get_xticks()
ax3.set_xticks(bv_ticks)
ax3.set_xticklabels([f"{int(4600/(0.92+bv))}K" if 0.92+bv > 0.1 else "" for bv in bv_ticks])
ax3.set_xlabel(r'Estimated Effective Temperature $T_{eff}$ (K)', fontsize=12, labelpad=15)

# Легенда с именами
legend_text = "\n".join([f"{i:2d}: {g['name']}" for i, g in enumerate(data, 1)])
plt.text(1.18, 0.5, legend_text, transform=ax1.transAxes, fontsize=8.5, family='monospace', verticalalignment='center', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

plt.tight_layout()
plt.show()
