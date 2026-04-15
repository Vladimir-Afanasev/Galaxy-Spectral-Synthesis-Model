import math
import matplotlib.pyplot as plt

# КОНСТАНТЫ
H, C, K = 6.62607e-34, 2.99792e8, 1.38065e-23
R_SUN, T_SUN, L_SUN = 6.957e8, 5778, 3.828e26
M_SUN_BOL = 4.74
# Логарифмическая сетка длин волн для поддержания постоянного разрешения
WL_GRID = [100 * (10000/100)**(i/199) for i in range(200)]

def get_star_params(m):
    """
    Эмпирические зависимости Масса-Светимость и Масса-Радиус.
    Основано на аппроксимации данных звезд Группы Главной Последовательности.
    Ref: Eker et al. (2018) "Main-sequence luminosity-mass relation", MNRAS, 479, 5491
    https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5491E/abstract
    """
    log_m = math.log10(m)
    # Кусочно-степенной закон светимости L ~ M^alpha
    if m < 0.45: log_l = 2.307 * log_m - 0.705
    elif m < 1.05: log_l = 4.551 * log_m + 0.054
    elif m < 2.40: log_l = 4.351 * log_m + 0.024
    else: log_l = 2.745 * log_m + 1.002
    # Зависимость радиуса от массы (учет расширения звезд на ГП)
    log_r = 0.917 * log_m - 0.020 if m < 1.05 else 0.641 * log_m + 0.011
    # Возвращаем радиус и эффективную температуру (через закон Стефана-Больцмана)
    return 10**log_r, (10**log_l / (10**log_r)**2)**0.25

# Кэширование закона Планка для ускорения расчетов    
# Planck's Law (1900): B_lambda(T) = (2hc^2/lambda^5) * (1/(exp(hc/lambda kT)-1))
planck_cache = {t: [((2*math.pi*H*C**2)/((wl*1e-9)**5 * (math.exp((H*C)/(wl*1e-9*K*t))-1)) * 1e-9 if (H*C)/(wl*1e-9*K*t)<700 else 0) for wl in WL_GRID] for t in range(2000, 50050, 50)}

def get_interpolated_val(wl, row):
    # Линейная интерполяция спектра в логарифмическом пространстве
    if wl <= WL_GRID[0]: return row[0]
    if wl >= WL_GRID[-1]: return row[-1]
    pos = (math.log10(wl)-math.log10(WL_GRID[0]))/(math.log10(WL_GRID[-1])-math.log10(WL_GRID[0]))*199
    idx = int(pos); f = pos - idx
    return row[idx]*(1-f) + row[idx+1]*f

def get_filter_flux(wl_grid, spec, center, sigma):
    """
    Интегрирование потока через Гауссову кривую отклика фильтра.
    Аппроксимация фотометрической системы Johnson-Cousins (1990).
    https://ui.adsabs.harvard.edu/abs/1990PASP..102.1181B/abstract
    """
    flux = 0
    for i, wl in enumerate(wl_grid):
        response = math.exp(-((wl - center)**2) / (2 * sigma**2))
        dw = (wl_grid[i] - wl_grid[i-1]) if i > 0 else (wl_grid[1] - wl_grid[0])
        flux += spec[i] * response * dw
    return flux

def get_pop_spec(age_layer, m_target, feh=0.0, z_red=0.0, ebv_int=0.0, ebv_mw=0.0, h_type='Sb'):
    # Масса точки поворота (Main Sequence Turn-off). Звезды тяжелее m_max уже погибли.
    # Ref: Iben (1967) "Stellar Evolution Within and Off the Main Sequence", ARA&A, 5, 571
    m_max = min(50.0, max(0.1, (10**10 / age_layer)**(1/2.5)))

    # Интегральная начальная функция масс (IMF) Солпитера (Salpeter, 1955)
    # dN/dM ~ M^-2.35. https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S/abstract
    A = m_target / ((50**-0.35 - 0.1**-0.35) / -0.35)

    m_grid = [0.1 * (m_max/0.1)**(i/79) for i in range(80)]
    spec, absorbed = [0.0]*200, 0.0
    for i, M in enumerate(m_grid):
        dM = m_grid[i] - m_grid[i-1] if i > 0 else 0.001
        r_rel, t_rel = get_star_params(M)

        # Учет стадии Красных Гигантов (увеличение радиуса и падение температуры)
        # Kippenhahn, R., & Weigert, A. "Stellar Structure and Evolution"
        if M > 0.95 * m_max:
            r_rel *= 8.0
            """
            #Металличность ([Fe/H]) увеличивает непрозрачность (opacity) звездного вещества.
            Больше металлов ⟹ больше столкновений фотонов ⟹ больше давление излучения⟹
            звезда расширяется и её эффективная температура Worthey, G. (1994). "Comprehensive
            Stellar Population Models", Astrophysical Journal Supplement, 95, 107.
            """
            t_rel = (3500 - 350 * feh) / T_SUN

        t_key = max(2000, min(50000, int(round(t_rel * T_SUN / 50.0) * 50)))
        row, dn = planck_cache[t_key], A * (M**-2.35) * dM
        s_area = 4 * math.pi * (r_rel * R_SUN)**2
        for j, wl_obs in enumerate(WL_GRID):
            # Космологический Redshift: определение z через растяжение длин волн.
            wl_rest = wl_obs / (1 + z_red)
            val = get_interpolated_val(wl_rest, row) * s_area * dn

            # SPECTRAL SIEVE (Решето)
            #  Бальмеровский скачок (364.7 nm): Предел ионизации со 2-го уровня водорода.
            # Фундаментальный признак спектрального класса (Knigge et al. 2011).
            jump = 1.0 / (1.0 + math.exp(-(wl_rest - 364.7) / 2.0))
            val *= (0.80 + 0.20 * jump) # 20% обрыв

            """
            Линии поглощения (бланкетирование). Metal Blanketing: Поглощение в УФ из-за обилия линий металлов.
            падение потока коррелирует с [Fe/H]. Sandage (1969), Tinsley, B. M. (1980). "Evolution of the Stars and Gas in
            Galaxies", Fundamentals of Cosmic Physics, 5, 287.
            """
            if wl_rest < 450:
                val *= (1.0 - 0.10 * (10**feh) * math.exp(-((wl_rest - 400)**2) / 800))

            # Эмиссионные линии HII областей (H-alpha, OIII).
            # Kennicutt (1998) "Star Formation in Galaxies Across the Hubble Sequence", ARA&A, 36
            if h_type not in ['E', 'S0'] and age_layer < 4e9:
                for line_wl, amp in [(656.3, 0.6), (500.7, 0.4), (486.1, 0.2)]:
                    val += (val * amp * math.exp(-((wl_rest - line_wl)**2) / 4.0))

            # Закон межзвездного покраснения
            # Cardelli, Clayton, & Mathis (1989), ApJ, 345, 245
            ext = math.exp(-3.1 * ebv_int * (550/wl_rest)) if ebv_int > 0 else 1.0
            if j > 0: absorbed += val * (1-ext) * (WL_GRID[j]-WL_GRID[j-1])

            # Космологическое красное смещение.Эффект Доплера (спектральный)
            v_out = (val * ext) / (1 + z_red)
            """
            Значения ebv_mw (покраснение нашей Галактики) берутся из карт пыли
            Schlegel, D. J., Finkbeiner, D. P., & Davis, M. (1998).
            "Maps of Dust Infrared Emission for Cosmic Infrared Background Radiation",
            Astrophysical Journal, 500, 525. + Закон межзвездного покраснения
            Cardelli, Clayton, & Mathis (1989), ApJ, 345, 245
            """
            if ebv_mw > 0: v_out *= math.exp(-3.1 * ebv_mw * (550/wl_obs))
            spec[j] += v_out
    return spec, absorbed

def get_galaxy_spectrum(total_mass, age, mode, z, ebv_i, ebv_mw, tau, feh, h_type):
    spec, absorbed_total = [0.0]*200, 0.0
    if mode == 'elliptical':
        return get_pop_spec(age, total_mass, feh, z, 0.01, ebv_mw)
    steps = 15
    dt = age / steps
    """
    Моделирование истории звездообразования (SFH).
    Tau-model (экспоненциальный спад). Ref: Tinsley (1980) "Evolution of stars and gas in galaxies"
    https://ui.adsabs.harvard.edu/abs/1980FCPh....5..287T/abstract
    Веса для каждого временного слоя согласно экспоненциальному закону SFR ~ exp(-t/tau)
    """
    weights = [math.exp(-(i*dt)/tau) for i in range(steps)]
    w_s = sum(weights)
    for i in range(steps):
        l_age = age - (i*dt) + 1e7
        # Химическая эволюция (увеличение металличности со временем).
        #Tinsley, B. M. (1980). "Evolution of the Stars and Gas in Galaxies", Fundamentals of Cosmic Physics, 5, 287.
        l_feh = -1.0 + (feh + 1.0) * (i/steps)
        l_spec, l_abs = get_pop_spec(l_age, total_mass*(weights[i]/w_s), l_feh, z, ebv_i, ebv_mw, h_type)
        for j in range(200): spec[j] += l_spec[j]
        absorbed_total += l_abs
    return spec, absorbed_total

def run_model(m_gal, h_type, inc, age, feh, z_gal, ebv_mw_gal, name="Galaxy"):
    # Типичные времена затухания SFR для разных типов Хаббла (Sandage, 1986)
    tau_dict = {'E': 0.5e9, 'S0': 1.5e9, 'Sa': 3.0e9, 'Sb': 5.0e9, 'Sc': 10.0e9, 'Sd': 15.0e9, 'Irr': 20.0e9}
    tau = tau_dict.get(h_type, 5.0e9)
    # Геометрический фактор наклонения галактики (увеличение оптического пути сквозь пыль)
    ebv_i = 0.01 if h_type == 'E' else 0.1 / max(math.cos(math.radians(inc)), 0.15)

    spec, _ = get_galaxy_spectrum(m_gal, age, 'elliptical' if h_type=='E' else 'spiral', z_gal, ebv_i, ebv_mw_gal, tau, feh, h_type)
   
    # Фотометрические расчеты. Показатель цвета B-V и Болометрическая величина
    # sigma ~ 40 нм для фильтров Джонсона (ширина полосы около 90-100 нм)
    f_b = get_filter_flux(WL_GRID, spec, 440, 40)
    f_v = get_filter_flux(WL_GRID, spec, 550, 38)
    bv = -2.5 * math.log10(f_b / f_v) + 0.45 # Константу чуть прибавили

    # Интегрирование полной светимости (Болометрический поток)
    total_l_bol = sum(spec) * (WL_GRID[-1]-WL_GRID[0])/200 / L_SUN
    m_bol = M_SUN_BOL - 2.5 * math.log10(total_l_bol)

    # Отношение Масса-Светимость
    mass_to_light = m_gal / total_l_bol

    max_energy = max(spec) # В Вт/нм
    l_max = WL_GRID[spec.index(max_energy)]

    print(f"--- {name.upper()} ({h_type}) ---")
    print(f"M_bol: {m_bol:.2f}, B-V: {bv:.2f}, M/L: {mass_to_light:.2f}")
    print(f"Пик: {int(l_max)} нм, Макс. энергия: {max_energy:.2e} Вт/нм")
    print("-" * 30)

    # Визуализация спектральной плотности потока (SED)
    plt.figure(figsize=(8, 4))
    plt.plot(WL_GRID, spec, label=f"{name}", color='crimson')
    plt.xscale('log'); plt.yscale('log'); plt.grid(alpha=0.2); plt.legend(); plt.show()

# ЗАПУСК
galaxies = [
    {"name": "M31 Andromeda", "m": 1.5e11, "type": "Sb", "inc": 77, "age": 1.0e10, "feh": 0.0, "z": -0.001, "ebv_mw": 0.062},
    {"name": "M33 Triangulum", "m": 5.0e9, "type": "Sc", "inc": 18, "age": 8.0e9, "feh": -0.2, "z": -0.0006, "ebv_mw": 0.041},
    {"name": "M87 Giant Elliptical", "m": 1.0e12, "type": "E", "inc": 0, "age": 1.2e10, "feh": 0.3, "z": 0.0044, "ebv_mw": 0.022},
    {"name": "M104 Sombrero", "m": 2.0e11, "type": "Sa", "inc": 80, "age": 1.0e10, "feh": 0.1, "z": 0.0034, "ebv_mw": 0.045},
    {"name": "I Zwicky 18", "m": 1.0e7, "type": "Irr", "inc": 30, "age": 1.0e9, "feh": -1.5, "z": 0.0025, "ebv_mw": 0.032},
    {"name": "Milky Way", "m": 6.0e10, "type": "Sb", "inc": 0, "age": 1.3e10, "feh": 0.0, "z": 0.0, "ebv_mw": 0.0},
    {"name": "M81 Bode", "m": 7.0e10, "type": "Sa", "inc": 35, "age": 1.1e10, "feh": 0.1, "z": -0.0001, "ebv_mw": 0.072},
    {"name": "M82 Cigar", "m": 1.0e10, "type": "Irr", "inc": 60, "age": 5.0e8, "feh": -0.3, "z": 0.0007, "ebv_mw": 0.140},
    {"name": "M49 (NGC 4472)", "m": 5.0e11, "type": "E", "inc": 0, "age": 1.1e10, "feh": 0.2, "z": 0.0033, "ebv_mw": 0.020},
    {"name": "SMC", "m": 7.0e8, "type": "Irr", "inc": 40, "age": 6.0e9, "feh": -0.7, "z": 0.0005, "ebv_mw": 0.035},
    {"name": "M101 Pinwheel", "m": 1.0e11, "type": "Sc", "inc": 18, "age": 8.0e9, "feh": 0.0, "z": 0.0008, "ebv_mw": 0.008},
    {"name": "Centaurus A", "m": 1.0e12, "type": "S0", "inc": 45, "age": 1.2e10, "feh": 0.2, "z": 0.0018, "ebv_mw": 0.11},
    {"name": "LMC", "m": 1.0e10, "type": "Irr", "inc": 35, "age": 6.0e9, "feh": -0.4, "z": 0.0009, "ebv_mw": 0.07},
    {"name": "Antennae", "m": 5.0e10, "type": "Irr", "inc": 60, "age": 1.0e8, "feh": 0.1, "z": 0.005, "ebv_mw": 0.046},
    {"name": "3C 273 Quasar", "m": 1.0e13, "type": "E", "inc": 0, "age": 1.0e10, "feh": 0.3, "z": 0.158, "ebv_mw": 0.02},
    {"name": "BX442", "m": 6.0e10, "type": "Sc", "inc": 25, "age": 3.0e9, "feh": -0.5, "z": 2.18, "ebv_mw": 0.01},
    {"name": "HFLS3", "m": 5.0e10, "type": "Irr", "inc": 45, "age": 8.0e8, "feh": -1.0, "z": 6.34, "ebv_mw": 0.01},
    {"name": "CR7 Galaxy", "m": 2.0e10, "type": "Irr", "inc": 0, "age": 5.0e8, "feh": -2.0, "z": 6.60, "ebv_mw": 0.01},
    {"name": "GN-z11", "m": 1.0e9, "type": "Irr", "inc": 0, "age": 4.0e8, "feh": -2.5, "z": 10.6, "ebv_mw": 0.01},
    {"name": "NGC 1052-DF2", "m": 2.0e8, "type": "E", "inc": 0, "age": 9.0e9, "feh": -1.2, "z": 0.004, "ebv_mw": 0.01},
    {"name": "M51 Whirlpool", "m": 1.6e11, "type": "Sc", "inc": 20, "age": 9.0e9, "feh": 0.1, "z": 0.0015, "ebv_mw": 0.03},
    {"name": "M63 Sunflower", "m": 1.0e11, "type": "Sb", "inc": 56, "age": 1.0e10, "feh": 0.0, "z": 0.0016, "ebv_mw": 0.01},
    {"name": "M74 Phantom", "m": 3.0e10, "type": "Sc", "inc": 0, "age": 1.0e10, "feh": -0.1, "z": 0.0022, "ebv_mw": 0.06},
    {"name": "M110 (NGC 205)", "m": 4.0e9, "type": "E", "inc": 0, "age": 1.2e10, "feh": -0.5, "z": -0.0008, "ebv_mw": 0.05},
    {"name": "NGC 4565 Needle", "m": 1.2e11, "type": "Sb", "inc": 86, "age": 1.2e10, "feh": 0.0, "z": 0.0041, "ebv_mw": 0.01},
    {"name": "M60 Giant Elliptical", "m": 1.0e12, "type": "E", "inc": 0, "age": 1.3e10, "feh": 0.3, "z": 0.0037, "ebv_mw": 0.02},
    {"name": "M106", "m": 1.9e11, "type": "Sb", "inc": 72, "age": 1.0e10, "feh": 0.1, "z": 0.0015, "ebv_mw": 0.01},
    {"name": "NGC 1300", "m": 5.0e10, "type": "Sb", "inc": 35, "age": 1.0e10, "feh": 0.0, "z": 0.0052, "ebv_mw": 0.02},
    {"name": "M100", "m": 1.6e11, "type": "Sc", "inc": 30, "age": 1.0e10, "feh": 0.1, "z": 0.0052, "ebv_mw": 0.02},
    {"name": "NGC 2403", "m": 1.0e10, "type": "Sc", "inc": 60, "age": 1.0e10, "feh": -0.2, "z": 0.0004, "ebv_mw": 0.03}]

for g in galaxies:
    run_model(g["m"], g["type"], g["inc"], g["age"], g["feh"], g["z"], g["ebv_mw"], g["name"])

def generate_report():
    # Данные для сравнения (реальность)
    real_data = {
        "M31 Andromeda": "M_bol: -21.5, B-V: 0.92, M/L: ~5-6, Peak: ~700-900 nm",
        "M33 Triangulum": "M_bol: -19.0, B-V: 0.44, M/L: ~1-2, Peak: ~450-550 nm",
        "M87 Giant Elliptical": "M_bol: -23.5..-25, B-V: 0.96, M/L: ~5-10, Peak: ~750 nm",
        "M104 Sombrero": "M_bol: -21.8, B-V: 1.05, M/L: ~8-10, Peak: > 850 nm",
        "I Zwicky 18": "M_bol: -15.5, B-V: -0.03, M/L: < 0.5, Peak: < 300 nm",
        "Milky Way": "M_bol: -20.9, B-V: 0.65, M/L: ~2-3, Peak: ~500-700 nm",
        "M81 Bode": "M_bol: -21.1, B-V: 0.80, M/L: ~3-5, Peak: ~700 nm",
        "M82 Cigar": "M_bol: -22.0, B-V: 0.4..0.9, M/L: ~0.1-0.5, Peak: ~300-400 nm",
        "M49 (NGC 4472)": "M_bol: -22.8, B-V: 0.95, M/L: ~6-8, Peak: ~750 nm",
        "SMC": "M_bol: -17.5, B-V: 0.35, M/L: ~0.7-1.0, Peak: ~400 nm",
        "M101 Pinwheel": "M_bol: -21.0, B-V: 0.45, M/L: ~1-2, Peak: ~450-500 nm",
        "Centaurus A": "M_bol: -22.2, B-V: 0.9-1.1, M/L: ~5-10, Peak: ~600-800 nm",
        "LMC": "M_bol: -18.5, B-V: 0.50, M/L: ~0.5-1.0, Peak: ~400-500 nm",
        "Antennae": "M_bol: -21..-22, B-V: 0.3..0.5, M/L: < 0.5, Peak: ~300-500 nm",
        "3C 273 Quasar": "M_bol: -26.7 (total), B-V: 0.2..0.4, M/L: N/A, Peak: Broad (UV/Visible)",
        "BX442": "M_bol: -21.5, B-V: ~1.0 (obs), M/L: ~0.5-1.0, Peak: ~900-1100 nm (obs)",
        "HFLS3": "M_bol: -22..-23, B-V: > 1.2 (obs), M/L: < 0.5, Peak: ~1600-2000 nm (obs)",
        "CR7 Galaxy": "M_bol: -22.0, B-V: ~1.0 (obs), M/L: < 0.3, Peak: ~1500-1700 nm (obs)",
        "GN-z11": "M_bol: -21.6, B-V: > 1.3 (obs), M/L: ~0.1-0.5, Peak: ~1500-2500 nm (obs)",
        "NGC 1052-DF2": "M_bol: -15.4, B-V: 0.7-0.8, M/L: ~2-3, Peak: ~600-700 nm",
        "M51 Whirlpool": "M_bol: -21.4, B-V: 0.60, M/L: ~2-4, Peak: ~500-600 nm",
        "M63 Sunflower": "M_bol: -21.2, B-V: 0.73, M/L: ~3-5, Peak: ~650-750 nm",
        "M74 Phantom": "M_bol: -20.2, B-V: 0.48, M/L: ~1-2, Peak: ~450-550 nm",
        "M110 (NGC 205)": "M_bol: -16.5, B-V: 0.75, M/L: ~3-4, Peak: ~650-700 nm",
        "NGC 4565 Needle": "M_bol: -21.1, B-V: 0.85, M/L: ~4-6, Peak: ~750-850 nm",
        "M60 Giant Elliptical": "M_bol: -22.5, B-V: 0.98, M/L: ~6-9, Peak: ~750-800 nm",
        "M106": "M_bol: -21.8, B-V: 0.70, M/L: ~3-5, Peak: ~600-700 nm",
        "NGC 1300": "M_bol: -21.0, B-V: 0.65, M/L: ~2-4, Peak: ~550-650 nm",
        "M100": "M_bol: -21.5, B-V: 0.55, M/L: ~1-3, Peak: ~500-600 nm",
        "NGC 2403": "M_bol: -19.2, B-V: 0.45, M/L: ~1-2, Peak: ~450-500 nm"
    }

    tau_dict = {'E': 0.5e9, 'S0': 1.5e9, 'Sa': 3.0e9, 'Sb': 5.0e9, 'Sc': 10.0e9, 'Sd': 15.0e9, 'Irr': 20.0e9}

    with open("galaxy_analysis_report.txt", "w", encoding="utf-8") as f:
        f.write("ФИНАЛЬНЫЙ СРАВНИТЕЛЬНЫЙ АНАЛИЗ ГАЛАКТИК: МОДЕЛЬ VS РЕАЛЬНОСТЬ\n")
        f.write("="*75 + "\n\n")

        for g in galaxies: # Использует список словарей galaxies
            tau = tau_dict.get(g["type"], 5.0e9)
            mode = 'elliptical' if g["type"] == 'E' else 'spiral'
            ebv_i = 0.01 if g["type"] == 'E' else 0.1 / max(math.cos(math.radians(g["inc"])), 0.15)

            # Вызов основной функции
            spec, _ = get_galaxy_spectrum(g["m"], g["age"], mode, g["z"], ebv_i, g["ebv_mw"], tau, g["feh"], g["type"])

            # Расчет параметров
            # sigma ~ 40 нм для фильтров Джонсона (ширина полосы около 90-100 нм)
            f_b = get_filter_flux(WL_GRID, spec, 440, 40)
            f_v = get_filter_flux(WL_GRID, spec, 550, 38)
            bv = -2.5 * math.log10(f_b / f_v) + 0.45 # Константу чуть прибавил
            total_l = sum(spec) * (WL_GRID[-1]-WL_GRID[0])/200 / L_SUN
            m_bol = M_SUN_BOL - 2.5 * math.log10(total_l)
            ml = g["m"] / total_l
            l_max = WL_GRID[spec.index(max(spec))]

            f.write(f"{g['name']} ({g['type']})\n")
            f.write(f" МОДЕЛЬ: M_bol: {m_bol:.2f}, B-V: {bv:.2f}, M/L: {ml:.2f}, Peak: {int(l_max)} nm\n")
            f.write(f" РЕАЛЬНОСТЬ: {real_data.get(g['name'], 'Нет данных')}\n")
            f.write("-" * 65 + "\n")

    print("Отчет 'galaxy_analysis_report.txt' успешно создан.")

# Вызов функции
generate_report()
