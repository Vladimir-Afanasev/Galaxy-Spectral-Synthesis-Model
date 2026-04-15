# Galaxy Spectral Synthesis Model

An analytical engine for simulating the Spectral Energy Distribution (SED) of galaxies. This tool leverages fundamental physics to compute galactic spectra dynamically, offering a flexible alternative to modern grid-based look-up table methods (e.g., GALAXEV or PEGASE).

## Project Purpose
The goal of this project is to provide a lightweight, transparent, and physically grounded tool for modeling galaxy evolution. Unlike static, pre-computed look-up tables, this model allows for on-the-fly adjustment of evolutionary parameters, providing immediate insight into how individual stellar components, chemical enrichment, and dust integration contribute to the final observed light.

## Authors
* **Vladimir Afanasev (MSc MIPT)** — *Lead Developer & Architect*. 
  Alumnus of the **Moscow Institute of Physics and Technology (National Research University) (MIPT)**. Conceptual design, physical model implementation, and core algorithms.
* **Nikita Korsunov** — *Assistant & Researcher*. 
  10th-grade student at the **Specialized Education and Scientific Center of the North-Caucasus Federal University (SESC NCFU)**. Data preparation, testing, and comparative analysis.

## Data Structure & Key Results

### Input Data (Parameters)
The researcher can define a galaxy by providing the following physical inputs in the **galaxies** list:
* **m**: Total stellar mass in Solar masses (M_sun).
* **type**: Morphological Hubble type (E, S0, Sa, Sb, Sc, Sd, Irr).
* **age**: Total age of the system (years).
* **feh**: Metallicity index **[Fe/H]**.
* **inc**: Inclination angle (degrees) for internal dust extinction.
* **z**: Cosmological redshift.
* **ebv_mw**: Foreground Milky Way extinction **E(B-V)**.

### Output Data
The model produces the following high-fidelity astrophysical results:
* **Absolute Bolometric Magnitude (M_bol)**.
* **Color Index (B-V)**: Johnson-Cousins photometric system.
* **Mass-to-Light Ratio (M/L)**.
* **Peak Wavelength**: The wavelength of maximum spectral energy (nm).
* **Wavelength Grid**: Logarithmic distribution for constant resolution.
* **Spectral Flux Density Plot**: High-resolution logarithmic SED visualization.

### Validation & Reporting: Theory vs. Observation
A unique feature of this code is the automated comparison module:
* **real_data Dictionary**: Researchers can input observed values (e.g., from NASA/IPAC Extragalactic Database) for specific galaxies.
* **generate_report() Function**: Automatically compares model outputs against the provided observations.
* **Output Report**: Generates a detailed text file **galaxy_analysis_report.txt** for audit and scientific verification.

## Usage Guide
1. **Model Customization**: Add new galaxy dictionaries to the **galaxies** list.
2. **Observational Data**: Update the **real_data** dictionary in the **generate_report()** function with known astronomical values.
3. **Execution**: Run the script to perform spectral integration via **Salpeter IMF** and the **Tau-model** history.
4. **Analysis**: Check the console for real-time results and open **galaxy_analysis_report.txt** for the full comparative study.

## Physical Factors & References
* **Stellar Physics:** Empirical relations from **Eker et al. (2018)**.
* **SFH & Evolution:** The **Tau-model** (exponentially declining SFR) as defined by **Tinsley (1980)**.
* **Metallicity:** Comprehensive **[Fe/H]** effects on opacity and temperature (**Worthey, 1994**).
* **Spectral Features:** Balmer Jump, Metal Blanketing (UV), and emission lines (**Kennicutt, 1998**).
* **Cosmology:** Redshift effects and interstellar extinction (**Cardelli et al., 1989; Schlegel et al., 1998**).

## License (MIT)
This project is licensed under the **MIT License**. This is a permissive license that allows for broad use while protecting the author's copyright.

**Key Permissions:**
* **Commercial use**: You may use this software and derivatives for commercial purposes.
* **Modification**: You may modify the software.
* **Distribution**: You may distribute the software.
* **Private use**: You may use and modify the software for private use.

**Conditions:**
* **Include Copyright**: A copy of the license and copyright notice must be included in all copies or substantial portions of the software.

**Limitations:**
* **No Warranty**: The software is provided "as is", without warranty of any kind.
* **No Liability**: The author cannot be held liable for any damages arising from the use of the software.

---

# Программная модель синтеза спектров галактик

Аналитический движок для динамического моделирования спектрального распределения энергии (SED) галактик, использующий фундаментальные физические принципы в качестве альтернативы сеточным методам (Look-up tables).

## Цель проекта
Предоставить легкий и физически прозрачный инструмент для моделирования эволюции галактик. В отличие от статичных таблиц, данная модель позволяет «на лету» изменять параметры эволюции, обеспечивая понимание вклада отдельных звездных компонентов и химического обогащения в итоговое излучение.

## Авторы
* **Владимир Афанасьев (магистр МФТИ)** — *Ведущий разработчик и Архитектор*. 
  Выпускник **Московского физико-технического института (национального исследовательского университета) (МФТИ)**.
* **Никита Корсунов** — *Ассистент и Исследователь*. 
  Ученик 10-го класса **Специализированного учебного научного центра Северо-Кавказского федерального университета (СУНЦ СКФУ)**.

## Данные и результаты

### Входные данные (Параметры)
Исследователь может задать галактику, указав следующие физические параметры в списке **galaxies**:
* **m**: Полная звездная масса в массах Солнца (M_sun).
* **type**: Морфологический тип Хаббла (E, S0, Sa, Sb, Sc, Sd, Irr).
* **age**: Полный возраст системы (в годах).
* **feh**: Индекс металличности **[Fe/H]**.
* **inc**: Угол наклонения (градусы) для расчета внутреннего поглощения пылью.
* **z**: Космологическое красное смещение.
* **ebv_mw**: Фоновое поглощение Млечного Пути **E(B-V)**.

### Выходные данные
Модель выдает следующие высокоточные астрофизические результаты:
* **Абсолютная болометрическая величина (M_bol)**.
* **Показатель цвета (B-V)**: Фотометрическая система Джонсона-Казинса.
* **Отношение Масса-Светимость (M/L)**.
* **Пиковая длина волны**: Длина волны максимальной спектральной энергии (нм).
* **Сетка длин волн**: Логарифмическое распределение для постоянного разрешения.
* **График спектральной плотности потока**: Высокодетализированная логарифмическая визуализация SED.

### Валидация и Отчетность: Теория  и наблюдения
Уникальной особенностью кода является модуль автоматического сравнения:
* **Словарь real_data**: Исследователи могут вводить наблюдаемые значения (напр. из базы NASA/IPAC NED) для конкретных объектов.
* **Функция generate_report()**: Автоматически сопоставляет выходные данные модели с предоставленными наблюдениями.
* **Итоговый отчет**: Генерирует подробный текстовый файл **galaxy_analysis_report.txt** для аудита и научной верификации.

## Инструкция по использованию
1. **Настройка модели**: Добавьте словари новых галактик в список **galaxies**.
2. **Данные наблюдений**: Обновите словарь **real_data** в функции **generate_report()** известными астрономическими значениями.
3. **Запуск**: Запустите скрипт для выполнения спектральной интеграции через **функцию масс Солпитера** и **Tau-модель** истории.
4. **Анализ**: Проверьте консоль для получения результатов в реальном времени и откройте **galaxy_analysis_report.txt** для полного сравнительного анализа.

## Физические факторы и ссылки
* **Звездная физика:** Эмпирические зависимости из **Eker et al. (2018)**.
* **SFH и Эволюция:** **Tau-модель** (экспоненциально затухающий SFR) согласно определению **Tinsley (1980)**.
* **Металличность:** Комплексное влияние **[Fe/H]** на непрозрачность и температуру (**Worthey, 1994**).
* **Спектральные особенности:** Бальмеровский скачок, металлическое бланкетирование (УФ) и эмиссионные линии (**Kennicutt, 1998**).
* **Космология:** Эффекты красного смещения и межзвездного поглощения (**Cardelli et al., 1989; Schlegel et al., 1998**).

## Лицензия (MIT)
Этот проект лицензирован под **лицензией MIT**. Это разрешительная лицензия, которая допускает широкое использование, защищая при этом авторские права автора.

**Основные разрешения:**
* **Коммерческое использование**: Вы можете использовать это программное обеспечение и его производные в коммерческих целях.
* **Модификация**: Вы можете изменять программное обеспечение.
* **Распространение**: Вы можете распространять программное обеспечение.
* **Частное использование**: Вы можете использовать и изменять программное обеспечение для частного использования.

**Условия:**
* **Включить авторские права**: Копия лицензии и уведомление об авторских правах должны быть включены во все копии или значительные части программного обеспечения.

**Ограничения:**
* **Отсутствие гарантий**: Программное обеспечение предоставляется «как есть», без каких-либо гарантий.
* **Отсутствие ответственности**: Автор не может нести ответственность за любые убытки, возникшие в результате использования программного обеспечения.

---
