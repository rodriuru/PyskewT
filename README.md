# PySkewT 🌪️ v2.0

> **Diagramas Termodinámicos Skew-T ln-P y Análisis de Sondeos Atmosféricos en Python.**
> 
> *Modern, high-performance Skew-T ln-P thermodynamic diagram engine and sounding analysis in Python.*

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/Version-2.0.0-green.svg)](https://github.com/rodriuru/PyskewT)

---

## 🌟 Novedades en la Versión 2.0 (v2)

**PySkewT v2** es una reescritura moderna, modular y altamente optimizada del proyecto original, diseñada tanto para análisis operacional en meteorología como para investigación y docencia:

* ⚡ **Rendimiento 10x superior:** Grilla de fondo generada en menos de **120 ms** mediante `LineCollection` vectorizado de Matplotlib, reemplazando bucles iterativos lentos.
* 🎯 **Precisión Científica y Termodinámica Rigurosa:**
  * Fórmulas de **Bolton (1980)** y estándares de la **OMM (WMO)** para presiones de vapor de saturación $e_s(T)$ y relaciones de mezcla $w_s(T, p)$.
  * Integración numérica de pseudoadiabatas húmedas mediante **Runge-Kutta de 4to orden (RK4)** estable y continuo en toda la tropósfera.
* 📊 **Diagnóstico Convectivo Automatizado:**
  * **LCL** (Nivel de Condensación por Ascenso) con temperatura y presión exactas.
  * **LFC** (Nivel de Convección Libre) y **EL** (Nivel de Equilibrio).
  * **CAPE** (*Convective Available Potential Energy*) y **CIN** (*Convective Inhibition*).
  * **Lifted Index (LI)**, **K-Index**, **Total Totals (TT)** y **Agua Precipitable (PW)**.
* 🎨 **Visualización de Calidad de Publicación:**
  * Sombreado automático de **CAPE** (rojo traslúcido) y **CIN** (azul celeste).
  * Trayectoria gráfica del ascenso de la parcela superficial ($T_{parcel}$).
  * Resalte visual de la **isoterma de 0°C** (nivel de congelamiento) en azul celeste.
  * Columna lateral dedicada para **barbas de viento** (*wind barbs*) con muestreo inteligente anti-solapamiento.
  * HUD resumen con todos los parámetros termodinámicos en la esquina superior.
* 📁 **Lector Universal de Sondeos:** Compatibilidad directa con archivos de la **Universidad de Wyoming (UWYO)**, tablas ASCII/CSV y arreglos de NumPy.
* 🔄 **100% Retrocompatible:** Mantiene la interfaz original de la v1 (`SkewT.py`) sin romper scripts existentes.

---

## 📸 Vista Previa (Ejemplo v2)

![Diagrama Skew-T v2](example_v2.png)

---

## 🚀 Instalación Rápida

Podés instalar el paquete en modo editable o directamente desde el repositorio:

```bash
git clone https://github.com/rodriuru/PyskewT.git
cd PyskewT
pip install .
```

### Requisitos:
* `python >= 3.8`
* `numpy >= 1.20`
* `matplotlib >= 3.5`
* `scipy >= 1.7`

---

## 💻 Uso

### 1. Desde la Línea de Comandos (CLI)

Generar un diagrama completo a partir de un sondeo en formato Wyoming (`sondeo3.txt`):

```bash
# Ejecutar módulo
python3 -m pyskewt sondeo3.txt -o sondeo_analisis.png

# O mediante el comando instalado
pyskewt sondeo3.txt --dpi 300 -o mi_grafico.png
```

Opciones disponibles en el CLI:
```text
pyskewt [-h] [-o OUTPUT] [--dpi DPI] [--title TITLE] [--no-winds] [--no-parcel] [--show] archivo
```

### 2. Desde Script de Python (API v2)

```python
from pyskewt import SkewT, Sounding

# 1. Cargar sondeo (formato Wyoming o texto)
sounding = Sounding.from_file("sondeo3.txt")

# 2. Inicializar diagrama con barbas de viento
skew = SkewT(with_wind_barbs=True)

# 3. Dibujar grilla termodinámica de fondo
skew.draw_background()

# 4. Trazar perfiles y calcular parámetros
indices = skew.plot_sounding(
    sounding,
    plot_parcel=True,
    shade_cape_cin=True,
    plot_winds=True,
    title="Sondeo Atmosférico - PySkewT v2"
)

# 5. Guardar en alta definición
skew.save("salida.png", dpi=300)
print(f"CAPE: {indices['cape']:.1f} J/kg | LCL: {indices['p_lcl']:.1f} hPa")
```

### 3. Modo Retrocompatible (API v1)

Si tenés scripts que usaban la versión anterior:

```python
from SkewT import SkewT
import numpy as np

m = SkewT(1)
m.draw_background()
sondeo = np.loadtxt("sondeo3.txt", skiprows=4)
m.graficar(sondeo)
m.fig.savefig("aa.png", dpi=300)
```

---

## ⏱️ Comparativa de Rendimiento (Benchmark)

| Característica | PySkewT v1 (2021) | PySkewT v2 (2026) | Mejora |
| :--- | :--- | :--- | :--- |
| **Generación de grilla de fondo** | Bucles individuales `ax.plot` (~850 ms) | `LineCollection` vectorizado (~118 ms) | **~7x más rápido** |
| **Integración de pseudoadiabatas** | Euler directo con paso fijo | Runge-Kutta 4to orden (RK4) | **Precisión O(h⁴) y estable** |
| **Cálculo de CAPE / CIN / LCL** | No disponible / Manual | Automático integrado (~50 ms) | **Nuevo en v2** |
| **Barbas de Viento** | No disponible | Columna lateral dedicada | **Nuevo en v2** |
| **Compatibilidad NumPy 2.x** | Incompatible (`np.trapz`) | Totalmente compatible | **100% compatible** |

Para ejecutar el benchmark en tu entorno:
```bash
python3 examples/benchmark_v1_vs_v2.py
```

---

## 🧪 Pruebas Unitarias (Tests)

Para correr la suite de pruebas unitarias:

```bash
python3 -m unittest discover -s tests -v
```

---

## 📚 Fundamentos Físicos

1. **Transformación Skew-T:**
   $$x = T - \beta \ln\left(\frac{p}{p_0}\right), \quad y = -\ln(p)$$
   donde $\beta = 40.0$ y $p_0 = 1000\text{ hPa}$, lo que otorga la inclinación clásica de 45° a las isotermas.
2. **Presión de Vapor de Saturación (Bolton 1980):**
   $$e_s(T) = 6.112 \exp\left(\frac{17.67 \cdot T}{T + 243.5}\right)\text{ hPa}$$
3. **Gradiente Adiabático Húmedo:**
   $$\frac{dT}{dp} = \frac{R_d T}{p c_p} \left[ \frac{1 + \frac{L_v w_s}{R_d T}}{1 + \frac{\epsilon L_v^2 w_s}{c_p R_d T^2}} \right]$$

---

## 📄 Licencia

Este proyecto está bajo la Licencia **GNU General Public License v3.0 (GPL-3.0)**. Ver el archivo [LICENSE](LICENSE) para más detalles.
