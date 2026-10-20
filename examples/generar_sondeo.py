# -*- coding: utf-8 -*-
"""Ejemplo de uso de PySkewT v2: Carga de sondeo y generación de diagrama."""

from pathlib import Path
import sys

# Permitir ejecución directa sin haber instalado el paquete
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pyskewt import SkewT, Sounding


def main():
    # 1. Cargar el sondeo meteorológico desde el archivo
    archivo_sondeo = Path(__file__).resolve().parent.parent / "sondeo3.txt"
    print(f"Leyendo sondeo: {archivo_sondeo.name}...")
    sounding = Sounding.from_file(archivo_sondeo)

    # 2. Inicializar el diagrama Skew-T con columna de vientos
    skew = SkewT(with_wind_barbs=True)

    # 3. Dibujar la grilla termodinámica de fondo
    # (isobaras, isotermas con 0°C destacado, adiabatas secas, pseudoadiabatas y líneas de mezcla)
    print("Dibujando fondo termodinámico...")
    skew.draw_background()

    # 4. Trazar el sondeo (temperatura, rocío, parcela de superficie, sombreado CAPE/CIN y barbas)
    print("Trazando perfiles atmosféricos y calculando índices convectivos...")
    indices = skew.plot_sounding(
        sounding,
        plot_parcel=True,
        shade_cape_cin=True,
        plot_winds=True,
        title="Sondeo Atmosférico - Skew-T ln-P (PySkewT v2)"
    )

    # 5. Imprimir resumen meteorológico
    print("\n" + "=" * 45)
    print("  ÍNDICES TERMODINÁMICOS CALCULADOS")
    print("=" * 45)
    print(f"• Presión en Superficie:  {indices['p_sfc']:.1f} hPa")
    print(f"• Temp. Superficie / Td:  {indices['t_sfc']:.1f}°C / {indices['td_sfc']:.1f}°C")
    print(f"• Nivel Condensación (LCL): {indices['p_lcl']:.1f} hPa ({indices['t_lcl']:.1f}°C)")
    if indices['p_lfc']:
        print(f"• Nivel Convección Libre: {indices['p_lfc']:.1f} hPa")
    if indices['p_el']:
        print(f"• Nivel de Equilibrio (EL): {indices['p_el']:.1f} hPa")
    print(f"• CAPE:                   {indices['cape']:.1f} J/kg")
    print(f"• CIN:                    {indices['cin']:.1f} J/kg")
    print(f"• Lifted Index (LI):      {indices['lifted_index']:.1f}°C")
    print(f"• K-Index:                {indices['k_index']:.1f}°C")
    print(f"• Total Totals:           {indices['total_totals']:.1f}")
    print(f"• Agua Precipitable (PW): {indices['precipitable_water']:.1f} mm")
    print("=" * 45)

    # 6. Guardar la figura en alta resolución
    salida = "ejemplo_sondeo_v2.png"
    skew.save(salida, dpi=300)
    print(f"\n¡Diagrama guardado exitosamente en '{salida}' (DPI=300)!")


if __name__ == "__main__":
    main()
