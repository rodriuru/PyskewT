# -*- coding: utf-8 -*-
"""Command-line interface for PySkewT v2."""

from __future__ import annotations
import argparse
import sys
from .skewt import SkewT
from .sounding import Sounding


def main() -> int:
    parser = argparse.ArgumentParser(
        description="PySkewT v2: Generador de Diagramas Termodinámicos Skew-T ln-P para Meteorología"
    )
    parser.add_argument("archivo", help="Ruta al archivo de sondeo (formato Wyoming o columnas numéricas)")
    parser.add_argument("-o", "--output", default="skewt_output.png", help="Nombre del archivo de imagen de salida (default: skewt_output.png)")
    parser.add_argument(
        "--theme",
        choices=["modern", "dark", "classic", "accessible"],
        default="modern",
        help="Paleta de colores para el diagrama: modern, dark, classic, accessible (default: modern)"
    )
    parser.add_argument("--dpi", type=int, default=300, help="Resolución DPI de la imagen (default: 300)")
    parser.add_argument("--title", default="Sondeo Atmosférico - Skew-T ln-P", help="Título del diagrama")
    parser.add_argument("--no-winds", action="store_true", help="Desactivar trazado de barbas de viento")
    parser.add_argument("--no-parcel", action="store_true", help="Desactivar ascenso de parcela y cálculo de CAPE/CIN")
    parser.add_argument("--show", action="store_true", help="Mostrar ventana interactiva")

    args = parser.parse_args()

    print(f"Cargando sondeo desde: {args.archivo}...")
    try:
        sounding = Sounding.from_file(args.archivo)
    except Exception as e:
        print(f"Error al leer el archivo de sondeo: {e}", file=sys.stderr)
        return 1

    print(f"Sondeo cargado: {len(sounding)} niveles (desde {sounding.pressure[0]:.1f} hPa hasta {sounding.pressure[-1]:.1f} hPa)")
    
    skew = SkewT(with_wind_barbs=not args.no_winds, theme=args.theme)
    skew.draw_background()

    indices = skew.plot_sounding(
        sounding,
        plot_parcel=not args.no_parcel,
        shade_cape_cin=not args.no_parcel,
        plot_winds=not args.no_winds,
        title=args.title,
    )

    print("\nResultados Termodinámicos:")
    print(f"  • CAPE:  {indices['cape']:.1f} J/kg")
    print(f"  • CIN:   {indices['cin']:.1f} J/kg")
    print(f"  • LCL:   {indices['p_lcl']:.1f} hPa ({indices['t_lcl']:.1f}°C)")
    if indices['p_lfc']:
        print(f"  • LFC:   {indices['p_lfc']:.1f} hPa")
    if indices['p_el']:
        print(f"  • EL:    {indices['p_el']:.1f} hPa")
    print(f"  • Lifted Index (LI): {indices['lifted_index']:.1f}°C")
    print(f"  • K-Index:           {indices['k_index']:.1f}°C")
    print(f"  • Total Totals:      {indices['total_totals']:.1f}")
    print(f"  • Agua Precipitable: {indices['precipitable_water']:.1f} mm")

    print(f"\nGuardando diagrama en: {args.output} (DPI={args.dpi})...")
    skew.save(args.output, dpi=args.dpi)
    print("¡Listo con éxito!")

    if args.show:
        skew.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
