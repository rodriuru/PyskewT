# -*- coding: utf-8 -*-
"""Physical and meteorological constants for PySkewT v2."""

# Gas constants (J / (kg * K))
RD = 287.053           # Specific gas constant for dry air
RV = 461.51            # Specific gas constant for water vapor
EPSILON = RD / RV      # Ratio of gas constants (~0.62198)

# Specific heat capacities (J / (kg * K))
CP_D = 1004.67         # Specific heat of dry air at constant pressure
CV_D = CP_D - RD       # Specific heat of dry air at constant volume (~717.6)
CP_V = 1875.0          # Specific heat of water vapor
C_L = 4186.0           # Specific heat of liquid water

# Poisson constant
KAPPA = RD / CP_D      # Rd / cp (~0.2857)

# Latent heat of vaporization (J / kg)
LV0 = 2.501e6          # Latent heat of vaporization at 0 deg C

# Standard thermodynamic reference points
T_ZERO = 273.15        # 0 deg C in Kelvin
P_STD = 1000.0         # Standard reference pressure in hPa (1000 hPa)
E_ZERO = 6.112         # Saturation vapor pressure at 0 deg C in hPa

# Standard gravity
GRAVITY = 9.80665      # m / s^2

# Skew-T diagram default parameters
DEFAULT_SKEW_ANGLE = 45.0  # Angle of isotherms (degrees)
DEFAULT_P_BOT = 1050.0     # Bottom pressure in hPa
DEFAULT_P_TOP = 100.0      # Top pressure in hPa
