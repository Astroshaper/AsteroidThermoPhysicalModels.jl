

const AU = SPICE.convrt(1, "au", "m")  # Astronomical Unit [m]
const G = 6.6740831e-11                # Gravitational constant [m^3/kg/s^2]
const GM☉ = 1.32712440018e20
const M☉ = GM☉ / G
const SOLAR_CONST = 1366.0             # Solar constant, Φ☉ [W/m^2]
const c₀ = 299792458.0                 # Speed of light [m/s]
const σ_SB = 5.670374419e-8            # Stefan–Boltzmann constant [W/m^2/K^4]


# ****************************************************************
#                              Planets
# ****************************************************************

# TODO: replace these dictionaries with structs
"""
# Mercury
"""
MERCURY      = Dict{Symbol, Float64}()
MERCURY[:GM] = 2.2032e13
MERCURY[:M]  = MERCURY[:GM] / G

"""
# Venus
"""
VENUS      = Dict{Symbol, Float64}()
VENUS[:GM] = 3.24859e14
VENUS[:M]  = VENUS[:GM] / G

"""
# Earth
"""
EARTH      = Dict{Symbol, Float64}()
EARTH[:GM] = 3.986004418e14
EARTH[:M]  = EARTH[:GM] / G

"""
# Mars
"""
MARS      = Dict{Symbol, Float64}()
MARS[:GM] = 4.282837e13
MARS[:M]  = MARS[:GM] / G

"""
# Jupiter
"""
JUPITER      = Dict{Symbol, Float64}()
JUPITER[:GM] = 1.26686534e17
JUPITER[:M]  = JUPITER[:GM] / G

"""
# Saturn
"""
SATURN      = Dict{Symbol, Float64}()
SATURN[:GM] = 3.7931187e16
SATURN[:M]  = SATURN[:GM] / G

"""
# Uranus
"""
URANUS      = Dict{Symbol, Float64}()
URANUS[:GM] = 5.793939e15
URANUS[:M]  = URANUS[:GM] / G

"""
# Neptune
"""
NEPTUNE      = Dict{Symbol, Float64}()
NEPTUNE[:GM] = 6.836529e15
NEPTUNE[:M]  = NEPTUNE[:GM] / G


# ****************************************************************
#                            Dwarf planets
# ****************************************************************

"""
# Ceres
"""
CERES      = Dict{Symbol, Float64}()
CERES[:GM] = 6.26325e10
CERES[:M]  = CERES[:GM] / G

"""
# Pluto
"""
PLUTO      = Dict{Symbol, Float64}()
PLUTO[:GM] = 8.71e11
PLUTO[:M]  = PLUTO[:GM] / G

"""
# Eris
"""
ERIS      = Dict{Symbol, Float64}()
ERIS[:GM] = 1.108e12
ERIS[:M]  = ERIS[:GM] / G


# ****************************************************************
#                            Satellites
# ****************************************************************

"""
# Moon
"""
MOON      = Dict{Symbol, Float64}()
MOON[:GM] = 4.9048695e12
MOON[:M]  = MOON[:GM] / G


# ****************************************************************
#                            Asteroids
# ****************************************************************

"""
# Asteroid 162173 Ryugu

## Physical parameters
- `:GM` : GM
- `:M`  : Mass
- `:μ`  : Standard gravitational parameter about Sun and Ryugu

## Orbital elements
- `:a`  : Semi-mojor axis [AU]
- `:e`  : Eccentricity [-]
- `:I`  : Inclination [deg]
- `:Ω`  : Longitude of the ascending node [deg]
- `:ω`  : Argument of periapsis [deg]
- `:Φ`  : # Mean anomaly [deg]

## Spin parameters
- `:α`  : Right ascension (RA) in equatorial coordinate system [deg]
- `:δ`  : Declination (Dec) in equatorial coordinate system [deg]
- `:P`  : Rotation period [h]
"""
RYUGU      = Dict{Symbol, Float64}()
RYUGU[:GM] = 30.0              # GM
RYUGU[:M]  = RYUGU[:GM] / G    # Mass
RYUGU[:μ]  = GM☉ + RYUGU[:GM]  # Standard gravitational parameter about Sun and Ryugu

RYUGU[:a] = 1.18956373  # Semi-mojor axis [AU]
RYUGU[:e] = 0.19027921  # Eccentricity
RYUGU[:I] = 5.8840222   # Inclination [deg]
RYUGU[:Ω] = 251.589203  # Longitude of the ascending node [deg]
RYUGU[:ω] = 211.435963  # Argument of periapsis [deg]
RYUGU[:Φ] = 21.9353799  # Mean anomaly [deg]

RYUGU[:α] = 96.4        # Right ascension (RA)
RYUGU[:δ] = -66.4       # Declination (Dec)
RYUGU[:P] = 7.63262     # Rotation period


"""
# Asteroid 65803 Didymos (1996 GT)

## Physical parameters
- `:GM` : GM
- `:M`  : Mass
- `:μ`  : Standard gravitational parameter about Sun and Ryugu

## Orbital elements
- `:a`  : Semi-mojor axis [AU]
- `:e`  : Eccentricity [-]
- `:I`  : Inclination [deg]
- `:Ω`  : Longitude of the ascending node [deg]
- `:ω`  : Argument of periapsis [deg]
- `:Φ`  : # Mean anomaly [deg]

## Spin parameters
- `:α`  : Right ascension (RA) in equatorial coordinate system [deg]
- `:δ`  : Declination (Dec) in equatorial coordinate system [deg]
- `:P`  : Rotation period [h]

## References
- Small Body Database Lookup: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=didymos
- Naidu et al. (2020)
"""
DIDYMOS      = Dict{Symbol, Float64}()
DIDYMOS[:GM] = 0.             # GM
DIDYMOS[:M]  = DIDYMOS[:GM] / G    # Mass
DIDYMOS[:μ]  = GM☉ + DIDYMOS[:GM]  # Standard gravitational parameter about Sun and Ryugu

DIDYMOS[:a] = 1.644206917783226   # Semi-mojor axis [AU]
DIDYMOS[:e] = 0.3838482230328272  # Eccentricity
DIDYMOS[:I] = 3.407952385532195   # Inclination [deg]
DIDYMOS[:Ω] = 73.19580068754534   # Longitude of the ascending node [deg]
DIDYMOS[:ω] = 319.3229511490275   # Argument of periapsis [deg]
DIDYMOS[:Φ] = 232.0090776515427   # Mean anomaly [deg]

DIDYMOS[:α] = 78.     # Right ascension (RA)
DIDYMOS[:δ] = -71.    # Declination (Dec)
# DIDYMOS[:λ] = 290.
# DIDYMOS[:β] = -75.
DIDYMOS[:P] = 2.2593  # Rotation period


"""
# Dimorphos

- `:a` : Semi-major axis of the mutual orbit with Dimorphos [m]
- `:P` : Rotation periof of the mutual orbit with Dimorphos [h]
"""
DIMORPHOS     = Dict{Symbol, Float64}()
DIMORPHOS[:a] = 1190.

DIMORPHOS[:α] = 78.    # Right ascension (RA)
DIMORPHOS[:δ] = -71.   # Declination (Dec)
DIMORPHOS[:P] = 11.93  # Rotation period


# ****************************************************************
#                              Comets
# ****************************************************************

