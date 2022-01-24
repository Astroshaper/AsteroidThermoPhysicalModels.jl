

const AU = 1.495978700e11     # Astronomical Unit [m]
const G = 6.6740831e-11       # Gravitational constant [m^3/kg/s^2]
const GM☉ = 1.32712440018e20
const M☉ = GM☉ / G
const SOLAR_CONST = 1366.0    # Solar constant, Φ☉ [W/m^2]
const c₀ = 299792458.0        # speed of light [m/s]
const σ_SB = 5.670374419e-8   # Stefan–Boltzmann constant [W/m^2/K^4]


# ****************************************************************
#                              Planets
# ****************************************************************

MERCURY      = Dict{Symbol, Float64}()
MERCURY[:GM] = 2.2032e13
MERCURY[:M]  = MERCURY[:GM] / G

VENUS      = Dict{Symbol, Float64}()
VENUS[:GM] = 3.24859e14
VENUS[:M]  = VENUS[:GM] / G

EARTH      = Dict{Symbol, Float64}()
EARTH[:GM] = 3.986004418e14
EARTH[:M]  = EARTH[:GM] / G

MOON      = Dict{Symbol, Float64}()
MOON[:GM] = 4.9048695e12
MOON[:M]  = MOON[:GM] / G

MARS      = Dict{Symbol, Float64}()
MARS[:GM] = 4.282837e13
MARS[:M]  = MARS[:GM] / G

CERES      = Dict{Symbol, Float64}()
CERES[:GM] = 6.26325e10
CERES[:M]  = CERES[:GM] / G

JUPITER      = Dict{Symbol, Float64}()
JUPITER[:GM] = 1.26686534e17
JUPITER[:M]  = JUPITER[:GM] / G

SATURN      = Dict{Symbol, Float64}()
SATURN[:GM] = 3.7931187e16
SATURN[:M]  = SATURN[:GM] / G

URANUS      = Dict{Symbol, Float64}()
URANUS[:GM] = 5.793939e15
URANUS[:M]  = URANUS[:GM] / G

NEPTUNE      = Dict{Symbol, Float64}()
NEPTUNE[:GM] = 6.836529e15
NEPTUNE[:M]  = NEPTUNE[:GM] / G

PLUTO      = Dict{Symbol, Float64}()
PLUTO[:GM] = 8.71e11
PLUTO[:M]  = PLUTO[:GM] / G

ERIS      = Dict{Symbol, Float64}()
ERIS[:GM] = 1.108e12
ERIS[:M]  = ERIS[:GM] / G


# ****************************************************************
#                            Asteroids
# ****************************************************************

RYUGU      = Dict{Symbol, Float64}()
RYUGU[:GM] = 30.0
RYUGU[:M]  = RYUGU[:GM] / G