module Astroshaper

export Point, double
export load_obj

include("obj.jl")

"""
A point type
"""
struct Point
    x::Float64
end

"""
    double(p)

Double `p`.
"""
double(p) = Point(p.x * 2)

"""
Greeting.
"""
greet() = print("Hello World!")

end # module
