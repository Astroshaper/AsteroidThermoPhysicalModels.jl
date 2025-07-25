using AsteroidThermoPhysicalModels
using Documenter

DocMeta.setdocmeta!(AsteroidThermoPhysicalModels, :DocTestSetup, :(using AsteroidThermoPhysicalModels); recursive=true)

makedocs(;
    modules=[AsteroidThermoPhysicalModels],
    repo="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl/blob/{commit}{path}#{line}",
    sitename="AsteroidThermoPhysicalModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Astroshaper.github.io/AsteroidThermoPhysicalModels.jl",
        assets=["assets/favicon.ico"],
        repolink="https://github.com/Astroshaper/AsteroidThermoPhysicalModels.jl"
    ),
    pages = [
        "Introduction" => "index.md",
        "Physical model" => "physical_model.md",
        "Usage Examples" => "examples.md",
        "Performance" => "benchmarks.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/Astroshaper/AsteroidThermoPhysicalModels.jl",
)
