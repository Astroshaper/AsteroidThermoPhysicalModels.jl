using Astroshaper
using Documenter

DocMeta.setdocmeta!(Astroshaper, :DocTestSetup, :(using Astroshaper); recursive=true)

makedocs(;
    modules=[Astroshaper],
    repo="https://github.com/MasanoriKanamaru/Astroshaper.jl/blob/{commit}{path}#{line}",
    sitename="Astroshaper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MasanoriKanamaru.github.io/Astroshaper.jl",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MasanoriKanamaru/Astroshaper.jl",
)
