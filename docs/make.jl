using MieScattering
using Documenter

DocMeta.setdocmeta!(MieScattering, :DocTestSetup, :(using MieScattering); recursive=true)

makedocs(;
    modules=[MieScattering],
    authors="Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo="https://github.com/JuliaRemoteSensing/MieScattering.jl/blob/{commit}{path}#{line}",
    sitename="MieScattering.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaRemoteSensing.github.io/MieScattering.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/JuliaRemoteSensing/MieScattering.jl",
    devbranch="main"
)
