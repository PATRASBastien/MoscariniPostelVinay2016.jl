using MoscariniPostelVinay2016
using Documenter

DocMeta.setdocmeta!(MoscariniPostelVinay2016, :DocTestSetup, :(using MoscariniPostelVinay2016); recursive=true)

makedocs(;
    modules=[MoscariniPostelVinay2016],
    authors="PATRASBastien <bastien.patras@sciencespo.fr> and contributors",
    repo="https://github.com/PATRASBastien/MoscariniPostelVinay2016.jl/blob/{commit}{path}#{line}",
    sitename="MoscariniPostelVinay2016.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://PATRASBastien.github.io/MoscariniPostelVinay2016.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PATRASBastien/MoscariniPostelVinay2016.jl",
    devbranch="main",
)
