using PopSimIBX
using Documenter

DocMeta.setdocmeta!(PopSimIBX, :DocTestSetup, :(using PopSimIBX); recursive=true)

makedocs(;
    modules=[PopSimIBX],
    authors="Peter Arndt <arndt@molgen.mpg.de> and contributors",
    sitename="PopSimIBX.jl",
    format=Documenter.HTML(;
        canonical="https://ArndtLab.github.io/PopSimIBX.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArndtLab/PopSimIBX.jl",
    devbranch="main",
)
