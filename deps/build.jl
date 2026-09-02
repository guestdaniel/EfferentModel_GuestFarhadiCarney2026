using Libdl

const srcdir = joinpath(@__DIR__, "model")
const libname = "libgfc2026." * Libdl.dlext
const outpath = joinpath(@__DIR__, libname)

cd(srcdir) do
    run(`zsh compile.sh`)
end
println("Built $outpath")
