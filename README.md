# AFFIRM.jl 
Alberta Farm Fertilizer Information Recommendation Manager (AFFIRM). AFFIRM has 3 concurrent versions - AFFIRM v3.0, AFFIRM-R and AFFIRM.jl. The AFFIRM v3.0 is the production version of AFFIRM. The AFFIRM-R is the nitrogen sub-model of AFFIRM v3.0. The AFFIRM.jl provides batch run utilites for multiple scenarios and provides similar functionalities as in AFFIRM-R.

## Get started
Install AFFIRM.jl
```julia
using("Pkg");
Pkg.add("AFFIRM.jl");
```
Create AFFIRM folders
```julia
julia> using AFFIRM;
julia> cd("your/project/folder")
julia> AFFIRM.create_affirm()
```
