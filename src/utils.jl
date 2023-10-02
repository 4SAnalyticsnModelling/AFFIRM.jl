"""
```julia
writeoutputs(xs...)
```
This function organizes and writes the outputs in order. ```xs...``` represents any sets of arguments.
"""
function writeoutputs(xs...)
    j = ""
    if length(xs) > 1 xs = copy(vcat(xs...)) end
    for x in xs
        if j == ""
            j = string(x)
        else
            j = j * "," * string(x)
        end
    end
    outputs_ :: String = replace(replace(replace(j * "\n", r"\[|\]" => ""), r"\.0\,|\.00\," => ","), r"\.0\n|\.00\n" => "\n")
    return outputs_
end
"""
```julia
get_distribution(text_in)
```
This function organizes the inputs for composite simulations (e.g. Monte-Carlo or Step-wise) for numerical variables by using user selections.
"""
function get_distribution(text_in)
    line = split(text_in, "|")
    if length(line) > 1
        distribution_id = parse(Int, line[length(line) - 1])
        if distribution_id == 1
            low_lim, up_lim = parse.(Float32, line[1:2])
            step_size = parse(Float32, line[4])
            out_arr = Array{Float32}(collect(low_lim:step_size:up_lim))
        elseif distribution_id == 2
            low_lim, up_lim = parse.(Float32, line[1:2])
            n_sim = parse(Int, line[4])
            out_arr = Array{Float32}(rand(Uniform(low_lim, up_lim), n_sim))
        elseif distribution_id == 3
            mean_, stdev_ = parse.(Float32, line[1:2])
            n_sim = parse(Int, line[4])
            out_arr = Array{Float32}(rand(Normal(mean_, stdev_), n_sim))
        end
    else
        out_arr = [parse(Float32, line[1])]
    end
    return round.(out_arr, digits = 2)
end
"""
```julia
get_combined_simulation(text_in)
```
This function organizes the inputs for composite simulations for categorical variables by using user selections.
"""
function get_combined_simulation(text_in)
    line = split(text_in, "|")
    out_arr = parse.(Int, line)
    return out_arr
end
