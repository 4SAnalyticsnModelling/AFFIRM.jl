"""
```julia
writeoutputs(xs...)
```
This function organizes and writes the outputs in order. ```xs...``` represents any sets of arguments.
"""
function writeoutputs(xs :: Vararg) :: String
    outputs_ = replace(join(xs, ','), r"\[|\]" => "")
    return outputs_
end
Uniform_(low_lim :: Float32, up_lim :: Float32) :: Vector{Float32} = Uniform(low_lim, up_lim)
Normal_(mean_ :: Float32, stdev_ :: Float32) :: Vector{Float32} = Normal(mean_, stdev_)
"""
```julia
get_distribution(text_in)
```
This function organizes the inputs for composite simulations (e.g. Monte-Carlo or Step-wise) for numerical variables by using user selections.
"""
function get_distribution(text_in :: Union{String, SubString{String}}) :: Vector{Float32}
    line = split(text_in, "|")
    if length(line) > 1
        distribution_id = parse(Int64, line[length(line) - 1])
        if distribution_id == 1
            low_lim, up_lim = parse.(Float32, line[1:2])
            step_size = parse(Float32, line[4])
            out_arr = Array{Float32}(collect(low_lim:step_size:up_lim))
        elseif distribution_id == 2
            low_lim, up_lim = parse.(Float32, line[1:2])
            n_sim = parse(Int64, line[4])
            out_arr = rand(Float32, Uniform_(low_lim, up_lim), n_sim)
        elseif distribution_id == 3
            mean_, stdev_ = parse.(Float32, line[1:2])
            n_sim = parse(Int64, line[4])
            out_arr = rand(Float32, Normal_(mean_, stdev_), n_sim)
        end
    else
        out_arr = [parse(Float32, line[1])]
    end
    out_arr .= round.(out_arr, digits = 2)
    return out_arr
end
"""
```julia
get_combined_simulation(text_in)
```
This function organizes the inputs for composite simulations for categorical variables by using user selections.
"""
function get_combined_simulation(text_in :: Union{String, SubString{String}}) :: Vector{Int64}
    line = split(text_in, "|")
    out_arr = parse.(Int64, line)
    return out_arr
end