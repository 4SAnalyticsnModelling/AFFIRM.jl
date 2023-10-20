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
"""
```julia
get_distribution(text_in)
```
This function organizes the inputs for composite simulations (e.g. Monte-Carlo or Step-wise) for numerical variables by using user selections.
"""
function get_distribution(text_in :: Union{String, SubString{String}}) :: Union{SArray, Array}
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
            out_arr = Array{Float32}(rand(Uniform(low_lim, up_lim), n_sim))
        elseif distribution_id == 3
            mean_, stdev_ = parse.(Float32, line[1:2])
            n_sim = parse(Int64, line[4])
            out_arr = Array{Float32}(rand(Normal(mean_, stdev_), n_sim))
        end
    else
        out_arr = [parse(Float32, line[1])]
    end
    out_arr .= round.(out_arr, digits = 2)
    if length(out_arr) <= 100
        out_arr = SArray{Tuple{size(out_arr)...}, eltype(out_arr)}(out_arr)
    end
    return out_arr
end
"""
```julia
get_combined_simulation(text_in)
```
This function organizes the inputs for composite simulations for categorical variables by using user selections.
"""
function get_combined_simulation(text_in :: Union{String, SubString{String}}) :: Union{SArray, Array}
    line = split(text_in, "|")
    out_arr = parse.(Int64, line)
    if length(out_arr) <= 100
        out_arr = SArray{Tuple{size(out_arr)...}, eltype(out_arr)}(out_arr)
    end
    return out_arr
end
