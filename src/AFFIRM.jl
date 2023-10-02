module AFFIRM

using Git
using BSON
using Distributed
using OffsetArrays
using Distributions
const git = Git.git()

const residue_management_multiplyer = [1.0f0, 0.0f0, 0.0f0]
const b0irrig = 519.83f0
const b1irrig = 1.56f0
const low_n_rate = 0.0f0
const n_rate_step_size = 10.0f0
const ns_max = 350.0f0
const kg_ha_n_lb_ac = 1.12f0
    
export run_affirm, create_affirm

include("utils.jl")
include("affirm-batch.jl")
end
