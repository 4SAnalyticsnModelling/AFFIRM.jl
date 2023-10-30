struct AFFIRM_coeff
    crop_name :: Vector{String}
    previous_crop :: Vector{String}
    previous_crop_yld_unit :: Vector{String}
    residue_management :: Vector{String}
    soil_zone :: Vector{String}
    n_source :: Vector{String}
    n_time :: Vector{String}
    n_place :: Vector{String}
    soil_texture :: Vector{String}
    spring_moisture_condition :: Vector{String}
    irrigation_flag :: Vector{String}
    n_source_percent_n :: Vector{Float32}
    wue :: Array{Float32, 5}
    epsilon :: Array{Float32, 5}
    nminus1 :: Array{Float32, 5}
    crop_unit_conv_coef :: Vector{Float32}
    spring_soil_moisture :: Matrix{Float32}
    b0ph :: Vector{Float32}
    b1ph :: Vector{Float32}
    b2ph :: Vector{Float32}
    phmax :: Vector{Float32}
    phmin :: Vector{Float32}
    b0ec :: Vector{Float32}
    b1ec :: Vector{Float32}
    b0precip :: OffsetArray{Float32, 3, Array{Float32, 3}}
    b1precip :: OffsetArray{Float32, 3, Array{Float32, 3}}
    soil_zone_id :: OffsetArray{Int64, 3, Array{Int64, 3}}
    b0ag :: Matrix{Float32}
    b0bg :: Matrix{Float32}
end

function AFFIRM_coeff(data_file_path :: String)
    crop_name = BSON.load(data_file_path * "crop_name.bson")[Symbol("crop_name")]
    previous_crop = BSON.load(data_file_path * "previous_crop.bson")[Symbol("previous_crop")]
    previous_crop_yld_unit = BSON.load(data_file_path * "previous_crop_yld_unit.bson")[Symbol("previous_crop_yld_unit")]
    residue_management = BSON.load(data_file_path * "residue_management.bson")[Symbol("residue_management")]
    soil_zone = BSON.load(data_file_path * "soil_zone.bson")[Symbol("soil_zone")]
    n_source = BSON.load(data_file_path * "n_source.bson")[Symbol("n_source")]
    n_time = BSON.load(data_file_path * "n_time.bson")[Symbol("n_time")]
    n_place = BSON.load(data_file_path * "n_place.bson")[Symbol("n_place")]
    soil_texture = BSON.load(data_file_path * "soil_texture.bson")[Symbol("soil_texture")]
    spring_moisture_condition = BSON.load(data_file_path * "spring_moisture_condition.bson")[Symbol("spring_moisture_condition")]
    irrigation_flag = BSON.load(data_file_path * "irrigation_flag.bson")[Symbol("irrigation_flag")]
    n_source_percent_n = BSON.load(data_file_path * "n_source_percent_n.bson")[Symbol("n_source_percent_n")]
    wue = BSON.load(data_file_path * "wue.bson")[Symbol("wue")]
    epsilon = BSON.load(data_file_path * "epsilon.bson")[Symbol("epsilon")]
    nminus1 = BSON.load(data_file_path * "nminus1.bson")[Symbol("nminus1")]
    crop_unit_conv_coef = BSON.load(data_file_path * "crop_unit_conv_coef.bson")[Symbol("crop_unit_conv_coef")]
    spring_soil_moisture = BSON.load(data_file_path * "spring_soil_moisture.bson")[Symbol("spring_soil_moisture")]
    b0ph = BSON.load(data_file_path * "b0ph.bson")[Symbol("b0ph")]
    b1ph = BSON.load(data_file_path * "b1ph.bson")[Symbol("b1ph")]
    b2ph = BSON.load(data_file_path * "b2ph.bson")[Symbol("b2ph")]
    phmax = BSON.load(data_file_path * "phmax.bson")[Symbol("phmax")]
    phmin = BSON.load(data_file_path * "phmin.bson")[Symbol("phmin")]
    b0ec = BSON.load(data_file_path * "b0ec.bson")[Symbol("b0ec")]
    b1ec = BSON.load(data_file_path * "b1ec.bson")[Symbol("b1ec")]
    b0precip = BSON.load(data_file_path * "b0precip.bson")[Symbol("b0precip")]
    b1precip = BSON.load(data_file_path * "b1precip.bson")[Symbol("b1precip")]
    soil_zone_id = BSON.load(data_file_path * "soil_zone_id.bson")[Symbol("soil_zone_id")]
    b0ag = BSON.load(data_file_path * "b0ag.bson")[Symbol("b0ag")]
    b0bg = BSON.load(data_file_path * "b0bg.bson")[Symbol("b0bg")]
    return AFFIRM_coeff(crop_name, previous_crop, previous_crop_yld_unit, residue_management, soil_zone, n_source, n_time,
    n_place, soil_texture, spring_moisture_condition, irrigation_flag, n_source_percent_n, wue, epsilon, nminus1, crop_unit_conv_coef,
    spring_soil_moisture, b0ph, b1ph, b2ph, phmax, phmin, b0ec, b1ec, b0precip, b1precip, soil_zone_id, b0ag, b0bg)
end