"""
```julia
create_affirm()
```
This function creates all the folders for AFFIRM.jl runs at the user specified folder location. No argument needed
"""
function create_affirm()
    try
        rm("AFFIRM-data", force = true, recursive = true)
    catch
    end
    run(`$git clone https://github.com/4SAnalyticsnModelling/AFFIRM-data`)
    for folder_ in ["data", "input", "output", "src"]
        try
            rm(folder_, force = true, recursive = true)
        catch
        end
        mkdir(folder_)
        if folder_ != "output"
            folder_to_copy_from = "AFFIRM-data/" * folder_
            for file_ in readdir(folder_to_copy_from)
                mv(folder_to_copy_from * "/" * file_, folder_ * "/" * file_)
            end
        end
    end
    rm("AFFIRM-data", force = true, recursive = true)
end
"""
```julia
get_coefficients(data_file_path :: String = "../data/")
```
This function gets the coefficients for AFFIRM.jl model runs.
"""
function get_coefficients(data_file_path :: String = "../data/")
    for item in ["crop_name", "previous_crop", "previous_crop_yld_unit", "residue_management", "soil_zone", "n_source", "n_source_percent_n", "n_time", "n_place", "soil_texture", "spring_moisture_condition", "irrigation_flag", "wue", "epsilon", "nminus1", "crop_unit_conv_coef", "spring_soil_moisture", "b0ph", "b1ph", "b2ph", "phmax", "phmin", "b0ec", "b1ec", "b0precip", "b1precip", "soil_zone_id", "b0ag", "b0bg"]
        file_config = BSON.load(data_file_path * item * ".bson")[Symbol(item)]
        if length(file_config) <= 100 && !occursin("OffsetArray", string(typeof(file_config)))
            file_config = SArray{Tuple{size(file_config)...}, eltype(file_config)}(file_config)
        end
        item_name = Symbol(item)
        @eval (($item_name) = ($file_config))
    end
end
"""
```julia
run_affirm(input_file_path :: String = "../input/AFFIRM-batch-inputs.csv", output_file_path :: String = "../output/AFFIRM-batch-outputs.csv", output_dir :: String = "../output/")
```
This function executes the AFFIRM.jl model runs. Key equations in AFFIRM.jl model was built on the research from the following publications:
Yield response modelling:https://www.sciencedirect.com/science/article/pii/S1573521400800160 
                         https://www.sciencedirect.com/science/article/pii/S1573521400800172
Nitrogen mineralization modelling: https://cdnsciencepub.com/doi/10.4141/cjss84-035 
Previous residue N-credit modelling: https://link.springer.com/article/10.1023/A:1025195826663
"""
function run_affirm(input_file_path :: String = "../input/AFFIRM-batch-inputs.csv", output_file_path :: String = "../output/AFFIRM-batch-outputs.csv", output_dir :: String = "../output/")
    get_coefficients()
    line_count_ = 0
    line_count = 0
    try
        map(file_ -> rm(output_dir * file_, force = true), readdir(output_dir))
    catch
    end

    Base.printstyled("\nPreparing ", color = :green, bold = true)
    print("AFFIRM.jl batch runs...\n")

    if line_count_ == 0
        input_f_ = open(input_file_path, "r")
        while !eof(input_f_)
            readline(input_f_)
            line_count_ += 1
        end
        close(input_f_)
        line_count_ = line_count_ - 1
    end

    write_f = open(output_file_path, "w")
    println(write_f, writeoutputs("Index", "Township", "Range", "Meridian", "Soil Zone", "Soil organic matter (0-6\") (%)", "Soil texture", "Spring soil moisture", "Soil pH (0-6\" or 0-12\")", "Soil EC (0-6\" or 0-12\") (mS/cm)", "Crop", "Irrigation", "Growing season moisture flag", "Growing season precipitation (May-Aug) + irrigation (if any) (mm)", "Nitrogen fertilizer product", "Nitrogen fertilizer application timing", "Nitrogen fertilizer application placement", "Soil test nitrogen (0-24\") (lb N/ac)", "Previous crop", "Previous crop yield", "Previous crop yield unit", "Residue management", "Crop available nitrogen from applied manure (lb N/ac)", "Expected crop price (\$/bu)", "Fertilizer price (\$/tonne)", "User chosen investment ratio", "Estimated N release from N mineralization over the growing season (lb N/ac)", "N credit from previous crop residue (lb N/ac)", "Total plant available nitrogen from soil (lb N/ac)", "Fertilizer N application rate (lb N/ac)", "Predicted crop yield (bu/ac)", "Predicted yield increase (bu/ac)", "Added yield increase (bu/ac)", "Estimated revenue from fertilizer N (\$/ac)", "Marginal return or Gross margin change (\$/ac)", "Total cost of fertilizer N (\$/ac)", "Marginal cost of fertilizer N (\$/ac)", "Estimated Investment Ratio", "Recommended?"))
    close(write_f)

    logf_ = open(output_dir * "AFFIRM-batch-logfile", "w")
    close(logf_)

    input_f = open(input_file_path, "r")
    line_ = readline(input_f)

    Base.printstyled("\nProcessing ", color = :green, bold = true)
    print("AFFIRM.jl batch runs...\n")

    write_f = open(output_file_path, "a")
    logf_ = open(output_dir * "AFFIRM-batch-logfile", "w")

    while !eof(input_f)
        line_count += 1
        progress_meter_ = round(100.0f0 * line_count / line_count_, digits = 2)
        progress_meter = string(progress_meter_)
        progress_meter_int = lpad(split(progress_meter, ".")[1], 3, "0")
        progress_meter_dec = rpad(split(progress_meter, ".")[2], 2, "0")
        progress_meter = progress_meter_int * "." * progress_meter_dec * "%"
        Base.printstyled("Info: ", color = :blue, bold = true)
        print("Progressing $progress_meter\u001b[1000D")   
        try
            line_ = readline(input_f)
            line = split(line_, ",")
            inp_index, inp_township, inp_range = parse.(Int, line[1:3])
            inp_meridian = parse(Int, line[4][2])
            inp_meridian_ = line[4]
            inp_som = get_distribution(line[5])
            inp_soil_texture, inp_spring_soil_moisture = get_combined_simulation.(line[6:7])
            inp_soil_ph, inp_soil_ec = get_distribution.(line[8:9])
            inp_current_crop, inp_irrig_flg = parse.(Int, line[10:11])
            inp_my_precip = line[12] == "" || line[12] == "-" || line[12] == "na" ? 0.0f0 : get_distribution(line[12])
            inp_my_irrig = line[13] == "" || line[13] == "-" || line[13] == "na" ? 0.0f0 : get_distribution(line[13])
            inp_n_source = parse(Int, line[14])
            inp_n_time, inp_n_place = get_combined_simulation.(line[15:16])
            inp_soil_test_n = get_distribution(line[17])
            inp_prev_crop = parse(Int, line[18]) 
            inp_prev_crop_yld = get_distribution(line[19])
            inp_prev_crop_yld_unit = parse(Int, line[20]) 
            inp_res_mgmt_flg = get_combined_simulation(line[21])
            inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio = get_distribution.(line[22:25])   

            if inp_township in axes(soil_zone_id, 1) && inp_range in axes(soil_zone_id, 2) && inp_meridian in axes(soil_zone_id, 3) && soil_zone_id[inp_township, inp_range, inp_meridian] > 0
                soil_zone_ = soil_zone[soil_zone_id[inp_township, inp_range, inp_meridian]]

                if inp_irrig_flg == 1 # no irrigation / dryland
                    growing_season_precip_optimum_moisture = round(b0precip[inp_township, inp_range, inp_meridian] - 10.0f0 * b1precip[inp_township, inp_range, inp_meridian], digits = 0)
                    growing_season_precip_intermediate_moisture = round(b0precip[inp_township, inp_range, inp_meridian] - 50.0f0 * b1precip[inp_township, inp_range, inp_meridian], digits = 0)
                    growing_season_precip_low_moisture = round(b0precip[inp_township, inp_range, inp_meridian] - 90.0f0 * b1precip[inp_township, inp_range, inp_meridian], digits = 0)
                    if sum(inp_my_precip) > 0.0f0
                        growing_season_precip = vcat(inp_my_precip, [growing_season_precip_low_moisture, growing_season_precip_intermediate_moisture, growing_season_precip_optimum_moisture])
                        growing_season_precip_flag = vcat(["Growing season precipitation - User input" for _ in 1:length(inp_my_precip)], ["Growing season precipitation - Low moisture condition", "Growing season precipitation - Intermediate moisture condition", "Growing season precipitation - Optimum moisture condition"])
                    else
                        growing_season_precip = [growing_season_precip_low_moisture, growing_season_precip_intermediate_moisture, growing_season_precip_optimum_moisture]
                        growing_season_precip_flag = ["Growing season precipitation - Low moisture condition", "Growing season precipitation - Intermediate moisture condition", "Growing season precipitation - Optimum moisture condition"]
                    end
                else  # irrigation
                    growing_season_precip_optimum_moisture = round(b0irrig - 10.0f0 * b1irrig, digits = 0)
                    growing_season_precip_intermediate_moisture = round(b0irrig - 50.0f0 * b1irrig, digits = 0)
                    growing_season_precip_low_moisture = round(b0irrig - 90.0f0 * b1irrig, digits = 0)
                    if sum(inp_my_irrig) > 0.0f0
                        growing_season_precip = vcat(inp_my_irrig, [growing_season_precip_low_moisture, growing_season_precip_intermediate_moisture, growing_season_precip_optimum_moisture])
                        growing_season_precip_flag = vcat(["Irrigation water amount - User input" for _ in 1:length(inp_my_irrig)], ["Irrigation water amount - Low moisture condition", "Irrigation water amount - Intermediate moisture condition", "Irrigation water amount - Optimum moisture condition"])
                    else
                        growing_season_precip = [growing_season_precip_low_moisture, growing_season_precip_intermediate_moisture, growing_season_precip_optimum_moisture]
                        growing_season_precip_flag = ["Irrigation water amount - Low moisture condition", "Irrigation water amount - Intermediate moisture condition", "Irrigation water amount - Optimum moisture condition"]
                    end
                end

                growing_season_precip = SArray{Tuple{length(growing_season_precip),}, eltype(growing_season_precip)}(growing_season_precip)
                growing_season_precip_flag = SArray{Tuple{length(growing_season_precip_flag),}, eltype(growing_season_precip_flag)}(growing_season_precip_flag)

                conv_yld_kgha_buac = crop_unit_conv_coef[inp_current_crop] * 0.000405f0

                for inp_soil_texture in inp_soil_texture, inp_spring_soil_moisture in inp_spring_soil_moisture
                    spring_soil_water_content = spring_soil_moisture[inp_spring_soil_moisture, inp_soil_texture]
                    growing_season_precip_list = growing_season_precip
                    for growing_season_precip in growing_season_precip
                        total_plant_available_moisture = growing_season_precip + spring_soil_water_content
                        growing_season_precip_id = findall(q -> q == growing_season_precip, growing_season_precip_list)
                        growing_season_precip_flag_ = growing_season_precip_flag[growing_season_precip_id]
                        for inp_soil_ph in inp_soil_ph
                            soil_ph = min(phmax[inp_current_crop], max(phmin[inp_current_crop], inp_soil_ph))
                            ph_adjust = max(0.0f0, min(1.0f0, b0ph[inp_current_crop] + b1ph[inp_current_crop] * soil_ph + b2ph[inp_current_crop] * soil_ph ^ 2.0f0))
                            # pH warning message when crop yield starts being affected by soil pH
                            if ph_adjust < 0.9f0
                                row_number = line_count + 1
                                msg = "Crop yield is severely affected by adverse soil pH. This warning is for the scenario at row $row_number of your input file."
                                println(logf_, writeoutputs("Warning: $msg"))
                                Base.printstyled("\nWarning: ", color = :yellow, bold = true)
                                print(msg, "\n")
                            elseif ph_adjust < 0.75f0
                                row_number = line_count + 1
                                msg = "Crop yield is severely affected by adverse soil pH. This warning is for the scenario at row $row_number of your input file."
                                println(logf_, writeoutputs("Warning: $msg"))
                                Base.printstyled("\nWarning: ", color = :yellow, bold = true)
                                print(msg, "\n")
                            end
                            for inp_soil_ec in inp_soil_ec
                                ec_adjust = max(0.0f0, min(1.0f0, b0ec[inp_current_crop] + b1ec[inp_current_crop] * inp_soil_ec))
                                # Salinity warning message when crop yield starts being affected by soil pH
                                if ec_adjust < 0.9f0
                                    row_number = line_count + 1
                                    msg = "Crop yield is moderately affected by soil salinity. This warning is for the scenario at row $row_number of your input file."
                                    println(logf_, writeoutputs("Warning: $msg"))
                                    Base.printstyled("\nWarning: ", color = :yellow, bold = true)
                                    print(msg, "\n")
                                elseif ec_adjust < 0.75f0
                                    row_number = line_count + 1
                                    msg = "Crop yield is severely affected by soil salinity. This warning is for the scenario at row $row_number of your input file."
                                    println(logf_, writeoutputs("Warning: $msg"))
                                    Base.printstyled("\nWarning: ", color = :yellow, bold = true)
                                    print(msg, "\n")
                                end
                                for inp_som in inp_som
                                    enr = round((20.6f0 + 13.2f0 * inp_som - 0.1777f0 * inp_som ^ 2.0f0) / kg_ha_n_lb_ac, digits = 1)
                                    for inp_prev_crop_yld in inp_prev_crop_yld, inp_res_mgmt_flg in inp_res_mgmt_flg
                                        residue_n_credit = round(inp_prev_crop_yld * (residue_management_multiplyer[inp_res_mgmt_flg] * b0ag[inp_prev_crop, inp_prev_crop_yld_unit] + b0bg[inp_prev_crop, inp_prev_crop_yld_unit]), digits = 0)
                                        for inp_soil_test_n in inp_soil_test_n, inp_manure_n in inp_manure_n
                                            plant_available_soil_n = round(enr + residue_n_credit + inp_soil_test_n + inp_manure_n, digits = 0)
                                            for inp_n_time in inp_n_time, inp_n_place in inp_n_place, inp_crop_price in inp_crop_price, inp_fertilizer_price in inp_fertilizer_price, inp_investment_ratio in inp_investment_ratio
                                                WUE = wue[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop]
                                                ϵ = epsilon[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop]
                                                NMINUS1 = nminus1[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop]
                                                if WUE == 0.0f0 && ϵ == 0.0f0 && NMINUS1 == 0.0f0
                                                    set_calculation_flag = 1
                                                    comments_ = "Yield response information is not available for either the legal land location or the combination of fertilizer management that you have chosen for your field in this scenario; please try with a different combination"
                                                    println(write_f, writeoutputs(inp_index, inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, comments_))
                                                    continue
                                                else
                                                    set_calculation_flag = 0                                                    
                                                    high_n_rate = ns_max - plant_available_soil_n
                                                    if high_n_rate % n_rate_step_size > 0.0f0
                                                        high_n_rate = n_rate_step_size * (1.0f0 + high_n_rate ÷ 10.0f0)
                                                    end
                                                    n_rate_list = low_n_rate:n_rate_step_size:high_n_rate
                                                    n_rate_list = SArray{Tuple{size(n_rate_list)...}, eltype(n_rate_list)}(n_rate_list)
                                                    # initialize arrays for subsequent calculations
                                                    max_n_rate_id = length(n_rate_list)
                                                    predicted_crop_yield = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    added_yield_increase = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    predicted_yield_increase = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    total_cost_of_fertilizer_n = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    marginal_cost_of_fertilizer_n = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    estimated_revenue_from_fertilizer_n = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    marginal_return = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                    estimated_investment_ratio = MArray{Tuple{max_n_rate_id,}, Float32}(zeros(max_n_rate_id))
                                                end
                                                if set_calculation_flag == 0
                                                    @turbo thread = 20 for i in eachindex(n_rate_list)
                                                        n_rate = n_rate_list[i]
                                                        n_rate_id = 1 + floor(Int64, n_rate / n_rate_step_size)
                                                        plant_available_total_n = plant_available_soil_n + n_rate
                                                        predicted_crop_yield[n_rate_id] = ph_adjust * ec_adjust * conv_yld_kgha_buac * WUE * total_plant_available_moisture * (1.0f0 - 10.0f0 ^ (- ϵ * plant_available_total_n * kg_ha_n_lb_ac * (WUE * total_plant_available_moisture) ^ NMINUS1))
                                                        total_cost_of_fertilizer_n[n_rate_id] = n_rate * inp_fertilizer_price / 1000.0f0 * 0.4536f0 / (n_source_percent_n[inp_n_source] / 100.0f0)
                                                    end
                                                    @turbo thread = 20 for i in 2:length(n_rate_list)
                                                        n_rate = n_rate_list[i]
                                                        n_rate_id = floor(Int64, n_rate / n_rate_step_size)
                                                        n_rate_id_prev = n_rate_id - 1
                                                        added_yield_increase[n_rate_id] = predicted_crop_yield[n_rate_id] - predicted_crop_yield[n_rate_id_prev]
                                                        marginal_cost_of_fertilizer_n[n_rate_id] = total_cost_of_fertilizer_n[n_rate_id] - total_cost_of_fertilizer_n[n_rate_id_prev]
                                                    end
                                                    @fastmath @inbounds @simd for i in 2:length(n_rate_list)
                                                        n_rate = n_rate_list[i]
                                                        n_rate_id = floor(Int64, n_rate / n_rate_step_size)
                                                        sum_ids = n_rate_id:-1:2
                                                        predicted_yield_increase[n_rate_id] = sum(added_yield_increase[sum_ids])
                                                    end
                                                    @turbo thread = 20 for i in 2:length(n_rate_list)
                                                        n_rate = n_rate_list[i]
                                                        n_rate_id = floor(Int64, n_rate / n_rate_step_size)
                                                        n_rate_id_prev = n_rate_id - 1
                                                        estimated_revenue_from_fertilizer_n[n_rate_id] = predicted_yield_increase[n_rate_id] * inp_crop_price
                                                        marginal_return[n_rate_id] = (predicted_yield_increase[n_rate_id] - predicted_yield_increase[n_rate_id_prev]) * inp_crop_price
                                                        estimated_investment_ratio[n_rate_id] = marginal_return[n_rate_id] / marginal_cost_of_fertilizer_n[n_rate_id]
                                                    end
                                                    @fastmath @inbounds @simd for i in eachindex(n_rate_list)
                                                        predicted_crop_yield[i] = round(predicted_crop_yield[i], digits = 1)
                                                        added_yield_increase[i] = round(added_yield_increase[i], digits = 1)
                                                        predicted_yield_increase[i] = round(predicted_yield_increase[i], digits = 1)
                                                        total_cost_of_fertilizer_n[i] = round(total_cost_of_fertilizer_n[i], digits = 2)
                                                        marginal_cost_of_fertilizer_n[i] = round(marginal_cost_of_fertilizer_n[i], digits = 2)
                                                        estimated_revenue_from_fertilizer_n[i] = round(estimated_revenue_from_fertilizer_n[i], digits = 2)
                                                        marginal_return[i] = round(marginal_return[i], digits = 2)
                                                        estimated_investment_ratio[i] = round(estimated_investment_ratio[i], digits = 1)
                                                    end
                                                end
                                                recommend_flag_set = 0
                                                recommend_flag = ""
                                                if set_calculation_flag == 0 
                                                    for i in eachindex(n_rate_list)
                                                        if recommend_flag_set == 0 
                                                            if i > 1
                                                                if estimated_investment_ratio[i] <= inp_investment_ratio && estimated_investment_ratio[i - 1] > inp_investment_ratio
                                                                    recommend_flag = "Yes"
                                                                    recommend_flag_set = 1
                                                                elseif estimated_investment_ratio[i] == inp_investment_ratio
                                                                    recommend_flag = "Yes"
                                                                    recommend_flag_set = 1
                                                                else
                                                                    recommend_flag = ""
                                                                end
                                                            else
                                                                recommend_flag = "" 
                                                            end
                                                        else
                                                            recommend_flag = ""
                                                        end
                                                        lock(write_f)
                                                        try
                                                            if added_yield_increase[i] >= 0.5f0
                                                                if i > 1
                                                                    if recommend_flag == "Yes"
                                                                        println(write_f, writeoutputs(inp_index, inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate_list[i], predicted_crop_yield[i], predicted_yield_increase[i], added_yield_increase[i], estimated_revenue_from_fertilizer_n[i], marginal_return[i], total_cost_of_fertilizer_n[i], marginal_cost_of_fertilizer_n[i], estimated_investment_ratio[i], recommend_flag))
                                                                    else
                                                                        println(write_f, writeoutputs(inp_index, inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate_list[i], predicted_crop_yield[i], predicted_yield_increase[i], added_yield_increase[i], estimated_revenue_from_fertilizer_n[i], marginal_return[i], total_cost_of_fertilizer_n[i], marginal_cost_of_fertilizer_n[i], estimated_investment_ratio[i]))
                                                                    end
                                                                else
                                                                    println(write_f, writeoutputs(inp_index, inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate_list[i], predicted_crop_yield[i]))
                                                                end
                                                            end
                                                        finally
                                                            unlock(write_f)
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            else
                row_number = line_count + 1
                msg = "Legal land location is not valid for the scenario at row $row_number of your input file. No output was written for this scenario."
                println(logf_, writeoutputs("Error: $msg"))
                Base.printstyled("\nError: ", color = :red, bold = true)
                print(msg, "\n")   
            end
        catch
            row_number = line_count + 1
            msg = "The inputs are not valid for the scenario at row $row_number of your input file. No output was written for this scenario."
            println(logf_, writeoutputs("Error: $msg"))
            Base.printstyled("\nError: ", color = :red, bold = true)
            print(msg, "\n")
            continue
        end
    end
    msg = "Congratulations!! AFFIRM.jl model run has completed successfully! 
    Please check out the model outputs in the output folder. 
    Also check out the AFFIRM-batch-logfile in the output folder for any important message. 
    Thank you for using AFFIRM.jl!!"
        println(logf_, writeoutputs("Success: $msg"))
    Base.printstyled("\nCompleted ", color = :green, bold = true)
    print("AFFIRM.jl batch runs...")
    Base.printstyled("\nInfo: ", color = :blue, bold = true)
    print(msg, "\n")

    close(input_f)
    close(write_f)
    close(logf_)
    sleep(10)        
end
