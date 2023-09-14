using Git
using BSON
using Distributed
const git = Git.git()
# Function to create AFFIRM folder environments
function create_affirm()
    """
    This function creates all the folders for AFFIRM.jl runs at the user specified folder location. No argument needed
    """
    run(`$git clone https://github.com/4SAnalyticsnModelling/AFFIRM-data`)
    for folder_ in ["data", "input", "output", "src"]
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
# Function to import AFFIRM coefficients from config files in data folder
function get_coefficients(data_file_path :: String = "../data/")
    """
    get_coefficients(data_file_path :: String = "../data/")
    This function gets the coefficients for AFFIRM.jl model runs.
    """
    for item in ["crop_name", "previous_crop", "previous_crop_yld_unit", "residue_management", "soil_zone", "n_source", "n_source_percent_n", "n_time", "n_place", "soil_texture", "spring_moisture_condition", "irrigation_flag", "wue", "epsilon", "nminus1", "crop_unit_conv_coef", "spring_soil_moisture", "b0ph", "b1ph", "b2ph", "phmax", "phmin", "b0ec", "b1ec", "b0precip", "b1precip", "soil_zone_id", "b0ag", "b0bg"]
        file_config = BSON.load(data_file_path * item * ".bson")[Symbol(item)]
        item_name = Symbol(item)
        @eval (const ($item_name) = ($file_config))
    end
end
# Function to run AFFIRM batch script, process the inputs and write the outputs
function run_affirm(input_file_path :: String = "../input/AFFIRM-batch-inputs.csv", output_file_path :: String = "../output/AFFIRM-batch-outputs.csv", output_dir :: String = "../output/")
    """
    run_affirm(input_file_path :: String = "../input/AFFIRM-batch-inputs.csv", output_file_path :: String = "../output/AFFIRM-batch-outputs.csv", output_dir :: String = "../output/")
    This function executes the AFFIRM.jl model runs.
    """
    get_coefficients()
    try
        map(file_ -> rm(output_dir * file_, force = true), readdir(output_dir))
    catch
    end
    input_f = open(input_file_path, "r")
    line_ = readline(input_f)
    try
        while !eof(input_f)
            line_ = readline(input_f)
            line = split(line_, ",")
            inp_township, inp_range = parse.(Int, line[1:2])
            inp_meridian = parse(Int, line[3][2])
            inp_meridian_ = line[3]
            inp_som = get_distribution(line[4])
            inp_soil_texture, inp_spring_soil_moisture = get_combined_simulation.(line[5:6])
            inp_soil_ph, inp_soil_ec = get_distribution.(line[7:8])
            inp_current_crop, inp_irrig_flg = parse.(Int, line[9:10])
            inp_my_precip = line[11] == "" || line[11] == "-" || line[11] == "na" ? 0.0f0 : get_distribution(line[11])
            inp_my_irrig = line[12] == "" || line[12] == "-" || line[12] == "na" ? 0.0f0 : get_distribution(line[12])
            inp_n_source = parse(Int, line[13])
            inp_n_time, inp_n_place = get_combined_simulation.(line[14:15])
            inp_soil_test_n = get_distribution(line[16])
            inp_prev_crop = parse(Int, line[17]) 
            inp_prev_crop_yld = get_distribution(line[18])
            inp_prev_crop_yld_unit = parse(Int, line[19]) 
            inp_res_mgmt_flg = get_combined_simulation(line[20])
            inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio = get_distribution.(line[21:24])   
    
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
            conv_yld_kgha_buac = crop_unit_conv_coef[inp_current_crop] * 0.000405f0
    
            write_f = open(output_file_path, "a")
            write(write_f, writeoutputs("Township", "Range", "Meridian", "Soil Zone", "Soil organic matter (0-6\") (%)", "Soil texture, Spring soil moisture", "Soil pH (0-6\" or 0-12\")", "Soil EC (0-6\" or 0-12\") (meq/100g)", "Crop", "Irrigation", "Growing season moisture flag", "Growing season precipitation (May-Aug) + irrigation (if any) (mm)", "Nitrogen fertilizer product", "Nitrogen fertilizer application timing", "Nitrogen fertilizer application placement", "Soil test nitrogen (0-24\") (lb N/ac)", "Previous crop", "Previous crop yield", "Previous crop yield unit", "Residue management", "Crop available nitrogen from applied manure (lb N/ac)", "Expected crop price (\$/bu)", "Fertilizer price (\$/tonne)", "User chosen investment ratio", "Estimated N release from N mineralization over the growing season (lb N/ac)", "N credit from previous crop residue (lb N/ac)", "Total plant available nitrogen from soil (lb N/ac)", "Fertilizer N application rate (lb N/ac)", "Predicted crop yield (bu/ac)", "Predicted yield increase (bu/ac)", "Added yield increase (bu/ac)", "Estimated revenue from fertilizer N (\$/ac)", "Marginal return or Gross margin change (\$/ac)", "Total cost of fertilizer N (\$/ac)", "Marginal cost of fertilizer N (\$/ac)", "Estimated Investment Ratio", "Recommended?" ))
            close(write_f)
    
            for inp_soil_texture0 in Distributed.splitrange(1, length(inp_soil_texture), Threads.nthreads())
                soil_texture_par = Threads.@spawn begin 
                    for inp_soil_texture in inp_soil_texture[inp_soil_texture0]
                        for inp_spring_soil_moisture0 in Distributed.splitrange(1, length(inp_spring_soil_moisture), Threads.nthreads())
                            spring_moisture_par = Threads.@spawn begin
                                for inp_spring_soil_moisture in inp_spring_soil_moisture[inp_spring_soil_moisture0]
                                    spring_soil_water_content = spring_soil_moisture[inp_spring_soil_moisture, inp_soil_texture]
                                    growing_season_precip_list = growing_season_precip
                                    for growing_season_precip0 in Distributed.splitrange(1, length(growing_season_precip), Threads.nthreads())
                                        growing_season_precip_par = Threads.@spawn begin       
                                            for growing_season_precip in growing_season_precip_list[growing_season_precip0]
                                                total_plant_available_moisture = growing_season_precip + spring_soil_water_content
                                                growing_season_precip_id = findall(q -> q == growing_season_precip, growing_season_precip_list)
                                                growing_season_precip_flag_ = growing_season_precip_flag[growing_season_precip_id]
                                                for inp_soil_ph0 in Distributed.splitrange(1, length(inp_soil_ph), Threads.nthreads())
                                                    soil_ph_par = Threads.@spawn begin
                                                        for inp_soil_ph in inp_soil_ph[inp_soil_ph0]
                                                            soil_ph = min(phmax[inp_current_crop], max(phmin[inp_current_crop], inp_soil_ph))
                                                            ph_adjust = max(0.0f0, min(1.0f0, b0ph[inp_current_crop] + b1ph[inp_current_crop] * soil_ph + b2ph[inp_current_crop] * soil_ph ^ 2.0f0))
                                                            for inp_soil_ec0 in Distributed.splitrange(1, length(inp_soil_ec), Threads.nthreads())
                                                                soil_ec_par = Threads.@spawn begin
                                                                    for inp_soil_ec in inp_soil_ec[inp_soil_ec0]
                                                                        ec_adjust = max(0.0f0, min(1.0f0, b0ec[inp_current_crop] + b1ec[inp_current_crop] * inp_soil_ec))
                                                                        for inp_som0 in Distributed.splitrange(1, length(inp_som), Threads.nthreads())
                                                                            som_par = Threads.@spawn begin
                                                                                for inp_som in inp_som[inp_som0]
                                                                                    enr = round((20.6f0 + 13.2f0 * inp_som - 0.1777f0 * inp_som ^ 2.0f0) / kg_ha_n_lb_ac, digits = 0)
                                                                                    for inp_prev_crop_yld0 in Distributed.splitrange(1, length(inp_prev_crop_yld), Threads.nthreads())
                                                                                        prev_crop_yld_par = Threads.@spawn begin
                                                                                            for inp_prev_crop_yld in inp_prev_crop_yld[inp_prev_crop_yld0]
                                                                                                for inp_res_mgmt_flg0 in Distributed.splitrange(1, length(inp_res_mgmt_flg), Threads.nthreads())
                                                                                                    res_mgmt_flg_par = Threads.@spawn begin
                                                                                                        for inp_res_mgmt_flg in inp_res_mgmt_flg[inp_res_mgmt_flg0]
                                                                                                            residue_n_credit = round(inp_prev_crop_yld * (residue_management_multiplyer[inp_res_mgmt_flg] * b0ag[inp_prev_crop, inp_prev_crop_yld_unit] + b0bg[inp_prev_crop, inp_prev_crop_yld_unit]), digits = 0)
                                                                                                            for inp_soil_test_n0 in Distributed.splitrange(1, length(inp_soil_test_n), Threads.nthreads())
                                                                                                                soil_test_n_par = Threads.@spawn begin
                                                                                                                    for inp_soil_test_n in inp_soil_test_n[inp_soil_test_n0]
                                                                                                                        for inp_manure_n0 in Distributed.splitrange(1, length(inp_manure_n), Threads.nthreads())
                                                                                                                            manure_n_par = Threads.@spawn begin
                                                                                                                                for inp_manure_n in inp_manure_n[inp_manure_n0]
                                                                                                                                    plant_available_soil_n = round(enr + residue_n_credit + inp_soil_test_n + inp_manure_n, digits = 0)
                                                                                                                                    for inp_n_time0 in Distributed.splitrange(1, length(inp_n_time), Threads.nthreads())
                                                                                                                                        n_time_par = Threads.@spawn begin
                                                                                                                                            for inp_n_time in inp_n_time[inp_n_time0]
                                                                                                                                                for inp_n_place0 in Distributed.splitrange(1, length(inp_n_place), Threads.nthreads())
                                                                                                                                                    n_place_par = Threads.@spawn begin
                                                                                                                                                        for inp_n_place in inp_n_place[inp_n_place0]
                                                                                                                                                            for inp_crop_price0 in Distributed.splitrange(1, length(inp_crop_price), Threads.nthreads())
                                                                                                                                                                crop_price_par = Threads.@spawn begin
                                                                                                                                                                    for inp_crop_price in inp_crop_price[inp_crop_price0]
                                                                                                                                                                        for inp_fertilizer_price0 in Distributed.splitrange(1, length(inp_fertilizer_price), Threads.nthreads())
                                                                                                                                                                            fertilizer_price_par = Threads.@spawn begin
                                                                                                                                                                                for inp_fertilizer_price in inp_fertilizer_price[inp_fertilizer_price0]
                                                                                                                                                                                    for inp_investment_ratio0 in Distributed.splitrange(1, length(inp_investment_ratio), Threads.nthreads())
                                                                                                                                                                                        investment_ratio_par = Threads.@spawn begin
                                                                                                                                                                                            for inp_investment_ratio in inp_investment_ratio[inp_investment_ratio0]
                                                                                                                                                                                                high_n_rate = ns_max - plant_available_soil_n
                                                                                                                                                                                                if rem(high_n_rate, n_rate_step_size) > 0.0f0
                                                                                                                                                                                                    high_n_rate = n_rate_step_size * (1.0f0 + div(high_n_rate, 10.0f0))
                                                                                                                                                                                                end
                                                                                                                                                                                                predicted_crop_yield_previous_n_rate = 0
                                                                                                                                                                                                estimated_revenue_from_fertilizer_n_previous_n_rate = 0.0f0
                                                                                                                                                                                                predicted_yield_increase = 0
                                                                                                                                                                                                marginal_return = 0
                                                                                                                                                                                                total_cost_of_fertilizer_n_previous_n_rate = 0.0f0
                                                                                                                                                                                                estimated_investment_ratio_previous_n_rate = 0.0f0
                                                                                                                                                                                                recommend_flag = ""
                                                                                                                                                                                                for n_rate in low_n_rate:n_rate_step_size:high_n_rate
                                                                                                                                                                                                    plant_available_total_n = plant_available_soil_n + n_rate
                                                                                                                                                                                                    predicted_crop_yield = round(ph_adjust * ec_adjust * conv_yld_kgha_buac * wue[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop] * total_plant_available_moisture * (1.0f0 - 10.0f0 ^ (-epsilon[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop] * plant_available_total_n * kg_ha_n_lb_ac * (wue[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop] * total_plant_available_moisture) ^ nminus1[inp_n_place, inp_n_time, inp_n_source, soil_zone_id[inp_township, inp_range, inp_meridian], inp_current_crop])), digits = 2)
                                                                                                                                                                                                    if n_rate == low_n_rate 
                                                                                                                                                                                                        predicted_yield_increase = 0
                                                                                                                                                                                                        added_yield_increase = 0
                                                                                                                                                                                                    else
                                                                                                                                                                                                        added_yield_increase = round(predicted_crop_yield - predicted_crop_yield_previous_n_rate, digits = 2)
                                                                                                                                                                                                        predicted_yield_increase += round(added_yield_increase, digits = 2)
                                                                                                                                                                                                    end
                                                                                                                                                                                                    estimated_revenue_from_fertilizer_n = round(predicted_yield_increase * inp_crop_price, digits = 2)
                                                                                                                                                                                                    if n_rate == low_n_rate
                                                                                                                                                                                                        marginal_return = 0
                                                                                                                                                                                                    else
                                                                                                                                                                                                        marginal_return = round(estimated_revenue_from_fertilizer_n - estimated_revenue_from_fertilizer_n_previous_n_rate, digits = 1)
                                                                                                                                                                                                    end
                                                                                                                                                                                                    total_cost_of_fertilizer_n = round(n_rate * inp_fertilizer_price / 1000.0f0 * 0.4536f0 / (n_source_percent_n[inp_n_source] / 100.0f0), digits = 2)
                                                                                                                                                                                                    marginal_cost_of_fertilizer_n = round(total_cost_of_fertilizer_n - total_cost_of_fertilizer_n_previous_n_rate, digits = 2)
                                                                                                                                                                                                    if n_rate > low_n_rate
                                                                                                                                                                                                        estimated_investment_ratio = round(marginal_return / marginal_cost_of_fertilizer_n, digits = 2)
                                                                                                                                                                                                    else
                                                                                                                                                                                                        estimated_investment_ratio = 0.0f0
                                                                                                                                                                                                    end
                                                                                                                                                                                                    if n_rate > low_n_rate + n_rate_step_size
                                                                                                                                                                                                        if estimated_investment_ratio < inp_investment_ratio && estimated_investment_ratio_previous_n_rate >= inp_investment_ratio
                                                                                                                                                                                                            recommend_flag = "Yes"
                                                                                                                                                                                                        else
                                                                                                                                                                                                            recommend_flag = ""
                                                                                                                                                                                                        end
                                                                                                                                                                                                    else
                                                                                                                                                                                                        recommend_flag = ""
                                                                                                                                                                                                    end
                                                                                                                                                                                                    predicted_crop_yield_previous_n_rate = predicted_crop_yield
                                                                                                                                                                                                    estimated_revenue_from_fertilizer_n_previous_n_rate = estimated_revenue_from_fertilizer_n
                                                                                                                                                                                                    total_cost_of_fertilizer_n_previous_n_rate = total_cost_of_fertilizer_n
                                                                                                                                                                                                    estimated_investment_ratio_previous_n_rate = estimated_investment_ratio
                                                                                                                                                                                                    predicted_crop_yield = round(predicted_crop_yield, digits = 1)
                                                                                                                                                                                                    predicted_yield_increase = round(predicted_yield_increase, digits = 1)
                                                                                                                                                                                                    added_yield_increase = round(added_yield_increase, digits = 1)
                                                                                                                                                                                                    estimated_revenue_from_fertilizer_n = round(estimated_revenue_from_fertilizer_n, digits = 1)
                                                                                                                                                                                                    marginal_return = round(marginal_return, digits = 1)
                                                                                                                                                                                                    total_cost_of_fertilizer_n = round(total_cost_of_fertilizer_n, digits = 1)
                                                                                                                                                                                                    marginal_cost_of_fertilizer_n = round(marginal_cost_of_fertilizer_n, digits = 1)
                                                                                                                                                                                                    estimated_investment_ratio = round(estimated_investment_ratio, digits = 1)
                                                                                                                                                                                                    write_f = open(output_file_path, "a")
                                                                                                                                                                                                    if added_yield_increase >= 0.5f0
                                                                                                                                                                                                        if recommend_flag == "Yes"
                                                                                                                                                                                                            if n_rate > low_n_rate
                                                                                                                                                                                                                write(write_f, writeoutputs(inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate, predicted_crop_yield, predicted_yield_increase, added_yield_increase, estimated_revenue_from_fertilizer_n, marginal_return, total_cost_of_fertilizer_n, marginal_cost_of_fertilizer_n, estimated_investment_ratio, recommend_flag))
                                                                                                                                                                                                            else
                                                                                                                                                                                                                write(write_f, writeoutputs(inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate, predicted_crop_yield))
                                                                                                                                                                                                            end
                                                                                                                                                                                                        else
                                                                                                                                                                                                            if n_rate > low_n_rate
                                                                                                                                                                                                                write(write_f, writeoutputs(inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate, predicted_crop_yield, predicted_yield_increase, added_yield_increase, estimated_revenue_from_fertilizer_n, marginal_return, total_cost_of_fertilizer_n, marginal_cost_of_fertilizer_n, estimated_investment_ratio))
                                                                                                                                                                                                            else
                                                                                                                                                                                                                write(write_f, writeoutputs(inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate, predicted_crop_yield))
                                                                                                                                                                                                            end
                                                                                                                                                                                                        end
                                                                                                                                                                                                    else
                                                                                                                                                                                                        if n_rate == low_n_rate
                                                                                                                                                                                                            write(write_f, writeoutputs(inp_township, inp_range, inp_meridian_, soil_zone_, inp_som, soil_texture[inp_soil_texture], spring_moisture_condition[inp_spring_soil_moisture], inp_soil_ph, inp_soil_ec, crop_name[inp_current_crop], irrigation_flag[inp_irrig_flg], growing_season_precip_flag_, growing_season_precip, n_source[inp_n_source], n_time[inp_n_time], n_place[inp_n_place], inp_soil_test_n, previous_crop[inp_prev_crop], inp_prev_crop_yld, previous_crop_yld_unit[inp_prev_crop_yld_unit], residue_management[inp_res_mgmt_flg], inp_manure_n, inp_crop_price, inp_fertilizer_price, inp_investment_ratio, enr, residue_n_credit, plant_available_soil_n, n_rate, predicted_crop_yield))
                                                                                                                                                                                                        end
                                                                                                                                                                                                    end
                                                                                                                                                                                                    close(write_f)
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                        end
                                                                                                                                                                                        fetch(investment_ratio_par)
                                                                                                                                                                                    end
                                                                                                                                                                                end
                                                                                                                                                                            end
                                                                                                                                                                            fetch(fertilizer_price_par)
                                                                                                                                                                        end
                                                                                                                                                                    end
                                                                                                                                                                end
                                                                                                                                                                fetch(crop_price_par)
                                                                                                                                                            end
                                                                                                                                                        end
                                                                                                                                                    end
                                                                                                                                                    fetch(n_place_par)
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                        end
                                                                                                                                        fetch(n_time_par)
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end
                                                                                                                            fetch(manure_n_par)
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                                fetch(soil_test_n_par)
                                                                                                            end
                                                                                                        end
                                                                                                    end
                                                                                                    fetch(res_mgmt_flg_par)
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                        fetch(prev_crop_yld_par)
                                                                                    end
                                                                                end
                                                                            end
                                                                            fetch(som_par)
                                                                        end
                                                                    end
                                                                end
                                                                fetch(soil_ec_par)
                                                            end
                                                        end
                                                    end
                                                    fetch(soil_ph_par)
                                                end
                                            end
                                        end
                                        fetch(growing_season_precip_par)
                                    end
                                end
                            end
                            fetch(spring_moisture_par)
                        end
                    end
                end
                fetch(soil_texture_par)
            end
        end
        println("Congratulations!! AFFIRM.jl model run has completed successfully! Please check out the model outputs in the output folder. Thank you for using AFFIRM.jl!!")
    catch
        println("Oops...something went wrong!! Please check your input file to make sure the input formats are correct. With the right inputs please try again! Good luck for your next run!")
    end
    sleep(10)        
end
