![AFFIRM-logo](https://github.com/4SAnalyticsnModelling/AFFIRM.jl/blob/main/resource/affirm_logo.png) 
# AFFIRM.jl 
Alberta Farm Fertilizer Information Recommendation Manager (AFFIRM) has 3 concurrent versions - [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager), [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) and AFFIRM.jl. The [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager) is the production version of AFFIRM. The [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) is the nitrogen sub-model of [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager). The AFFIRM.jl provides batch run utilites for multiple scenarios and provides similar functionalities as in [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/).

## Getting started with AFFIRM.jl
**==> Install Julia Programming Language**
<br>
<br>[Getting started with Julia](https://docs.julialang.org/en/v1/manual/getting-started/)
<br>
<br>**==> Setting up for parallel computing**
<br>
<br>If you want to take advantage of Julia multithreading for parallel computations with AFFIRM.jl, you can start Julia with as many threads as you want:
<br><br>
--> Windows command prompt:
```
C:\Users\[your username]\Local\AppData\Programs\Julia-[version]\bin\julia --threads 10
```
<br>--> MacOS terminal:
```
/Applications/Julia-[version].app/Contents/Resources/julia/bin/julia --threads 10
```
Once Julia REPL is open, check the number of threads by following command:
```julia
julia> Threads.nthreads()
10
```
Alternatively, you can use ```--threads auto``` to start Julia with the number of threads available on the system:
<br><br>
--> Windows command prompt:
```
C:\Users\[your username]\Local\AppData\Programs\Julia-[version]\bin\julia --threads auto
```
<br>--> MacOS terminal:
```
/Applications/Julia-[version].app/Contents/Resources/julia/bin/julia --threads auto
```
**==> Install AFFIRM.jl**
```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/4SAnalyticsnModelling/AFFIRM.jl")
```
or
```julia
julia> ]
(@v1.9) pkg> add "https://github.com/4SAnalyticsnModelling/AFFIRM.jl"
```
**==> Create AFFIRM.jl folders**
```julia
julia> using AFFIRM
julia> cd("path/to/your/project/folder")
julia> create_affirm()
```
## AFFIRM.jl folder structure
```
> data
> input
  - AFFIRM-batch-inputs.csv
> output
> src
  - runAFFIRM.jl
```
***Do not change any folder or file names without understanding the effects of each change. Better leave the folder and file names as they are and only change the values for the input variables as required***
## Understanding AFFIRM.jl input variables
The input variables are provided in the default ```input/AFFIRM-batch-inputs.csv``` file. The input variables in that file with their corresponding data types are listed below; visit [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) to understand these variables better.
```
- Township => Alberta township id; integer values ranging from 1 to 126
- Range => Alberta range id; integer values ranging from 1 to 30
- Meridian => Alberta meridian id; categorical variables: "W4", "W5" and "W6"
- Soil organic matter (0-6")(%) => numerical variables in decimals
- Soil texture => categorical variables represented by identifiers described below
- Spring soil moisture => categorical variables represented by identifiers described below
- Soil pH (0-6" or 0-12") => numerical variables in decimals
- Soil EC (0-6" or 0-12")(meq/100g) => numerical variables in decimals
- Crop => categorical variables represented by identifiers described below
- Irrigation => categorical variables represented by identifiers described below
- Growing season precipitation (May-Aug)(mm) => numerical variables in decimals or as integers; note: this value is only required if the user has a known value for this parameter. Otherwise leave this input as blank. AFFIRM.jl will estimate the value for this parameter from long-term precipitation probability distributions.
- Irrigation water amount, if irrigated (mm) => numerical variables in decimals or as integers; note: this value is only required if the user has a known value for this parameter. Otherwise leave this input as blank. AFFIRM.jl will use typical default value for this parameter.
- Nitrogen fertilizer product => categorical variables represented by identifiers described below
- Nitrogen fertilizer application timing => categorical variables represented by identifiers described below
- Nitrogen fertilizer application placement => categorical variables represented by identifiers described below
- Soil test nitrogen (0-24") (lb N/ac) => numerical variables in decimals
- Previous crop => categorical variables represented by identifiers described below
- Previous crop yield => numerical variables in decimals
- Previous crop yield unit => categorical variables represented by identifiers described below
- Residue management => categorical variables represented by identifiers described below
- Crop available nitrogen from applied manure (lb N/ac) => numerical variables in decimals
- Expected crop price ($/bu) => numerical variables in decimals
- Fertilizer price ($/tonne) => numerical variables in decimals
- Investment ratio => numerical variables in decimals
```
## User input identifiers for categorical variables
**==> Soil texture:**
```
"Very Coarse" => 1
"Coarse" => 2
"Medium" => 3
"Fine" => 4
"Very Fine" => 5
"Muck" => 6
"Peaty Muck" => 7
"Mucky Peat" => 8
"Peat" => 9
```
**==> Spring soil moisture:**
```
"Low" => 1
"Intermediate" => 2
"Optimum" => 3
```
**==> Crop:**
```
"Barley (Feed and Food)" => 1
"Barley (Hulless)" => 2
"Barley (Malt) " => 3
"Canola" => 4
"Canola (Argentine)" => 5
"Canola (Polish)" => 6
"Flax" => 7
"Oats" => 8
"Triticale (Spring)" => 9
"Wheat - Western Red Spring (WRS)" => 10
"Wheat - Northern Hard Red (NHR)" => 11
"Wheat - Western Amber Durum (WAD)" => 12
"Wheat - Western Extra Strong (WES)" => 13
"Wheat - Western Soft White Spring (WSWS)" => 14
```
**==> Irrigation:**
```
"No" => 1
"Yes" => 2
```
**==> Nitrogen fertilizer product:**
```
"ESN" => 1 
"ESN - Urea Blend (25:75)" => 2
"ESN - Urea Blend (50:50)" => 3
"ESN - Urea Blend (75:25)" => 4
"SuperU" => 5
"Urea" => 6
"Urea + eNtrench" => 7
"UAN (28-0-0) + Agrotain" => 8
"UAN (28-0-0)" => 9
"Anhydrous Ammonia" => 10
"Ammonium Nitrate" => 11
```
**==> Nitrogen fertilizer application timing:**
```
"Fall" => 1
"Spring" => 2
```
**==> Nitrogen fertilizer application placement:**
```
"Banded" => 1 
"Seed Placed" => 2
"Broadcast/incorporated (Surface banded)" => 3
"Broadcast" => 4
```
**==> Previous crop:**
```
"Alfalfa (Hay)" => 1
"Barley (Feed and Food)" => 2
"Barley (Hulless)" => 3
"Barley (Malt)" => 4
"Buckwheat" => 5
"Canary seed" => 6
"Canola" => 7
"Canola (Argentine)" => 8
"Canola (Hybrid)" => 9
"Canola (Juncea)" => 10
"Canola (Polish)" => 11
"Chickpeas" => 12
"Corn (Forage/Silage)" => 13
"Corn (Grain)" => 14
"Cowpeas" => 15
"Dry Bean (Black)" => 16
"Dry Bean (Great Northern)" => 17
"Dry Bean (Navy)" => 18
"Dry Bean (Pinto)" => 19
"Dry Bean (Shiny Black)" => 20
"Dry Bean (Small Red)" => 21
"Dry Bean (Yellow)" => 22
"Faba Bean" => 23
"Field Peas (Dun)" => 24
"Field Peas (Forage)" => 25
"Field Peas (Green)" => 26
"Field Peas (Maple)" => 27
"Field Peas (Processing)" => 28
"Field Peas (Red)" => 29
"Field Peas (Winter)" => 30
"Field Peas (Yellow)" => 31
"Flax" => 32
"Hay and Forage for Seed" => 33
"Lentils" => 34
"Lentils (Winter)" => 35
"Mustard (Brown)" => 36
"Mustard (Oriental)" => 37
"Mustard (Yellow)" => 38
"Oats (Feed)" => 39
"Oats (Forage/Silage)" => 40
"Oats (Hulless)" => 41
"Oats (Milling)" => 42
"Other Oilseed" => 43
"Other Pulse (Grain)" => 44
"Potatoes" => 45
"Rye (Fall)" => 46
"Rye (Spring)" => 47
"Safflower" => 48
"Soybeans" => 49
"Sugar beets" => 50
"Sunflowers" => 51
"Tame Hay (Legumes and Mix)" => 52
"Tame Hay (Other)" => 53
"Triticale (Spring)" => 54
"Triticale (Winter)" => 55
"Wheat - Northern Hard Red (NHR)" => 56
"Wheat - Prairie Spring Red (PSR)" => 57
"Wheat - Prairie Spring White (PSW)" => 58
"Wheat - Western Amber Durum (WAD)" => 59
"Wheat - Western Extra Strong (WES)" => 60
"Wheat - Western Hard White Spring (WHWS)" => 61
"Wheat - Western Red Spring (WRS)" => 62
"Wheat - Western Red Winter (WRW)" => 63
"Wheat - Western Soft White Spring (WSWS)" => 64
"Wheat - Western Special Purpose (WSP)" => 65
```
**==> Previous crop yield unit:**
```
"tons/ac" => 1 (must be used with previous crop ids - 1, 13, 45, 50, 52 and 53)
"bu/ac" => 2 (must be used with all previous crop ids except for the ones with above and below units)
"lb/ac" => 3 (must be used with previous crop ids - 6, 33 and 51)
```
**==> Residue management:**
```
"Soil Incorporated" => 1
"Removed from Field" => 2
"Removed by Burning" => 3
```
## Composite simulations for categorical variables
In composite simulations, a user is allowed to run sub-scenarios for a given variable in a scenario in a single row of the input file by inserting more than one input for that given variable within that row. Composite simulations are allowed for the following categorical variables:
```
- Soil texture
- Spring soil moisture
- Nitrogen fertilizer application timing
- Nitrogen fertilizer application placement
- Residue management
```
The composite simulations for the above categorical variables can be achieved by inserting more than one pipe separated identifier values of interests. For example, if a user wants to run both ```Fall``` and ```Spring``` scenarios for the variable ```Nitrogen fertilizer application timing```, the user can insert the input for that variable as ```1|2```
<br>
## Composite simulations for numerical variables
Composite simulations are also allowed for the following numerical variables:
```
- Soil organic matter (0-6")(%)
- Soil pH (0-6" or 0-12")
- Soil EC (0-6" or 0-12")(meq/100g)
- Growing season precipitation (May-Aug)(mm)
- Irrigation water amount, if irrigated (mm)
- Soil test nitrogen (0-24") (lb N/ac) => numerical variables in decimals
- Previous crop yield
- Crop available nitrogen from applied manure (lb N/ac)
- Expected crop price ($/bu)
- Fertilizer price ($/tonne)
- Investment ratio
```
The composite simulations for the above numerical variables can be performed in any of the following 3 ways:
<br>
<br>**1) Step-wise simulations:**
<br>
<br>Insert four (4) pipe-separated input values representing - lower limit, upper limit, composite simulation id (1 - for step-wise simulation), and step size or interval respectively; i.e., ```[lower limit]|[upper limit]|1|[interval]```
<br>For example, if a user is willing to run a step-wise simulation for a fertilizer price scenario between ```$650``` and ```$700``` per tonnes at a ```$10``` per tonne price interval, the input for the variable ```Fertilizer price ($/tonne)```will be ```650|700|1|10``` 
<br>
<br>**2) Monte-Carlo simulations with random uniform sampling:**
<br>
<br>Insert four (4) pipe separated input values representing - lower limit, upper limit, composite simulation id (2 - for Monte-Carlo simulations with random uniform samplings), and number of iterations respectively; i.e., ```[lower limit]|[upper limit]|2|[number of iterations]```
<br>If a user is willing to run the same fertilizer price scenario in the above example, but this time with a Monte-Carlo simulation with random uniform sampling repeated ```10```times, the input for the variable ```Fertilizer price ($/tonne)```will be ```650|700|2|10``` 
<br>
<br>**3) Monte-Carlo simulations with random normal sampling:**
<br>
<br>Insert four (4) pipe separated input values representing - average value, standard deviation value, composite simulation id (3 - for Monte-Carlo simulations with random normal samplings), and number of iterations respectively; i.e., ```[average]|[standard deviation]|3|[number of iterations]```
<br>If a user is willing to run a fertilizer price scenario with a Monte-Carlo simulation with random normal sampling repeated ```10```times, with an average fertilizer price of ```$675``` per tonne with a standard deviation of ```$20``` per tonne, the input for the variable ```Fertilizer price ($/tonne)```will be ```675|20|3|10```
<br>
## Execute AFFIRM.jl batch runs
```julia
julia> include("src/runAFFIRM.jl")
```
## Understanding AFFIRM.jl outputs
The output variables are written in the default ```output/AFFIRM-batch-outputs.csv``` file. The output variables in that file are listed below; visit [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) to understand these variables better.
```
- Township
- Range
- Meridian
- Soil Zone => Name of the agricultural soil zone of Alberta which the township falls under   
- Soil organic matter (0-6") (%)
- Soil texture
- Spring soil moisture
- Soil pH (0-6" or 0-12")
- Soil EC (0-6" or 0-12") (meq/100g)
- Crop
- Irrigation
- Growing season moisture flag => A description expressing whether it is a user input or long-term precipitation probability estimates representing Low, Intermediate or Optimum moisture conditions
- Growing season precipitation (May-Aug) + irrigation (if any) (mm)  => A description expressing whether it is a user input or typical precipitation + irrigation estimates representing Low, Intermediate or Optimum irrigation levels
- Nitrogen fertilizer product
- Nitrogen fertilizer application timing
- Nitrogen fertilizer application placement
- Soil test nitrogen (0-24") (lb N/ac)
- Previous crop
- Previous crop yield
- Previous crop yield unit
- Residue management
- Crop available nitrogen from applied manure (lb N/ac)
- Expected crop price ($/bu)
- Fertilizer price ($/tonne)
- User chosen investment ratio
- Estimated N release from N mineralization over the growing season (lb N/ac)
- N credit from previous crop residue (lb N/ac)
- Total plant available nitrogen from soil (lb N/ac) => Sum of Estimated N release from N mineralization over the growing season (lb N/ac), N credit from previous crop residue (lb N/ac), Soil test nitrogen (0-24") (lb N/ac), and Crop available nitrogen from applied manure (lb N/ac)
- Fertilizer N application rate (lb N/ac)
- Predicted crop yield (bu/ac)
- Predicted yield increase (bu/ac)
- Added yield increase (bu/ac)
- Estimated revenue from fertilizer N ($/ac)
- Marginal return or Gross margin change ($/ac)
- Total cost of fertilizer N ($/ac)
- Marginal cost of fertilizer N ($/ac)
- Estimated Investment Ratio
- Recommended? => A flag "Yes" indicates that the nitrogen rate in that row is an economically optimum nitrogen application rate and the predicted crop yield in that row is an economically optimum crop yield for the given scenario in that row.
```
