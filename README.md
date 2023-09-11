![AFFIRM-logo](https://github.com/4SAnalyticsnModelling/AFFIRM.jl/blob/main/resource/affirm_logo.png) 
# AFFIRM.jl 
Alberta Farm Fertilizer Information Recommendation Manager (AFFIRM) has 3 concurrent versions - [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager), [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) and AFFIRM.jl. The [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager) is the production version of AFFIRM. The [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) is the nitrogen sub-model of [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager). The AFFIRM.jl provides batch run utilites for multiple scenarios and provides similar functionalities as in [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/).

## Get started with AFFIRM.jl
Install AFFIRM.jl
```julia
using Pkg;
Pkg.add("https://github.com/4SAnalyticsnModelling/AFFIRM.jl");
```
Create AFFIRM.jl folders
```julia
using AFFIRM;
cd("your/project/folder");
create_affirm();
```
## AFFIRM.jl folder structure
AFFIRM.jl folder structure
```
> data
> input
  - AFFIRM-batch-inputs.csv
> output
> src
  - runAFFIRM.jl
```
***Do not change any folder or file names***
## Execute AFFIRM.jl batch runs
Execute ```src/runAFFIRM.jl``` file
## Understanding AFFIRM.jl input variables
List of user input variables in ```AFFIRM-batch-inputs.csv``` file with their data types; visit [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) to understand the variables more.
```
- Township => Alberta township id; integer values ranging from 1 to 126
- Range => Alberta range id; integer values ranging from 1 to 30
- Meridian => Alberta meridian id; categorical variables: "W4", "W5" and "W6"
- Soil organic matter (0-6")(%) => numerical values in decimals
- Soil texture => categorical values represents by identifiers described below
- Spring soil moisture => categorical values represents by identifiers described below
- Soil pH (0-6" or 0-12") => numerical values in decimals
- Soil EC (0-6" or 0-12")(meq/100g) => numerical values in decimals
- Crop => categorical values represents by identifiers described below
- Irrigation => categorical values represents by identifiers described below
- Growing season precipitation (May-Aug)(mm) => numerical values in decimals or as integers
- Irrigation water amount, if irrigated (mm) => numerical values in decimals or as integers
- Nitrogen fertilizer product => categorical values represents by identifiers described below
- Nitrogen fertilizer application timing => categorical values represents by identifiers described below
- Nitrogen fertilizer application placement => categorical values represents by identifiers described below
- Soil test nitrogen (0-24") (lb N/ac) => numerical values in decimals
- Previous crop => categorical values represents by identifiers described below
- Previous crop yield => numerical values in decimals
- Previous crop yield unit => categorical values represents by identifiers described below
- Residue management => categorical values represents by identifiers described below
- Crop available nitrogen from applied manure (lb N/ac) => numerical values in decimals
- Expected crop price ($/bu) => numerical values in decimals
- Fertilizer price ($/tonne) => numerical values in decimals
- Investment ratio => numerical values in decimals
```
## User input identifiers for categorical variables
Soil texture:
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
