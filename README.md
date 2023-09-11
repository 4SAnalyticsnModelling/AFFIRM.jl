![Model](https://github.com/4SAnalyticsnModelling/AFFIRM.jl/blob/main/resource/affirm_logo.png) # AFFIRM.jl 
Alberta Farm Fertilizer Information Recommendation Manager (AFFIRM). AFFIRM has 3 concurrent versions - [AFFIRM v3.0](https://www.alberta.ca/alberta-farm-fertilizer-information-and-recommendation-manager), [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) and AFFIRM.jl. The AFFIRM v3.0 is the production version of AFFIRM. The AFFIRM-R is the nitrogen sub-model of AFFIRM v3.0. The AFFIRM.jl provides batch run utilites for multiple scenarios and provides similar functionalities as in AFFIRM-R.

## Get started
Install AFFIRM.jl
```julia
using("Pkg");
Pkg.add("https://github.com/4SAnalyticsnModelling/AFFIRM.jl");
```
Create AFFIRM folders
```julia
using AFFIRM;
cd("your/project/folder");
create_affirm();
```
## AFFIRM folder structure
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
## Execute AFFIRM batch runs
Execute ```src/runAFFIRM.jl``` file
## Understanding AFFIRM input variables
List of user input variables in ```AFFIRM-batch-inputs.csv``` file; visit [AFFIRM-R](https://mezbahu.shinyapps.io/AFFIRM_R_version_yield_response_nitrogen/) to understand the variables more.
```
- Township
- Range
- Meridian
- Soil organic matter (0-6")(%)
- Soil texture
- Spring soil moisture
- Soil pH (0-6" or 0-12")
- Soil EC (0-6" or 0-12")(meq/100g)
- Crop
- Irrigation
- Growing season precipitation (May-Aug)(mm)
- Irrigation water amount, if irrigated (mm)
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
- Investment ratio
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
