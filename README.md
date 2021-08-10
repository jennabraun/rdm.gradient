# Arthropod associations with residual dry matter and foundations shrubs in California desert ecosystems

Pitfall traps were deployed along an aridity gradient in California to measure associations of ground-dwelling arthropods with foundation shrubs, and to assess the ability of residual dry matter (RDM) as an indicator of arthropod community structure.

![Panoche Hills](/panoche.jpg)

## Getting Started

These instructions describe the contents of this repository and will get you a copy of the project up and running on your local machine. 

See https://jennabraun.github.io/rdm.gradient for final work up 

### Data descriptions

* Raw Data folder holds excel workbooks used for data entry
* Clean Data folder holds .csv files for analysis along with metadata describing the variables
* Figures are the final figures
* gpx holds the .gpx for sampling locations
* scripts hold R scripts for processing the datasets


### Script descriptions
* cleaning.R takes the raw data workbooks as inputs and writes to clean data folder
* RII.R calculates RII and outputs results as .csv which are read into index.rmd
* ordination.R is the CCA and indicator species code
* worldclim.R downloads worldclim data for coordinates


