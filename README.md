Data repository for stabilising selection meta-analysis.

# Contents

`data.xlsx` - Datasheet with information necessary to perform analyses

`code.R` - R code for selection analyses using study data from `data/`

`meta-analysis.R` - R code for meta-analysis using result files in `results/` generated with `code.R`

#### `data/`

Public domain datasets from the studies analysed. Some columns may have been added to enable analyses, see notes column in `data.xlsx`.

#### `results/`

Folder containing generated RDS and XLSX files from `code.R` and `meta-analysis.R`

# Data Management Plan
### Types of data generated
RDS files containing results from selection analyses and meta-analysis.

### Types of data preserved
All of the above and the components used to generate it.

### Software and metadata implications
The datasheet is readable in both R and Excel and analysis code will be readable in R.

### Length of data preservation
Indefinitely.

### Value of data to others
Code and generated data can be used and modified for future analyses.

### How data will be shared
All data and code is available on GitHub at [https://github.com/s2090/selection](https://github.com/s2090/selection) with a README file explaining each file's purpose.
