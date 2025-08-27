# Ignoring within-flower self-fertilization and inbreeding depression biases estimates of selection on floral traits

This README file was generated on 2024-07-07 by Kai-Hsiu Chen

General information: 

Title: Ignoring within-flower self-fertilization and inbreeding depression biases estimates of selection on floral traits
Journal of Ecology, 2024
Authors: Kai-Hsiu Chen and John R. Pannell
Affiliations: Department of Ecology and Evolution, University of Lausanne, Switzerland
Correspondance: Kai-Hsiu Chen ([kai-hsiu.chen@unil.ch](mailto:kai-hsiu.chen@unil.ch))
Date of data collection: 2022.5-8
Location: Solalex, Vaud, Switzerland

Data and R Codes for analyses in the paper
Within-flower selfing rate
R file: Within-flower_selfing.R
Data used: Data_selfing.csv
Output: Figure 2, 3

Phenotypic selection
R file: Phenotypic_selection.R
Data used: Data_phenotypic_selection.csv
Output: Figure 4

List and Details of Datasets

**Data_selfing.csv**
1.Number of variables: 11

2.Number of cases/rows: 60 individuals with single a bisexual flower

3.Variable List:
a)Ind_ID: individual identity
b)Flower_ID: flower identity of the individual
c)Stamen_number_func: Number of stamens after removal
d)Stamen_removal: to which of the three treatments subjected (Intact: intact; SR50: 50% removal; SR100: 100%)
e)Petal_length: tepal length (cm)
f)Stalk_height: height of stalk (cm)
g)Flowering_date : date of the begining of the female function (Julian date)
h)Pistil_number: number of pistils of the flower
i)Total_typed_seed: number of genotyped seeds
j)Total_selfed_seed: number of selfed seeds
k)Selfing_rate: within-flower selfing rate

**Data_phenotypic_selection.csv**: each flower has two rows each for the two scenarios of inbreeding depression (d= 0 or 0.95)
1.Number of variables: 15

2.Number of cases/rows: 120 (60 individuals with single a bisexual flower)

3.Variable List:
a)Ind_ID: individual identity
b)Flower_ID: flower identity of the individual
c)Stamen_number_func: Number of stamens after removal
d)Stamen_removal: to which of the three treatments subjected (Intact: intact; SR50: 50% removal; SR100: 100%)
e)Petal_length: tepal length (cm)
f)Stalk_height: height of stalk (cm)
g)Flowering_date: date of the begining of the female function (Julian date) *NA due to missing data
h)Pistil_number: number of pistils of the flower
i)Total_typed_seed: number of genotyped seeds
j)Total_selfed_seed: number of selfed seeds
k)Mature_seed_number: Female reproductive success using seed number as a proxy (selfed progeny contributes equally to both female and male RS)
l)Relative_RS: Female reproductive success using relative fitness as a proxy; k) relative to the mean of the population (selfed progeny contributes equally to both female and male RS)
m)Inbreeding_derpession: inbreeding depression scenarios (0 or 0.95)
n)Mature_seed_number_S: Female reproductive success using seed number as a proxy (selfed progeny contributes to only female functions)
o)Relative_RS_S: Female reproductive success using relative fitness as a proxy; n) relative to the mean of the population (selfed progeny contributes to only female functions)
