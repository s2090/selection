# A test of the Grant-Stebbins pollinator shift model of floral evolution

[https://doi.org/10.5061/dryad.2280gb637](https://doi.org/10.5061/dryad.2280gb637)

## Description of the data and file structure

This data was collected using floral measurements at nine localities in the Eastern Cape of South Africa. At seven of these localities, we documented visitation rates and contact with floral reproductive parts to determine pollinator importance. We additionally caught and measured pollinators for their functional body traits. At three of these localities, we also conducted phenotypic selection experiments, in which we measured floral traits from several hundred individuals and female fitness.

### Files and variables

#### File: 18.\_Bees\_and\_Bee\_Flies\_Contact.csv

**Description:** Data used to analyse the contact of widespread pollinators with floral reproductive parts at localities where they were most important and also less so

##### Variables

* Ecotype: floral variety
* Population: Locality 
* Pollinator: functional pollinator
* Contact: presence-absence data (1: contact; 0: no contact)

#### File: 17.\_AM\_Selection\_Data.csv

**Description:** 

##### Variables

* population: Amakhala
* ind.: self explanatory
* TW: tepal width (mm)
* TL: tepal length (mm)
* SL: style length (mm)
* PL: petal length (mm)
* HT: inflorescence height (mm)
* FL: flower number (mm)
* seeds: seed set per individual

#### File: 16.\_RH\_Selection\_Data.csv

**Description:** 

##### Variables

* population: Rhini
* ind.: self explanatory
* TW: tepal width (mm)
* TL: tepal length (mm)
* SL: style length (mm)
* PL: petal length (mm)
* HT: inflorescence height (mm)
* FL: flower number (mm)
* seeds: seed set per individual

#### File: 15.\_GE\_Selection\_Data.csv

**Description:** 

##### Variables

* population: Good Earth
* ind.: self explanatory
* TW: tepal width (mm)
* TL: tepal length (mm)
* SL: style length (mm)
* PL: petal length (mm)
* HT: inflorescence height (mm)
* FL: flower number
* seeds: seed set per individual

#### File: 13.\_Short-Styled\_BS\_Seed\_Set.csv

**Description:** Data obtained from breeding system experiments associated with the short style site

##### Variables

* Species: B.gregaria
* Individual: Individual worked with
* Ecotype: floral variety
* Treatment: experimental treatment -
  * control treatment: Control
  * outcross treatment: OutCross
  * geitonogamy treatment:geit
  * facilitated selfing:FacSelf 
* Locality: site where experimental plants originated from
* Fruit_set: binary (1: present; 0: absent)
* fertilized_ovules: fertilized ovules per fruit capsule
* Aborted_ovules: aborted ovules per fruit capsule
* Total_ovules: all ovules, aborted and fertilized

#### File: 14.\_Intermediate\_BS\_Seed\_Set.csv

**Description:** Data obtained from breeding system experiments associated with the intermediate style site

##### Variables

* Species: B.gregaria
* Individual: Individual worked with
* Ecotype: floral variety
* Treatment: experimental treatment -
  * control treatment: Control,
  * outcross treatment: OutCross,
  * geitonogamy treatment:geit,
  * facilitated selfing treatment: FacSelf 
* Locality: site where experimental plants originated from
* Fruit_set: binary (1: present; 0: absent)
* fertilized_ovules: fertilized ovules per fruit capsule
* Aborted_ovules: aborted ovules per fruit capsule
* Total_ovules: all ovules, aborted and fertilized



#### File: 12.\_Long-Styled\_BS\_Seed\_Set.csv

**Description:** Data obtained from breeding system experiments associated with the long style site

##### Variables

* Species: B.gregaria
* Individual: Individual worked with
* Ecotype: floral variety
* Treatment: experimental treatment -
  * control treatment:Control
  * outcross treatment: OutCross
  * geitonogamy treatment: geit
  * facilitated selfing treatment: FacSelf 
* Locality: site where experimental plants originated from
* Fruit_set: binary (1: present; 0: absent)
* fertilized_ovules: fertilized ovules per fruit capsule
* Aborted_ovules: aborted ovules per fruit capsule
* Total_ovules: all ovules, aborted and fertilized

#### File: 5.\_Individual\_PI\_and\_PFL.csv

**Description:** Data set from which we have calculated pollinator importance

##### Variables

* locality: any of the three short, intermediate and long style sites
* Pollinator: pollinator functional group
* PI: Pollinator importance
* PI adj: Pollinator importance adjusted
* prop.PI: Proportion pollinator importance
* length: functional length of pollinators (mm)

#### File: 3.\_Visitation\_Rate.csv

**Description:** Used in visitation rate analysis as the number of visits per flower per hour (i.e. network)

##### Variables

* Pollinator: pollinator taxa
* MK: Makhanda
* AM: Amakhala
* GE: Good Earth
* KA: Kariega
* NN: Nanaga
* PA: Paterson
* RH: Rhini

#### File: 4.\_Contact\_Frequency.csv

**Description:**  Dataset used in contact analysis (i.e. contact frequency network - not shown)

##### Variables

* Pollinator: pollinator taxa
* MK: Makhanda
* AM: Amakhala
* GE: Good Earth
* KA: Kariega
* NN: Nanaga
* PA: Paterson
* RH: Rhini

#### File: 2.\_Pollinator\_Functional\_Length.csv

**Description:** Dataset associated with the functional lengths of pollinators in mm

##### Variables

* Pollinator:
* MK: Makhanda
* AM: Amakhala
* GE: Good Earth
* KA: Kariega
* NN: Nanaga
* PA: Paterson
* RH: Rhini
* PW: Petworth
* RE: Riebeeck East

#### File: 1.\_Floral\_Traits.csv

**Description:** dataset associated with floral measurements in mm

##### Variables

* locality: site
* Ecotype: floral variety
* SL: Style length (mm)
* TL: Tepal length (mm)
* PL: Pedicel length (mm)

#### File: 19.\_Floral\_Trait\_Means.csv

**Description:** Used in data analysis in which the means of floral traits are required, measurements are in mm

##### Variables

* locality: site
* Ecotype: floral variety
* SL: Style length (mm)
* TL: Tepal length (mm)
* PL: Pedicel length (mm)

#### File: 8.\_Tepals\_and\_Guides\_Intermediate.csv

**Description:** Dataset associated with colour spectra from an intermediate site, in which measurements were taken from the centre of the tepal: mid and end of tepal associated with the dominant colour of the tepal: side. 

##### Variables

* wavelength: wavelength in nanometer
* Side: see in description
* Mid: see in description

#### File: 9.\_Tepals\_and\_Guides\_Short-Styled.csv

**Description:** Dataset associated with colour spectra from a short style site, in which measurements were taken from the centre of the tepal: mid and end of tepal associated with the dominant colour of the tepal: side. 

##### Variables

* wavelength: wavelength in nanometer
* Side: see in description
* Mid: see in description

#### File: 6.\_Tepal\_Spectra\_All\_Ecotypes.csv

**Description:** tepal spectra compiled for the short-styled site, intermediate and long style sites

##### Variables

* wavelength: wavelength in nanometers
* Short-styled: refers to individuals from the short style site from which the spectra was measured
* Intermediate: refers to individuals from the intermediate style site from which the spectra was measured
* Long-styled: refers to individuals from the long style site from which the spectra was measured

#### File: 7.\_Tepals\_and\_Guides\_Long-Styled.csv

**Description:** Dataset associated with colour spectra from an long style site, in which measurements were taken from the centre of the tepal: mid and end of tepal associated with the dominant colour of the tepal: side. 

##### Variables

* wavelength: wavelength in nanometer
* Side: see in description
* Mid: see in description

#### File: 11.\_Eristalis\_Spectral\_Sensitivities.csv

**Description:** Hoverfly sensitivities that were imported into the fly colour vision model

##### Variables

* Wavelength: wavelength in nanometers, self explanatory
* R7p: photoreceptor sensitivity 
* R7y: photoreceptor sensitivity 
* R8p: photoreceptor sensitivity
* R8y: photoreceptor sensitivity

#### File: 10.\_Background\_Spectra.csv

**Description:** Background spectra averaged from several items occurring in the general area in which the plants grow

##### Variables

* 1 nm: wavelength in  1nm intervals
* Average: average spectra taken from several items in the background, from the area in which the plants grow.

