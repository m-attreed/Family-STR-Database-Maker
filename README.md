# family-str-db-maker

simulate a database of familial forensic STR profiles

## General Desciption
This program will take a database (CSV format) of DNA profiles (forensic STR) and output a database (CSV format) of simulated, familial STR profiles based on the input database. The program will have options such as number of children (e.g. 1, 3) and number of generations (e.g. -1, 1, 2, -2).

## Running the Program
`family-str-db-maker.py <input DB filename> <output DB filename> <children generation> <parent generation>`


## Model

### Allele Frequencies
The program uses the updated 2017 NIST Allele Frequencies. I added data to the table to fill in blank values between the lowest and highest observed alleles. This is to allow the program to use a low rate (0.01 times the lowest frequency recorded in the total set) instead of nothing for rare alleles.

### Children

### Generations
Here we will be using numerical values to build a range of generations. A negative sign indicates children, grandchildren, etc. produced with a simulated spouse. A positive integer indicates parents, grandparents, etc. Parents generation includes aunts and uncles based on the number of children specified plus the generation (parents 1, grandparents 2) as that seems fair based on changing family sizes.

Consider adding family archetype over each profile out of simplicity and repeatability. Maybe consider a family including children and grandchildren, and including parents and grandparents. Each generation allows for the making of siblings. To make aunts and uncles you need grandparents to get the alleles. All spouses will be randomly generated.

### Profile Generation
All profiles randomly generated will use NIST Allele Frequency Table as the weights for the alleles. Profiles using "rare" unseen alleler between the smallest recorded allele and the largest are possible at a rate of 0.01 times the lowest allele frequency recorded.

## Input
Input file formating:
Single column for each locus using comma separated alleles in phenotype or genotype format.
Multiple columns for each locus with single allele per column

## Output

## Assumptions
Indepedent assortment of alleles.
