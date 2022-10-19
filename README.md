# Introduction
**r-bats** (Bayesian Analysis of Tip Significance in R) is an R library for the analysis of phylogenetic data. It asks the question: 'given this distribution of traits on taxa, how likely is this pattern to have arisen by chance?'. It assumes discrete trait values, known for each taxon unambiguously. The phylogenetic relationship between the taxa must be represented as a posterior set of trees, e.g. a collection of trees assumed to correctly sample the posterior distribution of the phylogeny. Usually these will have been inferred from molecular data using [BEAST](https://www.beast2.org/) or [MrBayes](http://mrbayes.sourceforge.net/). 


# Requirements
R version 4.0.0 or above is required to download this package. The package "devtools" is also needed to be able to download packages off of github. 


# Instalation
To install this package, the devtools package needs to be installed and then used like this:
```
devtools::install_github("slequime/r-bats")
```

# Usage
The only function a user needs to call to use this package is `bats()`. This function can take four arguments:

| Argument  | Optional | Description                                                                            |
|-----------|----------|----------------------------------------------------------------------------------------|
| treefile  | No       | The path to the .trees or .tre file generated by BEAST                                 |
| xmlfile   | No       | The path to the .xml file                                                              |
| reps      | Yes      | An integer indicating how many shuffled trees should be generated for each normal tree |
| userinput | Yes      | Which attribute should be used as the state attribute                                  |

Example:
```
treefile <- "path/to/files/treefile.trees"
xmlfile <- "path/to/files/xmlfile.xml"

bats(treefile, xmlfile, reps=1)
```


# Contact
If something is wrong with the package or calculations, or something should be added to this package you can contact me at [email address]







