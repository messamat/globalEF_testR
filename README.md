## R code for 'Global and local estimates of environmental flow requirements to sustain river ecosystems are poorly correlated'

This repository contains R code associated with _Messager, M. L., Dickens, W. S. C., Eriyagama, N., Tharme, R. E., Stassen, R. (In review). Global and local estimates of environmental flow requirements to sustain river ecosystems are poorly correlated. [to be updated]_

## Abstract
Environmental flows (e-flows) are a central element of sustainable water resource management to mitigate the impacts of 
hydrological alteration on freshwater ecosystems and their benefits to people. Many nations protect e-flows, and thousands
of e-flow assessments have been conducted globally, leveraging local data and knowledge to quantify how much water must be
kept instream for sustaining healthy ecosystems. However, e-flow assessments and implementation are geographically uneven 
and cover a small fraction of rivers worldwide, which hinders globally consistent target-setting, monitoring and evaluation 
for international agreements to curb water scarcity and biodiversity loss like the UN Sustainable Development Goals. Therefore, 
global models have been developed to estimate the e-flow requirements of global rivers seamlessly across basins and administrative 
boundaries. But there has been little effort to benchmark these models against locally derived e-flow estimates, which may 
limit confidence in the relevance of global targets. The aim of this study was thus to assess whether current global methods 
reflect e-flow estimates used on the ground by comparing global and local estimates for 1194 sites across 25 countries. 
While global approaches can broadly approximate the bulk amount of water that should be precautionarily set aside to 
sustain aquatic ecosystems at the scale of large basins or countries, they explain a negligible 0-1% of the variability 
in locally derived estimates of the percentage of river flow that must be protected at a given site. Such a disconnect
between global and local assessments of e-flow requirements limits the credibility of current global e-flow estimates and 
associated targets for human and ecosystem water use. To accelerate the global implementation of e-flows requires a concerted
effort to compile and draw from the thousands of local e-flow assessments that have been conducted worldwide to bridge the 
gap from local to global scales.

## Introduction

This repository includes the portions of the analysis conducted in R, which encompass all data analysis for comparing global and local estimates
of Mean Annual Flow (MAF) and e-flows. Files needed to run this analysis are available by downloading the study's figshare permanent repository. 
The /data folder in the figshare repository contains raw data and the directory structure enables users to reproduce our 
study using the scripts herein. The pre-formatting of e-flow sites and the calculation of global MAF and e-flow estimates were performed
in Python with script available at: https://github.com/messamat/globalEF_testPy. 

These scripts are annotated but could be challenging to follow. If you encounter any trouble, please don't hesitate
to contact Mathis L. Messager for comments and clarifications by email or to log an issue in github.

## Analysis structure and underlying data

This analysis relies as much as possible on [good enough practices in scientific computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510), which users are encouraged to read.

**Structure**: the overall project directory is structured with the following sub-directories:  
data/ (raw data, read-only, not to be altered)  
results/ (results of the analysis, mostly reproduceable through code execution. However, also includes manually modified results)
src/ (code written for the project)  
|---- globalEF_testR (R project/source code for analysis in R)  

All scripts rely on this structure.

**R Workflow**: this project is setup with a [targets workflow](https://docs.ropensci.org/targets/), ensuring reproducibility.
In the `targets` philosophy, every action is a function, and every R object resulting from a workflow step is a "target" with dependencies.
Intermediate targets/objects are stored in a `_targets` directory. 

**Dependency management**: the R library of this project is managed by [renv](https://rstudio.github.io/renv/articles/renv.html).
This makes sure that the exact same package versions are used when recreating the project.
When calling `renv::restore()`, all required packages will be installed with their specific version. 
Please note that this project was built with R version 4.1.2. on a Windows 10 operating system.
The renv packages from this project **are not compatible with R versions prior to version 3.6.0.**

**Syntax**: this analysis relies on the [data.table](https://rdatatable.gitlab.io/data.table/) syntax, which provides a high-performance version of data.frame. It is concise, faster, and more memory efficient than conventional data.frames and the tidyverse syntax.

## Getting started
### Download the repository for R
In Git Bash, the following commands illustrate the procedure to make a local copy of the Github repository in a newly created directory at 
C://globalEF_testR/src :

```{r, engine = 'bash', eval = FALSE}
Mathis@DESKTOP MINGW64 /c/globalEF_testR/src
$ git clone https://github.com/messamat/globalEF_testR.git
```

In R Studio for Windows, the following procedure can be used:  

* Click on “File” in the menu ribbon  
* Select “New project…”  
* Choose the “Version control” option in the New Project Wizard window.
* Then, select “Git” in the next window.
* In the next window, fill the fields as follows:  
  * Repository URL: https://github.com/messamat/globalEF_testR
  * Project directory name: [will autofill as “globalEF_testR”]  
  * Create project as subdirectory of: [choose the parent directory of src, e.g., C://globalEF_testR//src]  
* Tick “Open in new session” and then click “Create project”.  


### Github repository structure
- [**R/**](https://github.com/messamat/globalEF_testR/tree/main/R) — core of the analysis
  - [*functions.R*](https://github.com/messamat/globalEF_testR/blob/main/R/functions.R) - all custom functions used in the data formatting and analysis. 
  - [*packages.R*](https://github.com/messamat/globalEF_testR/blob/main/R/packages.R) - all packages used in the workflow.
- [*.Rprofile*](https://github.com/messamat/globalEF_testR/blob/main/.Rprofile) — used to activate renv for new R sessions launched in the project.
- [*globalEF_testR.Rproj*](https://github.com/messamat/globalEF_testR/blob/main/GeneticScaling.Rproj) — R project file.
- [*LICENSE*](https://github.com/messamat/globalEF_testR/blob/main/LICENSE) - terms of use, modification and sharing for this software.
- [*README.md*](https://github.com/messamat/globalEF_testR/blob/main/README.md) — README for Github (this file)
- [*\_targets.R*](https://github.com/messamat/globalEF_testR/blob/main/_targets.R) — configuration script for targets workflow,  this specific file name is required by the targets package. Contains the targets “plan”, the high-level catalog of all the steps in the workflow (see the corresponding chapter in the targets user manual). This plan defines the order of functions to use, their inputs and outputs (usually, targets), and the relationship among targets and steps.
- [*renv.lock*](https://github.com/messamat/globalEF_testR/blob/main/renv.lock) — renv lockfile, describing the state of the project’s library (installed packages and their version).


## Running the analysis
Provided that your were given the necessary data, the entire analysis can simply be re-run with the following code found in 
```{r rmake, eval = FALSE}
source('_targets.R')
tar_make()
```
`tar_make()` is the central function of the targets approach. It runs all the steps of the workflow in the correct order, skipping any work that is already up to date. Because of how targets tracks global functions and objects as dependencies of targets, the use of `tar_make()`  is needed to run the analysis pipeline in a clean reproducible environment. If all targets are up to date in the caching directory, then nothing will be run.

## Inspecting results
If you were provided intermediate targets (i.e., a `_targets/` directory; or once you have re-run the analysis), you can load individual targets in the environment with the following commands (even if the targets are not up to date due to e.g. a change in source path). 
``` {r loadtarg, eval = FALSE}
tar_load(countrytab) #Load target in memory (R environment) with original target name as variable name 
countrytab <- tar_read(countrytab) #Load target in memory with new variable name
```

### Notes and resources 
* The [issues tracker](https://github.com/messamat/globalEF_testR/issues) is the place to report problems or ask questions 
* See the repository [history](https://github.com/messamat/globalEF_testR/issues/commits/master) for a fine-grained view of progress and changes.