
### Updated R code and data from Gasparrini Epidem 2011

--------------------------------------------------------------------------------

A comparative analysis of the main and added effects of temperature on mortality. The code originally reproduced the example included in the article:

Gasparrini A, Armstrong B. The impact of heat waves on mortality. *Epidemiology*. 2011;**22**(1):68-73. [[freely available here](http://www.ag-myresearch.com/2011_gasparrini_epidem.html)]

The original example included in the article was based on data for the 108 USA cities available from the National Mortality, Morbidity, and Air Pollution Study (NMMAPS), which at the time of the publication were available through the R package NMMAPSlite. Unfortunately, the data are not available any more and the package NMMAPSlite has been archived. This means that the analysis of the paper is not replicable, unless you have access to the original data and you modify the code accordingly.

Functions in the in the R package [dlnm](https://github.com/gasparrini/dlnm) are used for modelling the main effect and graphically representing the added effect, while functions in the R package [mvmeta](https://github.com/gasparrini/mvmeta) are applied for pooling the results from multiple studies.


--------------------------------------------------------------------------------

The material:

  * *01.main.R* reproduces the main results in the analysis
  * *02.additional1.R* and *03.additional2.R* reproduces additional results
  * *04.sensitivity.R* performs the sensitivity analyses

Download as a ZIP file using the green button *Clone or download* above