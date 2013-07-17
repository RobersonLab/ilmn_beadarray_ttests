ilmn_beadarray_ttests
=====================

## Introduction
It is common to perform transcriptome profiling by array and then analyze the differential expression between different conditions. A stanard way of calculating expression for R users is to apply the limma package to calculate which genes are differentially expressed. But what if you have only one test and one control sample to compare?

If the data are Illumina Beadarray derived there is one potential answer. Each probe for each array provides the Detection p-value, bead count, average intensity and standard error. Using these data there is sufficent information to perform a standard t-test between two samples. One challenge is how to normalize the data. The best solution I have discovered is to perform median quantile normalization making sure to treat all four of those measures as one data rather than just substituting averages.

The master function for t-test calculation is customIlmnTtest. Others are support functions. Doing a 1v1 expression test is not recommended when at all possible. The number of significant p-values is vastly overinflated for this method. Therefore I suggest focusing on absolute fold-changes > 2.00. Greater than 2.25 and 2.50 are probably even less likely to be false-positives. 

## Requirements

#### R
beadarray (Bioconductor)

## Input Data
The only summary data used as input should be non-normalized, as normalization is typically applied only to the average values. Data can be either background subtracted or non-background subtracted. The script will adjust the average intensity offset as necessary.

## Usage

#### customIlmnTest Options
raw_object: the beadarray imported (readBeadSummaryData) array data.

test: The NAME (colname) of the test experiment. NOTE: case-sensitive and will match only on the column name of the sample, not the full name of the column.

control: The NAME (colname) of the control experiment.

log2transform: Should the data (including errors) be appropriately log2 transformed.

qNormalize: Should the data be quantile normalized?

qMethod: The only current method is "median". Added for future flexibility.

annotate.frame: A data frame containing AT LEAST ProbeID and ILMN_GENE for the array. Only these two columns are used, though others may be included.

showBox: Should a boxplot be plotted to show pre- to post-normalization bar charts of average intensity?

#### Example code
    # Array data is in file "myData.txt"
    # Annotation of ProbeID and ILMN_GENE is in "myAnnotation.csv"
    # control column is "myControlExp"
    # test column is "myTestExp"
    
    myAnnotation = read.csv("myAnnotation.csv", header=TRUE)
    
    require(beadarray)
    source("ilmn_beadarray_ttests.R")
    
    raw_bead_summary_data = readBeadSummaryData("myData.txt", skip=0)
    
    ttest.results = customIlmnTtest( raw_bead_summary_data, test="myTestExp", control="myControlExp", removeUnexp=TRUE, log2transform=TRUE, qNormalize=TRUE, qMethod="median", annotate.frame=myAnnotation, showBox=FALSE)
    
    head(ttest.results)
      ProbeID ILMN_GENE FC myTestExp/myControlExp q.value Bonferroni AbsFC myTestExp myControlExp
    3 2570615 A1BG 18.206338 3.020714e-16 3.020714e-16 18.206338 1049.60567 57.65056
    9 6450255 7A5 -6.126913 5.076439e-16 1.015288e-15 6.126913 171.31068 1049.60567
    8 6370619 A1BG -2.040916 3.685755e-06 1.105726e-05 2.040916 83.93814 171.31068
    1 1580181 A2BP1 -1.293890 2.003013e-01 8.012052e-01 1.293890 64.87270 83.93814
    5 2650615 A1CF -1.125274 9.029617e-01 1.000000e+00 1.125274 57.65056 64.87270
    2 2000519 A26C3 1.000000 1.000000e+00 1.000000e+00 1.000000 86.13228 86.13228
      DegreesOfFreedom FC_Err FC_Lower FC_Upper T p.value
    3 22.19567 4.6274310 13.5789072 22.833769 23.6794147 3.020714e-17
    9 17.65227 -0.7682736 -5.3586396 -6.895187 -30.5023526 1.015288e-16
    8 37.26916 -0.5189858 -1.5219300 -2.559902 -5.8106518 1.105726e-06
    1 48.54424 -0.3844967 -0.9093933 -1.678387 -1.7873765 8.012052e-02
    5 42.58171 -0.3617858 -0.7634885 -1.487060 -0.7599345 4.514808e-01
    2 54.00000 0.3300649 0.6699351 1.330065 0.0000000 1.000000e+00```
    