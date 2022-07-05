# Figures
Scripts to draw the pretty or informative figures

## 1: SWNE Plot
An improved version of [SWNE](https://yanwu2014.github.io/swne/). More details can be seen in ./SWNE_Plot/Figures/TrajectoryInference0829.html
![Improved SWNE](https://github.com/chengwenxuan1997/Figures/blob/main/SWNE_Plot/Figures/SWNE.png)

### 1.1 Purposes
1: Run Geneset Enrichment Analysis (GSEA) on the gene loadings for each factor to find the top genesets associated with that factor after the SWNE embedding is calculated. This function is present in the SWNE article but not included in the SWNE R package.  
    **Correction: The function is also included in the SWNE packages, just not present in the tutorial.**  
2: Take a try to choose the root node via an entropy-based method in the trajectory analysis because both monocle2 and monocle3 can not automatically provide a root node. The idea was inspired by [StemID](https://doi.org/10.1016/j.stem.2016.05.010).  


## 2: Wheel plot
A wheel plot program rewritten with R.
![](https://github.com/chengwenxuan1997/Figures/blob/main/WheelPlot/Figures/wheelplot.png)

### 2.1 Purposes
Translate [linnarsson's python script](https://github.com/linnarsson-lab/ipynb-lamanno2016) into Rscript for the convenience of R users.

### 2.2 Notes
1: Each node represents a sample, and each vertex of the polygon represents a class label. The closer the sample point is to the vertex, the more likely the sample is to be classified into this category.  
2: The **WheelPlot.R** tried to reproduce the whole pipeline shown in the python script but failed because the classification algorithm in R is not the same as in python. The **WheelPlot_Simple_V2.R** only includes codes to draw a wheel plot from a provided predictive probability table.  

## 3: DEMaP
A cpp implement of **D**enoised **E**mbedding **MA**nifold **P**reservation ([DEMaP](https://github.com/scottgigante/DEMaP))

### 3.1 Purpose
The fast calculation of DEMaP is the premise of wide application, and can also bring convenience to the evaluation of data visualization methods.  

### 3.2 Notes
1: The parallel cpp version is available, and it only costs two minutes on a 30,000 cell dataset(threads=10).  
2: The cpp program can be sourced using **EvaluationMetrics.R**.
