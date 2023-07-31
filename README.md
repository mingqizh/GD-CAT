# GD-CAT
Genetically-derived correlations across tissues, or `GD-CAT`is a shiny-based endocrine genetic web app for tissue cross-talking. Now we have both human data in sex and mouse data in diets. <br><br>
![initial]([https://github.com/mingqizh/GD-CAT/blob/main/images/fgt.png](https://github.com/mingqizh/GD-CAT/blob/main/images/pipeline.png))<br><br>
You can select the gene and tissue combination that interests you, and have a look at the demographic characteristics from the cohort of choice.
and then do the enrichment, scatter plot, single cell and network visualization.  <br><br>
You can access our app via [GD-CAT](https://pipeline.biochem.uci.edu/gtex/demo2/).<br><br>
Now the [paper](https://pubmed.ncbi.nlm.nih.gov/37214953/) is preprint.<br><br>
## Tutorial 
### Initial settings
First start from the dashboard. Select the sex that you want and input an NBCI gene symbol, or you will see a pop-up warning. <br><br>
![initial](https://github.com/mingqizh/GD-CAT/blob/main/images/fgt.png)<br><br>
Hit the "process data" button after choosing a tissue. A progress bar will show at the bottom right of the screen. <br><br>
![pre](https://github.com/mingqizh/GD-CAT/blob/main/images/1689177951994.png)<br><br>
After the processing of the original data, you can see the pie chart where the gene is enriched in, age and sex of the cohort, cell type, and top genes. <br><br>
![pie](https://github.com/mingqizh/GD-CAT/blob/main/images/sfd.png)<br><br>
### Enrichment Analysis
To generate the enrichment analysis of the second part, users should click the pie chart body to select a tissue and threshold first. Your choice will be displayed then in this section. Then users can hit the "Start Analysis" button to see the results. <br><br>
![enrichment](https://github.com/mingqizh/GD-CAT/blob/main/images/enrich.png)<br><br>
### Scatter Plot
Same as the initial settings, make sure you put an official NBCI gene symbol here. <br><br>
![scatter](https://github.com/mingqizh/GD-CAT/blob/main/images/1689177689112.png)<br><br>
### Network Analysis
Move the slider bars to custom the gene numbers for within and peripheral tissue and then process the analysis. <br><br>
![net](https://github.com/mingqizh/GD-CAT/blob/main/images/net.png)
#### Attention
- All the images plus some of the data that generate them can be downloaded. The interactive images download buttons are on the top right and for other plots are from the bottom left. <br><br>
- You should hit the button again to update your settings. 

