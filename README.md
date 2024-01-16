# GD-CAT
Genetically-derived correlations across tissues, or [GD-CAT](https://pipeline.biochem.uci.edu/gtex/) is a shiny-based endocrine genetic web app for tissue cross-talking. Now we have both human data in sex from [GTEx](https://gtexportal.org/home/) and mouse data in diets from [HMDP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878195/). <br><br>
![initial](https://github.com/mingqizh/GD-CAT/blob/main/images/pipeline.png) <br><br>
First, select the gene and tissue combination of interest, and have a look at the demographic characteristics from the cohort of choice.
and then do the enrichment, scatter plot, single cell and network visualization.  <br><br>
It can also be run on your own device too. Only need to run "human app.R" or "mouse app.R" from the [drive folder](https://drive.google.com/drive/folders/1KvzT_USvOXiXbktSs-e3DCzwrqsfZDEl?usp=sharing).<br><br>
You are encouraged to cite this [paper](https://elifesciences.org/articles/88863) if you find it useful for your research. <br><br>
## Tutorial 
> [!NOTE]
> For mouse data is pretty much the same pipeline to use, just a few different. And be awarded that the color of the headboard for the mouse is green and the human one will be blue. 

> [!NOTE]
> All the images plus some of the data that generate them can be downloaded. The interactive images download buttons are on the top right and for other plots are from the bottom left. 

> [!IMPORTANT]
> You should hit the button again to update your settings. 
### Initial settings
First start from the dashboard. Select the sex that you want and input an [NBCI gene symbol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7494048/), or you will see a pop-up warning. <br><br>

![initial](https://github.com/mingqizh/GD-CAT/blob/main/images/hset.png) <br><br>
Hit the "process data" button after choosing a tissue. A progress bar will show at the bottom right of the screen. <br><br>
![pre](https://github.com/mingqizh/GD-CAT/blob/main/images/1689177951994.png) <br><br>
After the processing of the original data, you can see the pie chart where the gene is enriched in, age and sex of the cohort, cell type, and top genes. For the pie chart, it means how the gene you inputted from the selected tissue expressed across all tissue; and the `Original-removed` means how it looks like if the original tissue is removed. <br><br>
![pie](https://github.com/mingqizh/GD-CAT/blob/main/images/sfd.png) <br><br>
Click the `table` tab you can explore more about the data you are using. You can filter, rank, search, copy, or download the data in 2 formats.<br><br>
![table](https://github.com/mingqizh/GD-CAT/blob/main/images/table.png)<br><br>
### Enrichment Analysis
To generate the `GSEA` pathway analysis of the second part, users should click the pie chart body to select a tissue and threshold first. Your choice will be displayed then in this section. Then users can hit the `Start Analysis` button to see the results. After clicking the button, a progress bar will appear at the bottom right. <br><br>
![enrichment](https://github.com/mingqizh/GD-CAT/blob/main/images/enrich.png)<br><br>
### Scatter Plot
Same as the initial settings, make sure you put an official NBCI gene symbol here. <br><br>
![scatter](https://github.com/mingqizh/GD-CAT/blob/main/images/1689177689112.png)<br><br>
### Network Analysis
Move the slider bars to custom the gene numbers for within and peripheral tissue and then process the analysis. And you can also acquire the gene list generated the figure by clicking the `Download gene list` button at the bottom left corner. <br><br>
![net](https://github.com/mingqizh/GD-CAT/blob/main/images/net.png)


