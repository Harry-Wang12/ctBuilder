# ctBuilder: A framework for building pathway crosstalks by combining single cell data with bulk cell data using feature selection

## Abstract

### Background
Biomedical research scientists routinely use publicly available gene regulatory pathway databases (e.g., KEGG, Reactome) to extract interpretation from high-throughput multi-omics studies.  However, deficiencies of these existing pathway resources include: (i) the databases store “individualized” pathways as a collection of networks thus making discovery of crosstalk across the boundaries of curated pathways difficult; and (ii) curated pathways themselves are incomplete, particularly from the disease perspective. 
### Method
We propose a computational pathway extension framework, called ctBuilder (Fig. 1.), which aims to tackle the deficiencies. First, it explores if “individually” curated pathways can be meaningfully extended using curated gene regulatory networks (GRNs) (e.g., BioGrid, STRING, TRRUST) which are another type of public resources available for omics data analysis. Second, it uses single cell gene expression data to identify potential candidates by tuning the pathway extension “context-specific” with ReliefF feature selection method. Third, it uses bulk cell gene expression data sets to further refine and improve likelihoods of the regulatory relationships estimated to interlink constituents involved between pathways. We claim the resulting extension incorporates crosstalk between pathways that is “disease specific” when the extension is driven exclusively by a “cohort” of gene expression data sets from a specific disease. For demonstration, we applied ctBuilder to 6 non-alcoholic steatohepatis related studies including 1 single cell dataset and 4 bulk cell datasets publicly available from GEO.

![image](https://user-images.githubusercontent.com/114254986/215357548-27d7261d-c6e1-4980-a2f8-183ef496e44b.png)<br>
*Fig. 1. The overall process of ctBuilder*


### Result and Conclusion
We show how three pairs of known pathways, TNF signaling vs. NF-kB signaling, Toll-like signaling vs. NF-kB signaling, NF-kB signaling vs Insulin signaling, can be extended to include crosstalk subnetworks specifically related to the liver disease (Fig. 2, 3, 4). Statistical significance of the expended network is shown through null distribution test statistics. Most of newly identified genes in the crosstalk paths are already known in the liver disease literature. Discovered relationships among them merit wet-lab experiments for validation.

![image](https://user-images.githubusercontent.com/114254986/215357583-6233275c-7b4d-4761-8441-d65f9ada100b.png)<br>
*Fig. 2. Crosstalk routes from TNF signaling pathway to NF-kB signaling pathway*

![image](https://user-images.githubusercontent.com/114254986/215357604-6e37091f-9f83-40cb-b96f-1ea1eb945e1b.png)<br>
*Fig. 3. Crosstalk routes from TNF signaling pathway to NF-kB signaling pathway*

![image](https://user-images.githubusercontent.com/114254986/215357608-d9bfc749-bd26-4018-a8cd-82f819354a62.png)<br>
*Fig. 4. Crosstalk routes from TNF signaling pathway to NF-kB signaling pathway*


## Questions
<p><em>For more information please connecting honglin.wang@uconn.edu</em></p>
<p><em>For citation please check the <a href="https://ieeexplore.ieee.org/document/9669811">here</a></p>
