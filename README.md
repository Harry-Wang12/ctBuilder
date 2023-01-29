# ctBuilder: A framework for building pathway crosstalks by combining single cell data with bulk cell data

## Abstract

### Background
Biomedical research scientists routinely use publicly available gene regulatory pathway databases (e.g., KEGG, Reactome) to extract interpretation from high-throughput multi-omics studies.  However, deficiencies of these existing pathway resources include: (i) the databases store “individualized” pathways as a collection of networks thus making discovery of crosstalk across the boundaries of curated pathways difficult; and (ii) curated pathways themselves are incomplete, particularly from the disease perspective. 
### Method
We propose a computational pathway extension framework, called ctBuilder, which aims to tackle the deficiencies. First, it explores if “individually” curated pathways can be meaningfully extended using curated gene regulatory networks (GRNs) (e.g., BioGrid, STRING, TRRUST) which are another type of public resources available for omics data analysis. Second, it uses single cell gene expression data to identify potential candidates to tune the pathway extension “context-specific”. Third, it uses bulk cell gene expression data sets to further refine and improve likelihoods of the regulatory relationships estimated to interlink constituents involved between pathways. We claim the resulting extension incorporates crosstalk between pathways that is “disease specific” when the extension is driven exclusively by a “cohort” of gene expression data sets from a specific disease. For demonstration, we applied ctBuilder to 6 non-alcoholic steatohepatis related studies including 1 single cell dataset and 4 bulk cell datasets publicly available from GEO. 
### Result and Conclusion
We show how three pairs of known pathways, TNF signaling vs. NF-kB signaling, Toll-like signaling vs. NF-kB signaling, NF-kB signaling vs Insulin signaling, can be extended to include crosstalk subnetworks specifically related to the liver disease. Statistical significance of the expended network is shown through null distribution test statistics. Most of newly identified genes in the crosstalk paths are already known in the liver disease literature. Discovered relationships among them merit wet-lab experiments for validation.


## Questions
