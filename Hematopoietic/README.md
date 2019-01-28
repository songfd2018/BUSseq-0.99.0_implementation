# Generate the Figures of Hematopoietic Study

This folder contains scripts that generate results for the mouse hematopoietic study, including Fig 5, in the manuscript. To reproduce all the results, please

1. Run data_preprocessing.R to:
	- conduct quality control on the two public hematopoietic scRNA-seq datasets from the Gene Expression Omnibus;
	- select the shared highly variable genes;
	- save observed read count matrix in the "./RawCountData" folder.

2. Run run_BUSseq.R to apply BUSseq to the hematopoietic dataset. The script runs BUSseq with different numbers of cell types and plots the BIC plots.

3. Run run_MNN.R to apply MNN to the hematopoietic dataset.

4. Run run_Seurat.R to apply Seurat to the hematopoietic dataset.

5. Run run_ZINBWaVE.R to apply ZINB-WaVE to the hematopoietic dataset

6. Run method_evaluation.R to:
	- calculate the ARI and Silhouette coefficients and save them in the "./Results" folder; 
	- generate t-SNE plots for Fig 5 and save them in the "./Image/tSNE" folder.
	- draw the violin plot of Silhouette coefficients and save it in the "./Image/Other" folder.
