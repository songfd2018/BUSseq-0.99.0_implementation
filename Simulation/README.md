# Generate the Figures of Simulation Study

This folder contains scripts that generate results for the simulation study, including Fig 3 and 4, in the manuscript. To reproduce all the results, please

1. Run simulate_data.R to:
	- generate simulated read count matrices for all four bacthes and save them one by one in the "./RawCountData" folder;
	- generate the heatmap of the true parameters, the logarithm of the underlying true expression levels and the observed read counts and save them in the "./Image/Heatmap" folder;
	- store the true cell type indicators pi, which will be used to  calculate the ARI, and the cell-specific size factors delta, which will be used to draw scatter plot, in the "./True_para" folder.

2. Run run_BUSseq.R to apply BUSseq to the simulated dataset. The script runs BUSseq with different numbers of cell types and plots the BIC plots.

3. Run draw_heatmap.R to:
	- plot the heatmap of the estimated parameters, the logarithm of the estimated underlying true expressoin levels, and the logarithm of the corrected read counts and save them in the "./Image/Heatmap" folder;  
	- generate the scatter plot for the estimated versus true cell-specific size factors and save them in the "./Image/Other" folder.

4. Run run_MNN.R to apply MNN to the simulated dataset.

5. Run run_Seurat.R to apply Seurat to the simulated dataset.

6. Run run_ZINBWaVE.R to apply ZINB-WaVE to the simulated dataset

7. Run method_evaluation.R to:
	- calculate the ARI and Silhouette coefficients and save them in the "./Results" folder; 
	- generate t-SNE plots for Fig 4 and save them in the "./Image/tSNE" folder.
	- draw the violin plot of Silhouette coefficients and save it in the "./Image/Other" folder.