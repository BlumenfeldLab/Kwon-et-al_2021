# Kwon-et-al_2021
This includes the code files used for the analysis **"Early Cortical Signals in Visual Stimulus Detection"** - Kwon et al., 2021.

## Summary
A detailed description of all methods used in the paper is given in the Methods section of the manuscript. If you have further questions regarding this paper, we encourage you to contact Dr. Hal Blumenfeld(hal.blumenfeld@yale.edu)

## Software
1.	GammaPowerExtraction_preprocessing.m – Including following steps: extracting icEEG signal and information from raw data, artifact rejection and z-scored gamma power calculation.
2.	MappingGammaPower2Surface.m - Mapping z-scored gamma power to the brain surface.
3.	PipelineForPermutationTest.m - Cluster-based statistics: clust_perm1.m function in matlabmk Toolbox ( http://kutaslab.ucsd.edu/matlabmk_fn_docs ) was modified for the icEEG data set.
4.	ROIAnalysis.m 

## Data
Data used in this study are from the free-recall task in the RAM project, and are publicly available from http://memory.psych.upenn.edu/Data_Request.
