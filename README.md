# Uncovering the embodied dimension of the wandering mind (Body-Wandering)

Leah Banellis, Niia Nikolova, Malthe Brændholdt, Melina Vejlø, Francesca Fardo, Jonathan Smallwood, Micah G. Allen

---

## Abstract

When at rest, the mind becomes preoccupied with self-generated thoughts, commonly known as mind-wandering. While the social, autobiographical, and temporal features of these thoughts have been extensively studied, little is known about how frequently the wandering mind turns towards the interoceptive and somatic body. To map this under-explored component of "body-wandering," we conducted a large-scale neuroimaging study in 536 healthy participants, expanding a retrospective multidimensional experience sampling approach to include probes targeting visceral and somatomotor thoughts. Our findings reveal a robust inter-individual dimension of body-wandering characterized by negative affect, high autonomic arousal, and a reduction in socially oriented thoughts. Despite this negative tone, individual differences in the propensity for body-wandering were associated with lower self-reported symptoms of ADHD and depression. Multivariate functional connectivity analyses further revealed that affective, body-oriented thoughts are related to a pattern of thalamocortical connectivity interlinking somatomotor and interoceptive-allostatic cortical networks. Collectively, these results demonstrate that self-generated thoughts exhibit core embodied features which are linked to the ongoing physical and emotional milieu of the visceral body.

---

## Preprint

Access the preprint here: [https://www.biorxiv.org/content/10.1101/2024.10.25.620252v2.full](https://www.biorxiv.org/content/10.1101/2024.10.25.620252v2.full)

---

## Citation

> Banellis, L., Nikolova, N., Brændholdt, M., Vejlø, M., Fardo, F., Smallwood, J., & Allen, M. G. (2024). Body-wandering reveals an embodied dimension of thought with distinct affective and neural signatures. *bioRxiv*, 2024-10.

## Figures

![Body-Wandering & Cognitive-Wandering Relationships](figures/Final_Figures/Fig3_EFA.png)



![Cortical Fingerprint of Embodied Mind-Wandering](figures/Final_Figures/Fig4_NeuralCCA.png)

# Scripts for Body-Wandering Analyses

Note: data & scripts for raw/preprocessed fMRI, raw physio, and raw survey individual subject data stored in Visceral Mind Project BIDS structure:  

**1. Body/Mind-Wandering Item-level Analyses**:  
    /BodyWanderingCCA/scripts/body_wandering_itemplots.Rmd

**2. Mental Health and Physiological Arousal Relationships with Body/Cognitive-Wandering**:  
    /BodyWanderingCCA/scripts/ItemAnalyses_MHPhysio.Rmd

**3. Exploratory Factor Analysis of Body/Mind-Wandering Items**:  
    /BodyWanderingCCA/scripts/MDIES_validity_reliability.Rmd

--------------------------------------------------------------------------------------------------------

**4. Cortical Fingerprinting of Body-Wandering via CCA**:

    4a. Create matched input matrices for CCA (X = Brain Connectivity, Y = Body/Mind-Wandering Items, C = Confounds):  
    /BodyWanderingCCA/scripts/CCA/CCAprepData.m

    4b. X.mat, Y.mat & C.mat from (4a.) saved in:  
    '/BodyWanderingCCA/scripts/CCA/cca_pls_toolkit-master/_Project_BodyWandering/data/', 
    and framework folder created:  
    '/BodyWanderingCCA/scripts/CCA/cca_pls_toolkit-master/_Project_BodyWandering/framework/'  

    4c. Run CCA with CCA/PLS toolkit (run on cluster via BodyWanderingCCA/scripts/CCA/cca_pls_toolkit-master/cca_jobs_slurm.sh):  
    BodyWanderingCCA/scripts/CCA/cca_pls_toolkit-master/RunCCA.m

    4d. Plot CCA result (variate scatterplot, body/mind-wandering loading barplot & wordcloud, extraction of cca mode result details):  
    BodyWanderingCCA/scripts/CCA/cca_pls_toolkit-master/CCA_plots.m

    4e. Plot CCA result: connectivity loadings via Chordplot:
    /BodyWanderingCCA/scripts/CCA/CCA_plots.Rmd
     
    4f: Plot CCA result: connectivity loading sums projected on Schaefer-Subcortical216 parcellated brain (save nifti to plot with MRIcroGL):  
    /BodyWanderingCCA/scripts/CCA/CCA_Brain_Figures_abend.ipynb

--------------------------------------------------------------------------------------------------------

**5. Gradient Association of CCA Result**:

    5a. Create gradient-averages per Schaefer-Subcortical216 parcellated brain regions (store parcel numbers & labels):
    /BodyWanderingCCA/scripts/Gradients/gradients.ipynb

    5b. Convert Schaefer-Subcortical216 atlas to surface space using the freesurfer function ‘mri_vol2surf’ :
    /BodyWanderingCCA/scripts/Gradients/gradient2surface.sh

    5c. Gradient averages for top 1% regions from CCA & Spin test for creating null distribution (Need toolbox: https://github.com/spin-test/spin-test and FreeSurfer: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall):
    /BodyWanderingCCA/scripts/Gradients/SpinPermuFS_LB.m

    5d. Statistics and plots:
    /BodyWanderingCCA/scripts/Gradients/SpinTestResults_AverageTopCCA.m

--------------------------------------------------------------------------------------------------------

**6. Convergent Validity of EFA-related connectivity and CCA connectivity**:  
    /BodyWanderingCCA/scripts/ControlAnalyses/FactorConnectivity.ipynb
    (& end section of /BodyWanderingCCA/scripts/MDIES_validity_reliability.Rmd)

--------------------------------------------------------------------------------------------------------

**7. Additional analyses and control analyses**:

  6. Additional analyses and control analyses:
  ~/Git/BodyWanderingCCA/scripts/ControlAnalyses
  
--------------------------------------------------------------------------------------------------------

### JASP statistics:

Multiple Linear Regressions of Body/Cognitive-Wandering Average Relationships, Patial Correlations, and Wilcoxon Tests (had to hide currently due to protected mental health data - will be released with Visceral Mind Project): 
    /BodyWanderingCCA/data/JASP/BodyWanderingStats.jasp
  
