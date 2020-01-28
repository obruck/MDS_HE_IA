# MDS_HE_IA
Image analysis of H&amp;E stained slides of patients with myelodysplastic syndrome (MDS) using convolutional neural networks and elastic net.

Visual inspection of Hematoxylin&Eosin (H&E) stained tissue samples is a standard procedure for histopathological diagnosis. We collected formalin-fized paraffin-embedded (FFPE) bone marrow (BM) trephine biopsies from 236 MDS patients, 87 MDS/myeloproliferative neoplasia (MPN) patients and 10 healthy controls. These were transformed into tissue microarrays (TMA) and stained with a standard in-house H&E staining. Each image was split into 500 equally-sized tiles and colors converted into grayscale for more robust analysis. Only non-white (= non-empty) tiles were retained for image analysis. Visual features were extracted at the tile level with VGG16 and Xception networks using imagenet weights. In addition, we analyzed each TMA image at the pixel-level using a trained Weka model and at the cell level using Qupath. Hence, we are able to associate high-resolution image analysis data from CNNs with more comprehensible image analysis features from pixel and cell level image analysis.

Visual features from the CNNs were fed into several elastic net model with predefined alpha penalization (0.0, 0.5 and 1.0) and lambda values were optimized for each. The dataset was divided at the TMA level into training (2/3) and test (1/3) datasets. In total 3 different penalization from 2 different CNN models were used to predict
- overall survival
- risk of AML
- treatment response to azacitidine
- diagnosis (between MDS and MDS/MPN)
- mutation status for TP53, ASXL1, DNMT3A, SRSF2, TET2, RUNX1, SF3B1, NRAS/KRAS, IDH1/IDH2, STAG2
- mutations in cell cycle, cell differentiation, DNA chromatin structurea and spliceosome regulation pathways
- chromosome 3, 5q, 7q, 7, 8, 20q, complex karyotype
- patient age and gender
- de novo vs. secondary MDS
- ipssr score and ipssr cytopenia score
Only visual features from MDS patients were included for these analyses except for diagnosis prediction were also MDS/MPN patients included.

As visual features from the Xception network performed superiorly to VGG16, we used these features to analyze the natural embedding of MDS, MDS/MPN and control subjects in an unsupervised fashion. For this purpose, we used PCA and UMAP and clustered data with Phenograph and Kmeans clustering. Due to the non-parametric nature of visual data, we observed UMAP and Kmeans to perform more accurately than other combinations.

This repository included scripts to
1. preprocess image data (Python)
2. extract visual features (R)
3. build elastic net models (R)
4. generate AUC of prediction model ROCs and plot AUC data + generate most representative images for each prediction models to assess models (R)
5. combine tile-level prediction models with pixel-level and cell-level image analysis data (R)


Scripts have been by Oscar Brück.


Oscar Brück, MD  
Hematology Research Unit Helsinki, University of Helsinki &  
Data Admnistration, Helsinki University Hospital &  
Hematology, Helsinki University Hospital  
Helsinki Finland  
oscar.bruck@helsinki.fi
