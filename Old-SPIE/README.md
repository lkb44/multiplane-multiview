# Multi-Plane and Multi-Region Radiomic Features

## Conference Paper

The features and code in this repository were used to for the following 2022 SPIE submission: *Integrating Multi-Plane and Multi-Region Radiomic Features to Predict Pathologic Response to Neoadjuvant Chemoradiation in Rectal Cancers via Pre-Treatment MRI*

The conference submission is still under review.

Please read the 4-page abstract submission for details on background/purpose,  experimental design, and final results.

The goal is to determine whether multi-plane, multi-region shape and texture radiomic features are associated with pathological complete response (pCR).

## Dataset Description

- Pre-treatment T2-weighted axial and coronal scans from UH and VA
- The axial scans had the following regions annotated [*note, the label number of region in the masks are provided*]:
    - Tumor (label = 1)
    - Lumen (label = 2)
    - Perirectal Fat (label = 4)
  - The coronal scans had the following regions annotated [*note, the label number of region in the masks are provided*]:
    - Rectal Wall (label = 7, 8)
    - Lumen (label = 2)
    - Perirectal Fat (label = 4)
- Limited analyses to 2D slices with:
    - Largest tumor area on axial scans
    - Largest rectal wall area on coronal scans

## Experimental Design

- We extracted either 3D or 2D shape and texture features from a total of 10 regions across axial and coronal views:


| **Feature Family** | **View** |         **ROIs**        |
|:------------------:|:--------:|:-----------------------:|
|   Shape Features   |   Axial  |    Lumen, Tumor, Fat    |
|                    |  Coronal | Lumen, Rectal Wall, Fat |
|  Texture Features  |   Axial  |        Tumor, Fat       |
|                    |  Coronal |     Rectal Wall, Fat    |

*Note, the rectal lumen is hollow. Hence, it does not make sense to extract texture features from the lumen...because there is no texture*

- Using the shape and texture radiomic features, we classified each patient as pCR or non-pCR (i.e., partial/incomplete pathological response to treatment). To do this binary classification, we divided patients based on their pathological T-stage (ypT):
  - T0-T2 = pCR = Class 0
  - T3-T4 = non-pCR = Class 1

- For all experiments, feature selection was performed on the training cohort using 50 iterations of 5-fold cross validation. Wilcoxon Ranksum or MRMR was used as the feature selection scheme.

- For all experiments, QDA or Random Forest was used as the classifier. The classifier was evaluated on the training cohort and holdout testing cohort.

## Code Description
