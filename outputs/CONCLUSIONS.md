# Conclusions

This analysis examined **82 cleaned IC50 bioactivity records** for **PHD finger protein 23 (CHEMBL2424508)** retrieved from ChEMBL.

## Bioactivity profile

- Active compounds (IC50 <= 1,000 nM): **0**
- Intermediate compounds (1,000 < IC50 < 10,000 nM): **0**
- Inactive compounds (IC50 >= 10,000 nM): **82**

## Interpretation

The most common bioactivity class is **inactive**, representing **82/82 (100.0%)** of the processed observations.

The observed pIC50 values range from **4.52** to **5.00**, with a median of **5.00**.

## Main limitation

Only one bioactivity class is represented after preprocessing. Therefore, an active-vs-inactive classification model would not be scientifically meaningful for this data slice. This is a limitation of the available target data rather than a reason to manufacture or relabel observations.

## Overall conclusion

The workflow successfully converts raw ChEMBL measurements into a reproducible dataset containing activity classes, pIC50 values and molecular descriptors, together with plots that summarize the chemical and bioactivity space. The amount and diversity of available PHF23 bioactivity data should determine whether predictive modeling is justified. Future work could expand the dataset with additional compatible activity measurements or apply the same workflow to a target with a larger and better-balanced ChEMBL dataset.

## Generated plots

Plots are saved in `outputs/plots/` and include the bioactivity-class distribution, pIC50 distribution, molecular-weight/pIC50 relationship, LogP/pIC50 relationship and descriptor correlation matrix.