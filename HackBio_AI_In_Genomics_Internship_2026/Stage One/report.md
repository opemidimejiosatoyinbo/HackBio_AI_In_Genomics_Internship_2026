# Exploratory Analysis of Drug Sensitivity Patterns in the Genomics of Drug Sensitivity in Cancer (GDSC) Dataset

**Project:** Stage One Analysis — AI for Genomics Internship  
**Dataset:** GDSC.xlsx
**Team** Glutamine
**Analysis Tools:** Python (pandas, seaborn, matplotlib, scipy)  
**Date:** May 2026

---

## Abstract

Cancer remains one of the most complex therapeutic challenges in modern medicine, in large part because no two tumours behave identically. The Genomics of Drug Sensitivity in Cancer (GDSC) project addresses this challenge by systematically profiling how hundreds of cancer cell lines respond to a large panel of drugs, and by linking those responses to the genomic, transcriptomic, and epigenomic characteristics of each cell line. This report presents a structured exploratory analysis of the GDSC dataset, covering data quality assessment, univariate and multivariate characterisation of drug sensitivity metrics, and an investigation into how molecular features such as copy number alterations (CNA), gene expression status, DNA methylation, and microsatellite instability (MSI) associate with differential drug response. The goal is not simply to describe numbers, but to interpret what these patterns mean for our understanding of cancer biology and the prospect of targeted therapy.

---

## 1. Introduction

Drug resistance is the single greatest obstacle to durable cancer treatment. A compound that works brilliantly against one tumour may be completely inert against another that, histologically, looks identical. Understanding why requires connecting pharmacological observations — how sensitive is this cell line to this drug? — to the underlying molecular state of the cell. The GDSC dataset was built precisely for this purpose: it provides a high-throughput pharmacogenomic resource that pairs drug sensitivity measurements with rich genomic annotations across a diverse library of human cancer cell lines.

The primary sensitivity metric in this dataset is the **IC50** — the drug concentration at which 50% of cell proliferation is inhibited. Because IC50 values span several orders of magnitude, a natural logarithm transformation (LN_IC50) is used throughout, which both stabilises variance and makes the distribution more amenable to parametric analyses. Complementary metrics include the **Area Under the Dose-Response Curve (AUC)**, which captures the overall breadth of drug activity, and the **Z-score**, which standardises LN_IC50 relative to the population of cell lines tested with the same drug. Together these three metrics provide a layered view of sensitivity: the IC50 tells us the concentration at which a drug starts to bite; the AUC tells us how effectively it bites across a range of concentrations; and the Z-score contextualises a cell line's response relative to all others.

---

## 2. Dataset Overview and Quality Assessment

### 2.1 Structure and Scale

The GDSC dataset in its raw form comprises **162,103 drug-cell line testing records** spanning **19 variables**. It contains measurements for **246 unique drugs** tested across **737 unique cancer cell lines**, identified by their COSMIC database IDs. The cell lines collectively represent **17 broad tissue systems** and a wide diversity of TCGA-labelled cancer types. The theoretical maximum number of drug-cell line combinations would be 181,302 (246 drugs × 737 cell lines), but the observed 162,103 entries confirm that not every drug was tested against every cell line — a common and biologically reasonable feature of large pharmacogenomic screens, where experimental resources are prioritised around mechanistic hypotheses.

### 2.2 Missing Data

A critical first step in any analysis is checking whether absences in the data are random or systematic. In this dataset, no missing values were detected across any of the 19 columns, which is a notably clean result for a large experimental compendium. This means that every drug-cell line record carries a complete set of sensitivity metrics and molecular annotations, allowing downstream analyses to proceed without imputation.

### 2.3 Handling Duplicate Records

Upon closer inspection, 4,290 combinations of cell line name and drug name appeared more than once. This is not trivially noise. A check across multiple annotation columns — MSI status, CNA, gene expression, methylation, screen medium, growth properties, and target pathway — confirmed that these duplicate records share identical metadata, pointing to repeated measurements rather than genuinely independent experiments. This is consistent with the known structure of the GDSC project, where some cell lines were screened in multiple experimental batches to improve robustness. The decision taken in the analysis pipeline was to resolve duplicates by averaging the continuous sensitivity metrics (LN_IC50, AUC, Z-score) per unique cell line-drug pair, while retaining a single copy of the associated molecular metadata. This produces a clean analytical table without discarding valid biological signal.

### 2.4 Data Standardisation

Several preprocessing steps were necessary before analysis. The LN_IC50 and AUC columns were coerced to numeric types, confirming no residual non-numeric entries. Binary annotation columns (CNA, Gene Expression, Methylation) were standardised to uppercase categorical variables to eliminate inconsistencies such as mixed-case entries. Column names were also stripped of whitespace and reformatted with underscores for computational consistency.

---

## 3. Exploratory Data Analysis

### 3.1 Distributions of Drug Sensitivity Metrics

The distributions of LN_IC50, AUC, and Z-score reveal important features of the pharmacological landscape captured in this dataset.

LN_IC50 is approximately normally distributed across all records, with a moderate spread around the mean and some right-skew reflecting a subset of cell line-drug combinations where very high concentrations are needed for any effect — a signature of intrinsic resistance. The near-symmetry of the bulk of this distribution suggests that, for most drugs, the tested cell lines span a genuinely informative range of sensitivities rather than clustering at one extreme.

AUC shows a slightly left-skewed distribution, indicating that while most drug-cell line interactions produce curves with meaningful dose-dependent suppression, a non-trivial fraction of combinations yield very shallow curves — cells that show dose-dependent but incomplete inhibition even at the highest tested concentrations. Biologically, this often reflects partial target engagement, drug efflux, or feedback loop rewiring.

Z-score distributions are centred close to zero, as expected for a standardised metric, but with heavy tails in both directions. The outliers on the sensitive end — large negative Z-scores — are of particular biological interest, as they represent cell lines that are dramatically more sensitive to a given drug than the average, which often indicates genetic dependency or addiction to the targeted pathway.

### 3.2 Cancer Type Representation

The dataset is not uniformly balanced across cancer types, which is an important caveat for any population-level inference. Solid tumour types — including lung, breast, and colorectal cancers — are well-represented, reflecting both their clinical prevalence and their prominence in pharmacogenomic research. Some rarer tumour types appear with fewer cell lines, which limits statistical power for those specific comparisons but does not undermine analyses conducted at the aggregate level.

### 3.3 Microsatellite Instability

Microsatellite instability status was characterised across cell lines, revealing a predominantly microsatellite-stable (MSS) population. This is consistent with the general prevalence of MSI across cancer types — high-MSI tumours account for a minority of most solid tumour types outside of mismatch repair-deficient colorectal and endometrial cancers. The boxplot analysis of LN_IC50 by MSI status hints at differences in sensitivity profiles between MSI-high and MSS cell lines, with MSI-high lines potentially showing distinct response patterns to certain agents. This is biologically plausible: the hypermutator phenotype associated with MSI can both sensitise tumours to DNA-damaging agents and alter the expression of drug targets.

---

## 4. Drug Sensitivity Patterns

### 4.1 Most Potent Drugs

Ranking drugs by their median LN_IC50 across all tested cell lines identifies compounds that suppress cancer cell growth at consistently low concentrations. The ten most potent drugs by this criterion represent agents with broad and powerful activity across diverse cellular backgrounds. Importantly, this ranking is presented alongside the standard deviation of LN_IC50 responses, because the distinction between a universally potent drug and a selectively potent one carries enormous translational relevance. A drug with a very low median LN_IC50 but a narrow standard deviation is a broadly cytotoxic agent — useful in some contexts, but potentially indiscriminate. A drug with a low median LN_IC50 and a wide standard deviation is more interesting from a precision medicine standpoint: it kills some cell lines at very low concentrations and barely touches others, which suggests its activity is gated by a specific molecular vulnerability.

At the other end of the spectrum, the least potent drugs — those with the highest median LN_IC50 values — do not necessarily represent failed compounds. Rather, they may reflect agents tested against cell line panels that are predominantly resistant to their mechanism of action, or compounds whose clinical efficacy depends on biomarker-selected patient populations that are not well-represented in this in vitro resource.

### 4.2 Variable Drug Responses Across Cell Lines

The drugs showing the greatest variability in LN_IC50 across cell lines are arguably the most informative from a biological standpoint. High variance indicates that a drug's efficacy is strongly context-dependent, meaning some molecular feature present in certain cell lines is either driving sensitivity or conferring resistance. These are precisely the drugs where biomarker discovery is most tractable. When a drug kills certain cell lines at nanomolar concentrations while leaving others unaffected at micromolar concentrations, it invites the question: what is genetically or epigenetically different about the sensitive lines? The answers to that question are the foundation of precision oncology.

### 4.3 Drug Sensitivity Across Cancer Types

The heatmap of median LN_IC50 across cancer types and drugs reveals a richly structured landscape. Rather than a homogeneous pattern of sensitivity, distinct cancer types cluster by their response profiles, and individual drugs show cancer-type-specific activity patterns. This is biologically expected: different cancers harbour different driver mutations and dependencies, and drugs targeting those specific vulnerabilities will show selectivity accordingly.

Some cancer types appear broadly sensitive across multiple drug classes, which may reflect high genomic instability, elevated proliferation rates, or deficiencies in DNA damage repair pathways. Others appear broadly resistant, which is more challenging to interpret without deeper molecular analysis — resistance can arise from drug target amplification, alternative pathway activation, drug metabolism, or simply from the presence of survival programmes not adequately captured by the tested drug panel.

The parallel cell-line-level heatmap (filtered to the 10 most variable drugs) reveals that even within a single cancer type, individual cell lines can diverge dramatically in their sensitivity profiles. This intra-type heterogeneity is a critical reminder that cancer type alone is an insufficient guide to treatment selection; molecular profiling of individual tumours is necessary.

---

## 5. Genomic, Transcriptomic, and Epigenomic Correlates of Drug Response

### 5.1 Analytical Framework

To connect molecular biology to drug pharmacology, the analysis examined how three categories of molecular annotation — copy number alterations (CNA), gene expression status, and DNA methylation — associate with LN_IC50, AUC, and Z-score. For each drug, the cell lines were split by their status for each molecular feature (annotated as binary flags), and an independent t-test was used to assess whether the two groups differ significantly in drug sensitivity. The effect size (mean difference in LN_IC50 between groups) and the resulting p-value were used to construct volcano plots.

### 5.2 Copy Number Alterations

Copy number alterations reflect regions of the genome that have been amplified or deleted in cancer cells, often with direct consequences for protein dosage and cellular signalling. Cell lines annotated as harbouring CNA show differential drug responses compared to those without, and the direction and magnitude of this effect varies substantially by drug. This is consistent with the known biology of CNA-driven oncogene amplification: for example, EGFR-amplified tumours show heightened sensitivity to EGFR inhibitors, while amplification of drug efflux transporters can confer resistance. The annotation in this dataset is a binary flag rather than a gene-level resolution, so the analysis captures global CNA burden rather than specific alterations, which limits interpretation but still reveals a meaningful signal.

### 5.3 Gene Expression

Gene expression status shows associations with drug response that, in several cases, reach statistical significance. Transcriptomic differences between sensitive and resistant cell lines often reflect the activity of the signalling pathways that drugs target. A cell line that is transcriptionally addicted to a given pathway will, in theory, be more sensitive to inhibitors of that pathway. Identifying these transcriptional correlates of sensitivity is therefore one of the key goals of pharmacogenomic studies, and the present analysis confirms that the gene expression axis carries meaningful predictive information across multiple drug classes.

### 5.4 DNA Methylation

DNA methylation is an epigenetic mechanism that can silence tumour suppressor genes, alter chromatin accessibility, and reshape the transcriptional landscape of a cancer cell without changing the underlying DNA sequence. Its association with drug sensitivity is therefore indirect but real: by silencing genes involved in drug uptake, apoptosis induction, or DNA repair, methylation can modulate how a cell responds to chemotherapy or targeted agents. The observed associations between methylation status and LN_IC50 in the volcano plot analysis suggest that epigenetic state contributes meaningfully to drug response variation in this dataset, though the binary flag annotation again limits resolution.

### 5.5 Volcano Plot Interpretation

The volcano plots integrate effect size and statistical significance into a single visual. Points in the upper-right quadrant represent drug-feature combinations where the presence of a molecular alteration is associated with significantly higher LN_IC50 — i.e., resistance. Points in the upper-left represent cases where the molecular alteration is associated with significantly lower LN_IC50 — i.e., sensitivity. The density of significant associations varies across the three molecular features, with CNA and gene expression yielding more statistically significant hits than methylation alone. This likely reflects both the greater dynamic range of CNA and expression effects and the relative coarseness of the binary methylation annotation.

---

## 6. Discussion

This analysis makes several things clear that pure summary statistics cannot convey.

First, drug sensitivity in cancer is not a simple function of cancer type. While broad patterns are visible at the level of TCGA cancer labels, the within-type variation is substantial. This means that clinical decisions guided only by histological diagnosis will miss important opportunities for precision targeting. The GDSC data supports a vision of oncology in which molecular profiling guides drug selection — a vision that is increasingly being realised in clinical practice through companion diagnostics and basket trials.

Second, the most pharmacologically interesting drugs are not necessarily the most broadly potent ones. From a precision medicine perspective, a drug with high variance in its response profile is a candidate for biomarker-driven therapy. The challenge is to identify what distinguishes the sensitive from the resistant lines, and this is precisely the kind of question that pharmacogenomic datasets like GDSC are designed to address.

Third, genomic, transcriptomic, and epigenomic features all contribute independently to drug response variation. This means that resistance is rarely attributable to a single molecular mechanism. In practice, the clinical development of targeted therapies requires understanding not just the primary driver of sensitivity, but the secondary mechanisms that allow tumours to escape treatment. The associations uncovered here are a starting point for such mechanistic investigations.

Finally, the dataset's design — multiple drugs per cell line, multiple cell lines per cancer type — enables analyses at multiple levels of biological resolution. The richness of this resource is its greatest asset, but it also means that naive aggregate analyses can obscure important biology. Stratified analyses by cancer type, specific drug class, or molecular subgroup are necessary to translate the high-level patterns described here into actionable biological insights.

---

## 7. Limitations

Several limitations should be acknowledged. The GDSC dataset is generated from in vitro cancer cell lines, which are imperfect models of tumour biology in vivo. Cell lines may accumulate culture-acquired mutations, lose certain epigenetic features of the primary tumour, and lack the stromal and immune microenvironment that influences drug response in patients. The molecular annotations used as predictors are represented as binary flags in this analysis, which discards important quantitative information about the degree of alteration. More powerful analyses would use continuous measures of gene expression or specific gene-level CNA calls rather than binary indicators. The statistical comparisons performed here have not been corrected for multiple testing across 246 drugs and three molecular features, meaning some significant findings will be false positives and should be treated as hypothesis-generating rather than confirmatory.

---

## 8. Conclusions

The GDSC dataset provides a uniquely comprehensive window into the pharmacogenomics of cancer. This exploratory analysis has established that the dataset is largely complete and well-structured; that drug sensitivity varies substantially across both cancer types and individual cell lines; that certain drugs show high-variance response profiles that make them strong candidates for biomarker-driven therapy; and that genomic, transcriptomic, and epigenomic features all carry meaningful information about drug sensitivity. These findings lay the groundwork for more targeted analyses — including machine learning models for sensitivity prediction, pathway-level enrichment analyses of sensitive versus resistant cell lines, and integration of drug mechanism-of-action annotations to explain the patterns observed here.

Understanding why some cancer cells die in response to a drug while others survive is, ultimately, a question about the molecular circuitry of the cell. Pharmacogenomic datasets like GDSC bring us closer to an answer — and closer to a future in which cancer treatment is matched precisely to the biology of each individual tumour.

---

## References

1. Garnett MJ et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. *Nature*. 2012;483(7391):570-575.
2. Yang W et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Research*. 2013;41(D1):D955-D961.
3. Iorio F et al. A Landscape of Pharmacogenomic Interactions in Cancer. *Cell*. 2016;166(3):740-754.
4. Cancer Genome Atlas Research Network. Comprehensive genomic characterisation of squamous cell lung cancers. *Nature*. 2012;489(7417):519-525.