# Discussion (V62 Complete, ~1200 words)

## Spatial and Anatomical Specificity of NDUFB7

Our multi-platform analysis reveals that NDUFB7 is not a generic heart failure biomarker but a spatially restricted, anatomically variable, and disease-sensitive mitochondrial subunit. Five independent platforms (n=431 total) show no robust overall downregulation in heart failure (p=0.44-0.52), yet spatial transcriptomics identifies fibrotic-zone-specific depletion (FZ 40.2% vs IZ 64.9% positivity). This apparent paradox resolves when considering anatomical location and cell-type composition: normal left atrium shows higher NDUFB7 in non-cardiomyocytes (NCM median=4 vs CM median=1, GSE109816, n=9,994), while left ventricle shows cardiomyocyte-dominant expression (CM>FB, p=9.9×10⁻⁹⁵, GSE183852, n=220,752). Disease further attenuates this pattern (Donor CM mean=0.205 vs DCM CM mean=0.140, p=5×10⁻²⁵). These findings suggest that bulk-averaged analyses mask cell-type and location-specific dynamics, and that NDUFB7's role in heart failure is niche-dependent rather than universal.

The extreme zero-inflation observed across all datasets (87-92% zero values in single-cell/nucleus data, median=0) indicates that NDUFB7 is not a constitutive "housekeeping" subunit but rather a low-abundance, high-specificity accessory protein. This sparse expression profile may reflect metabolic heterogeneity within cell populations: only a subset of cells at a given time point require full Complex I assembly, and NDUFB7 may be dynamically regulated in response to local energy demands.

## Tissue-Specific Genetic Regulation

The opposite directional effects of blood versus heart eQTLs (eQTLGen rs11085898: T allele increases NDUFB7; GTEx Heart-LV rs8103021: T allele decreases NDUFB7) highlight the complexity of cis-regulatory architecture. This tissue discordance cautions against inferring cardiac causal effects from blood-based genetic instruments—a common practice in Mendelian randomization studies that may have contributed to inconsistent findings in the mitochondrial genetics literature. The heart-specific eQTL (GTEx) is particularly relevant given the cardiac phenotype, and the negative slope (slope = -0.080) suggests that the T allele reduces NDUFB7 expression in cardiac tissue, consistent with our observed disease-related attenuation.

Future cardiac-specific GWAS and eQTL resources with larger sample sizes (n>1,000) will be necessary to establish definitive genetic causality and to identify additional regulatory variants that may explain the inter-patient heterogeneity observed in our bulk validation cohorts.

## Mechanistic Implications and Therapeutic Hypotheses

In silico knockout predictions (PROGENy) suggest that NDUFB7 attenuation impairs OXPHOS (p=1.13×10⁻¹³), increases ROS (p=6.06×10⁻⁴), and activates fibrotic pathways (p=3.94×10⁻⁶). The iron death hypothesis—linking NDUFB7 loss to SLC25A28-mediated iron transport and subsequent ferroptosis—emerges from our correlation analysis (GSE57338, Spearman rho [待填入], p [待填入]) and provides a mechanistic bridge between mitochondrial dysfunction and cell death. However, we emphasize that all in silico predictions are hypotheses requiring experimental validation.

The therapeutic implications are twofold. First, niche-specific mitochondrial rescue—targeting NDUFB7 restoration specifically within the fibrotic zone—may be more effective than systemic Complex I activation, which could disrupt energy homeostasis in healthy tissue. Second, the "anatomical specificity" finding suggests that therapeutic strategies may need to be tailored to chamber-specific pathophysiology (atrial vs ventricular).

## Limitations

First, our spatial autocorrelation analysis confirms pseudoreplication in Visium data, and our conservative LMM approach may underestimate true spatial effects. Alternative spatial models (e.g., CAR, SAR) will be explored in future work. Second, anatomical specificity is based on two independent datasets (left atrium and left ventricle); multi-chamber sampling (right atrium, right ventricle, septum) is warranted. Third, genetic causal inference is limited by the lack of heart-failure-specific GWAS with sufficient power and by the tissue mismatch between blood eQTLs and cardiac phenotypes. Fourth, all mechanistic predictions are computational and require experimental validation (e.g., siNDUFB7 in H9C2 cardiomyocytes, AAV9-NDUFB7 in MI mouse models). Finally, the therapeutic hypotheses (coenzyme Q10, MitoQ, ferrostatin-1) are preliminary and require preclinical efficacy and safety testing.

## Conclusion

NDUFB7 emerges as a microenvironment-sensitive accessory subunit of Complex I, whose spatially restricted loss in myocardial infarction fibrotic zones reflects anatomical location, cell-type composition, and disease stage. Rather than a universal heart failure marker, NDUFB7 represents a precision-medicine target for niche-specific mitochondrial rescue, with implications for both diagnostic stratification and therapeutic development.
