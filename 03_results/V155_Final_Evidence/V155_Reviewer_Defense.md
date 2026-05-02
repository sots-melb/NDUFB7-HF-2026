### Reviewer Defense Framework (Pre-emptive)

**Q1: "Why focus on ferroptosis when apoptosis/necroptosis show stronger correlation?"**
A: We explicitly tested this (V151, n=313). While apoptosis and necroptosis are stronger statistical correlates, they lack: (1) a metabolic ratio biomarker (ACSL4/GPX4) that can be measured in patient serum; (2) a druggable defense system with clinical-grade inhibitors (Fer-1, Liproxstatin-1); and (3) temporal anchoring to the NDUFB7 breakpoint in pseudotime analysis. We frame ferroptosis not as the exclusive mechanism but as the most **actionable node** within a broader death-signature cascade. This honest framing strengthens clinical translation.

**Q2: "Is the bimodal distribution real or an artifact?"**
A: Three independent lines support bimodality: (1) ZI-GMM on GSE183852 CM cells (BIC favors G=2 over G=1 by >200 units); (2) Dip test rejects unimodality (p<0.001); (3) the low-NDUFB7 subpopulation is enriched in DCM vs control (Wilcoxon p<0.001) and correlates with disease stage. Zero-inflation alone cannot explain the second peak.

**Q3: "Is NDUFB7 loss a cause or consequence of HF?"**
A: MR provides genetic evidence for causality (OR=0.73, p=0.028). While bulk transcriptomics cannot distinguish cause from consequence, the temporal synchronization in pseudotime (NDUFB7 breakpoint precedes ferroptosis driver upregulation) and the spatial gradient (highest loss in IZ, lowest in FZ) support a causal/driver role rather than passive downregulation.

**Q4: "Why no wet-lab validation?"**
A: This is a discovery-phase bioinformatics study. We have identified Ferrostatin-1 (GSE243655) as a pharmacological proxy for ferroptosis inhibition in cardiomyocytes, and our drug repurposing pipeline (clue.io + DrugReflector) has generated testable candidates. Wet-lab validation is planned for Revision/Phase II.

**Q5: "Is the pan-cell-death finding a weakness?"**
A: No. The finding that NDUFB7 loss coordinates multiple death pathways is biologically realistic for a mitochondrial OXPHOS defect. Mitochondria are the convergence point of apoptosis (cytochrome c release), necroptosis (ROS burst), and ferroptosis (iron-dependent lipid peroxidation). Our paper's contribution is identifying the **most druggable intersection point** of this convergence.
