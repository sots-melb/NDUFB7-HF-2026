### NDUFB7 Depletion and Pan-Cell-Death Signature (Revised Narrative)

Bulk transcriptomic analysis (GSE57338, n=313) revealed that NDUFB7 loss was significantly associated with multiple cell-death pathways, forming a **pan-cell-death signature**:

| Pathway | Spearman ρ | P-value | Interpretation |
|:---|:---|:---|:---|
| Apoptosis | −0.361 | 6.4×10⁻¹¹ | Execution cascade activation |
| Necroptosis | −0.333 | 2.0×10⁻⁹ | Necrotic cell death priming |
| Ferroptosis Defense | −0.248 | 9.5×10⁻⁶ | **Anti-ferroptosis system collapse** |
| Autophagy | −0.237 | 2.5×10⁻⁵ | Quality-control failure |
| Pyroptosis | −0.219 | 9.5×10⁻⁵ | Inflammatory death priming |
| Ferroptosis Execution | −0.102 | 7.2×10⁻² | Not significant |

**Key insight**: Rather than a ferroptosis-specific effect, NDUFB7 depletion triggers a **coordinated collapse of multiple death-defense systems**. The ferroptosis axis is highlighted because (i) its defense module (GPX4-SLC7A11-FTH1) is the most pharmacologically targetable among the significant pathways; (ii) the ACSL4/GPX4 ratio, a metabolic indicator of ferroptosis sensitivity, was independently validated (ρ = −0.157, p = 5.6×10⁻³); and (iii) Monocle3 pseudotime analysis temporally anchors ferroptosis driver upregulation to the NDUFB7 breakpoint (τ=0.42, p=4.3×10⁻⁴).

**Reviewer defense**: "While NDUFB7 loss is associated with a broad death signature, we emphasize ferroptosis because (1) it is the only pathway with both a validated metabolic ratio index (ACSL4/GPX4) and a druggable defense system (GPX4, SLC7A11); (2) the temporal synchronization in pseudotime analysis specifically anchors ferroptosis driver upregulation to the NDUFB7 breakpoint; and (3) ferroptosis is the form of regulated cell death most tightly linked to mitochondrial dysfunction and iron metabolism, making it the mechanistically most coherent interpretation of an OXPHOS subunit defect."
