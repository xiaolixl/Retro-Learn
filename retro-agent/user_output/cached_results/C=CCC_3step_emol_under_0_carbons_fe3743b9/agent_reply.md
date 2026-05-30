## 1. Target Information

- **Name:** 1-Butene  
- **SMILES:** `C=CCC`  
- **Molecular weight:** 56.06 g/mol  
- **In stock:** No  

## 2. Result Summary

**Mode:** Multi‑step retrosynthesis (3 steps)  
**Best route score:** 0.6  

The recommended forward synthesis proceeds from ethyl bromide and ethylmagnesium chloride → butane → 2‑bromobutane → 1‑butene.

### Forward Synthesis Steps

| Step | Reactants | Conditions | Product | Step Score |
|------|-----------|------------|---------|------------|
| 1 | `CCBr` (ethyl bromide) + `CC[Mg]Cl` (ethylmagnesium chloride) | Grignard coupling (no specific conditions given) | `CCCC` (butane) | 0.4 |
| 2 | `CCCC` (butane) + `BrBr` (bromine) | Br₂ (radical bromination, e.g., light/heat) | `CCC(C)Br` (2‑bromobutane) | 0.2 |
| 3 | `CCC(C)Br` (2‑bromobutane) | (CH₃)₃COK (potassium tert‑butoxide) in (CH₃)₃COH (tert‑butanol) | `C=CCC` (1‑butene) | 0.0 |

**Rationale for the sequence:**  
- Step 1 couples two C2 fragments to form butane.  
- Step 2 introduces a bromine at the 2‑position (radical bromination, non‑selective but workable).  
- Step 3 eliminates HBr using a bulky base to favour the less substituted alkene (Hofmann elimination) → 1‑butene.

## 3. Key Observations

- **Preferred reactants:** None were specified by the user; no match was required.  
- **Stock status:** All three leaf reactants (`BrBr`, `CCBr`, `CC[Mg]Cl`) are marked as **not in stock**. They would need to be synthesised or purchased separately.  
- **Viability:** The route is computationally suggested. The low step scores (especially step 3 with 0.0) indicate the disconnections have low predicted plausibility; in practice, the bromination step would give a mixture of isomers and the elimination step would require careful control to obtain pure 1‑butene.  

⚠️ **The routes above are computational suggestions, not experimentally validated. Expert review is required before experimental use.**