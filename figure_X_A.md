# Figure X_A Generation Notes

## Current State
The figure_X_A.R script now includes a comprehensive analysis of INS expression, XBP1 splicing, and ER stress gene expression in beta cells, comparing protected and unprotected donors. The script generates a multi-panel figure and performs statistical tests.

## Objective
Create a comprehensive figure for the Cell Reports article that illustrates the relationship between INS expression, XBP1 splicing, and ER stress in the context of protected status (true/false).

## Key Files and Functions

### 1. AAA_defaults_and_constants.R
- Contains important constants, color palettes, and gene lists for the analysis.
- Key elements used:
  - `diabetes_palette` for consistent color coding
  - `custom_theme` for consistent figure styling

### 2. starsolo_xbp1_psi_per_cell.R
- Contains the function to calculate XBP1 PSI (Percent Spliced In) per cell

## Data Preparation
- Filters metadata for Beta cells and NODM status
- Calculates mature, total, and nascent INS counts
- Computes log-transformed INS expression values
- Calculates XBP1 PSI using the starsolo_xbp1_psi_per_cell function

## Figure Generation Strategy

1. **Violin Plots**
   - Mature INS Expression
   - XBP1 Splicing (PSI)

2. **ER Stress Heatmap**
   - Includes key ER stress genes: HSPA5, DDIT3, ATF4, XBP1, ERN1, EIF2AK3, ATF6
   - Compares expression between protected and unprotected groups

3. **XBP1 PSI vs ER Stress Scatter Plot**
   - X-axis: XBP1 PSI
   - Y-axis: Core UPR Stress Score
   - Color-coded by protected status
   - Includes linear regression lines

4. **Statistical Analysis**
   - Wilcoxon rank-sum tests for INS expression, XBP1 PSI, and ER stress scores
   - Correlation analysis between XBP1 PSI and ER stress score

5. **Metadata Summary**
   - Number of donors, mean age, and sex distribution for protected and unprotected groups

## Styling and Polish
- Consistent use of color palette (diabetes_palette)
- Custom theme applied to all plots
- Publication-quality figure layout and annotations

## Next Steps
1. Review the generated figure and statistical results
2. Make any necessary adjustments to improve clarity or address specific research questions
3. Prepare figure legend and methods description for the manuscript

## Important Files to Reference
- /scratch/domeally/DCD.tienhoven_scRNAseq.2024/figure_X_A.R (main script for figure generation)
- /scratch/domeally/DCD.tienhoven_scRNAseq.2024/figure_X_A.md (this file)
- /scratch/domeally/DCD.tienhoven_scRNAseq.2024/R/AAA_defaults_and_constants.R (for color palettes and themes)
- /scratch/domeally/DCD.tienhoven_scRNAseq.2024/R/starsolo_xbp1_psi_per_cell.R (for XBP1 PSI calculation)