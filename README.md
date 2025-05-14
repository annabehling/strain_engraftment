# strain_engraftment

## The Gut Bugs Trial

The Gut Bugs Trial was a double-blinded randomised placebo-controlled trial that assessed the efficacy of fecal microbiome transplantation (FMT) to treat adolescent (aged 14-18 years) obesity and improve metabolism. In total, 28 acid-resistant capsules containing gut microbiota from 4 sex-matched donors were administered to the FMT recipients over two consecutive days. Recipients were clinically assessed at baseline, and at 6-, 12-, and 26-weeks post-treatment. Donor and recipient stool samples collected at each clinical assessment underwent shotgun metagenomic sequencing. In total, 381 metagenomic sequencing files were analysed.

Protocol paper: https://doi.org/10.1136/bmjopen-2018-026174

Trial paper: https://doi.org/10.1001/jamanetworkopen.2020.30415

Metagenomic data: https://doi.org/10.1186/s40168-021-01060-7

## Donor strain engraftment analysis

This repository contains the following R scripts used to analyse metagenomic data from the Gut Bugs Trial for evidence of donor strain engraftment.

- 0_gb_data_formatting.R : format data for analysis
- 1_gb_optimal_strain_matching_threshold_selection.R : select optimal DNA distance threshold for strain matching
- 2_gb_unblind_strain_engraftment_analysis.R : perform strain engraftment analysis using optimal threshold
- 3_gb_unblind_donor_engraftment_efficiency_factors.R : identify factors influencing engraftment efficiency
- 4_gb_assembly_stats.R : generate assembly statistics
