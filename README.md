cpra
=============================================================================

Calculated Panel Reactive Antibody with R

This is an R package with functions to compute calculated panel reactive antibody (cPRA) values from given allelic and haplotypic HLA frequencies.

In here, we have three main functions.: cpras() which computes cPRAs values for a group a patients with their HLA sensitization; txComps() which computes a transplantability score assuming AB0 compatibility for patients' HLA sensitization; and txIdents() which computes a transplantability score assuming AB0 identical for patients' HLA sensitization.

Available HLA allelic and haplotipic frequencies in this package are from Portuguese bone marrow volunteer donors published previously (Bruno A. Lima, Helena Alves. HLA-A, -C, -B, and -DRB1 allelic and haplotypic diversity in bone marrow volunteer donors from Northern Portugal. Organs, Tissues & Cells. Volume 16(1), 2013, March: 19-26). This frequencies can be changed in input frequencies allowing for cPRA calculation population appropriated.
