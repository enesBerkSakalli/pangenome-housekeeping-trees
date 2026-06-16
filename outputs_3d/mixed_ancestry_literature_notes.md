# Mixed-Ancestry Overlay Notes

This visualization uses the tree as a backbone and adds reticulation chords for
variants whose secondary population allele frequency is substantial relative to
the dominant population allele frequency.

Interpretation:

- The tree topology is a constrained population-expanded UPGMA tree.
- The reticulation chords are AF-found mixed-ancestry candidates, not proven
  recombination events.
- A true recombination/ancestral-history reconstruction would require phased
  haplotypes, local ancestry calls, LD blocks, or ARG inference.

Literature basis:

- Human pangenome references are intended to represent diverse haplotypes and
  variation more completely than a single linear reference:
  https://www.nature.com/articles/s41586-023-05896-x
- Human population panels such as 1000 Genomes provide global variant-frequency
  context across populations:
  https://www.nature.com/articles/nature15393
- gnomAD stratifies allele frequencies by genetic ancestry group and has
  explicitly discussed local ancestry for admixed groups:
  https://gnomad.broadinstitute.org/news/2023-11-genetic-ancestry/
- Local ancestry inference can refine allele-frequency interpretation in
  admixed groups:
  https://www.nature.com/articles/s41467-025-63340-2
- Recombination means different genomic segments can have different genealogical
  histories, motivating ARG or reticulation-style views rather than a single
  strict tree:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10796009/

Current implementation:

- `MIX_RATIO_THRESHOLD = 0.35`
- `MIX_MIN_SECONDARY_AF = 0.001`
- Chords aggregate candidates by `dominant_population -> secondary_population`
  within each gene layer.
- Diamond markers identify individual variant leaves that drive the aggregate
  chords.
