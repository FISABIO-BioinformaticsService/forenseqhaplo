# Population frequency files for Familias

This directory is intended for population frequency files distributed with
`forenseqhaplo` for use in Familias.

The frequency files must contain population-level reference information only.
They must not contain individual-level genotypes, sample identifiers, case
information, or other confidential forensic data.

## Location after package installation

Files placed in this source directory are installed under:

```r
system.file(
  "population_frequencies",
  package = "forenseqhaplo",
  mustWork = TRUE
)
```

List all distributed files with:

```r
frequency_dir <- system.file(
  "population_frequencies",
  package = "forenseqhaplo",
  mustWork = TRUE
)

list.files(frequency_dir, full.names = TRUE)
```

Copy a selected file to the current working directory with:

```r
file.copy(
  from = file.path(frequency_dir, "FILENAME.txt"),
  to = "FILENAME.txt",
  overwrite = FALSE
)
```

## Required provenance

Before a frequency file is distributed, document at least:

- Filename and version.
- Population name and geographic scope.
- Marker panel and haplotype definition.
- Sample size and frequency-estimation method.
- Source publication, report, or dataset.
- DOI or stable URL when available.
- Licence or explicit redistribution permission.
- Date of preparation and responsible author.
- Familias version or import format used for validation.

A manifest describing all distributed files should be added when the final
frequency tables are available.
