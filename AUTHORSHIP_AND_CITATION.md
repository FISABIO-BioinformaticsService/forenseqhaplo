# Authorship and citation

This document records the confirmed contribution of each package author and
the fields that remain pending for the associated scientific paper.

## Confirmed package authors

| Person | Package role | Confirmed contribution |
|---|---|---|
| Vicente Soriano Chirona | `aut`, `cre` | Software development, package documentation, repository preparation, and maintenance |
| Sandra Carbó Ramírez | `aut` | Scientific research, algorithm development, filter definition, and threshold definition |
| Jorge Ruiz Ramírez | `aut` | Scientific research, algorithm development, filter definition, threshold definition, package testing, and validation |
| Alan Codoñer Alejos | `aut` | Scientific research, algorithm development, filter definition, and threshold definition |

## Package-role interpretation

- `aut`: substantial contribution to the package.
- `cre`: package maintainer and contact person.
- `ctb`: a smaller contribution that does not qualify for package authorship.
- `cph`: legal copyright holder. This must not be assigned without
  institutional confirmation.

The package author order currently places Vicente Soriano Chirona first because
he developed, documented, prepared, and maintains the software. Sandra Carbó
Ramírez, Jorge Ruiz Ramírez, and Alan Codoñer Alejos are package authors because
their scientific contributions define the algorithmic behaviour implemented by
the software.

## Associated scientific paper

The provisional title of the associated manuscript is:

```text
Nombre del Paper
```

Sandra Carbó Ramírez is the principal author.

The active provisional citation is stored in `inst/CITATION` as an
`Unpublished` manuscript and in `CITATION.cff` as a related reference. Both
records are explicitly marked as provisional and incomplete.

The final paper title, complete author list and order, journal, year, volume,
pages or article number, and DOI must be updated when they are definitive.

When the paper is accepted or published:

1. Replace the provisional `Unpublished` entry in `inst/CITATION` with a final
   `Article` entry.
2. Complete the reference under `references:` in `CITATION.cff`.
3. Add the final reference and DOI to `README.Rmd`, `README.md`, and the package
   website.
4. Ask users to cite both the software version and the scientific article.
5. Verify the CRediT contribution statement with all authors.

## Suggested CRediT contribution statements

These statements are a preparation aid and should be reviewed by the research
group before manuscript submission.

```text
Sandra Carbó Ramírez:
Conceptualization; Methodology; Investigation; Algorithm development;
Definition of filters and thresholds; Writing – original draft;
Writing – review & editing.

Jorge Ruiz Ramírez:
Methodology; Investigation; Algorithm development; Definition of filters and
thresholds; Validation; Software testing; Writing – review & editing.

Alan Codoñer Alejos:
Methodology; Investigation; Algorithm development; Definition of filters and
thresholds; Writing – review & editing.

Vicente Soriano Chirona:
Software – Lead; Documentation; Repository preparation and maintenance;
Writing – review & editing.
```

Only retain roles that accurately reflect each person's final contribution to
the manuscript.

## Files to update when the paper is final

- `CITATION.cff`
- `inst/CITATION`
- `README.Rmd`
- `README.md`
- Zenodo release metadata
- The manuscript's CRediT statement
