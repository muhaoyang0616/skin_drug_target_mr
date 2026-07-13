# Validated R Environment

The independent validation used:

| Component | Version |
| --- | --- |
| R | 4.6.0 |
| TwoSampleMR | 0.7.6 |
| coloc | 5.2.3 |
| data.table | 1.18.4 |
| dplyr | 1.2.1 |
| yaml | 2.3.12 |
| readr | 2.2.0 |
| jsonlite | 2.0.0 |
| ggplot2 | 4.0.3 |
| ieugwasr | 1.1.0 |

The primary MR scripts require `data.table`, `yaml`, and `readr`.
Colocalization requires `coloc`, `jsonlite`, and `ggplot2`. The optional
OpenGWAS refresh requires `ieugwasr` plus an `OPENGWAS_JWT`. `TwoSampleMR` is
retained for consistency with the validated environment, although the final
single-variant estimates are explicit Wald ratios and were independently
reproduced in Python.
