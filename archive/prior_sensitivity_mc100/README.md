# Archived fits of the prior sensitivity analysis over 100 datasets

This folder holds the 1,300 fits of the prior sensitivity analysis over 100
datasets (13 prior settings times 100 datasets) as RDS summary checkpoints in
`checkpoints/`, and the 39 tables of the original runs in `tables/`. The
posterior draws are not included. Script 02b rebuilds all tables and figures
from the saved summaries and does not need them.

`sensitivity_mc100/reference/ARCHIVE_MANIFEST.csv` lists the MD5 checksum of
every file, and script 02b checks all files against it before reading them.
Script 02b writes its results to `results/prior_sensitivity/mc100/` and never
modifies this folder.
