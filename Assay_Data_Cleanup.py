import pandas as pd
import numpy as np

# Data Extraction
assay_summary = pd.read_excel("Tox21CellBasedAssays.xlsx", sheet_name="Assays summary info")
compounds_code = pd.read_excel("Tox21CellBasedAssays.xlsx", sheet_name="Compounds code")
assay_endpoint_results = pd.read_excel("Tox21CellBasedAssays.xlsx", sheet_name="Assay endpoint results")
assay_full = pd.read_excel("Tox21CellBasedAssays.xlsx", sheet_name="Full assay info")
extra_assays = pd.read_excel("Tox21CellBasedAssays.xlsx", sheet_name="Extra assays")

drug_nasdili = pd.read_csv("medNAS_dili_tt.csv", encoding = "ISO-8859-1")

# Data Cleanup

# Set 2 - cpd2: Rosiglitazone Maleate (LessDILI)
# Set 2 - cpd3: Troglitazone (MoreDILI)

# Drop rows with both NAN values in endpoint results
rosig_trog_ac50_results = assay_endpoint_results[["aeid", "Set 2 - cpd2", "Set 2 - cpd3"]].copy()
rosig_trog_ac50_results = rosig_trog_ac50_results.dropna(thresh=2)

# Rename the column name
rosig_trog_ac50_results.rename(columns={"Set 2 - cpd2": "Rosiglitazone_Maleate",
                                   "Set 2 - cpd3": "Troglitazone"},
                                   inplace=True)

# Add assay endpoint component information
component = assay_endpoint_results[["aeid",
                                    "assay_component_endpoint_name"]].copy()


# Add assay target information
target = assay_full[["aeid",
                     "biological_process_target",
                     "intended_target_family", "intended_target_family_sub",
                     "intended_target_gene_symbol", "technological_target_gene_symbol",
                     "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                     "intended_target_gene_name", "technological_target_gene_name",
                     ]].copy()


rosig_trog_results = rosig_trog_ac50_results.merge(target, how="left", on="aeid").copy()
rosig_trog_results = rosig_trog_results.merge(component, how="left", on="aeid")


# Troglitazone

trog_results_filter = rosig_trog_results["Troglitazone"] < 1000000
trog_nasdili_filter = drug_nasdili["chnm"]=="troglitazone"

trog_results = rosig_trog_results[trog_results_filter].drop(columns=["Rosiglitazone_Maleate"]).copy()
trog_nasdili = drug_nasdili[trog_nasdili_filter].copy()

# Calculate NAS
trog_cmax = drug_nasdili.loc[drug_nasdili["chnm"]=="troglitazone", "CmaxStand"].values
# CmaxStand: 6.3867
trog_results["Trog_NAS"] = (6.3867-trog_results["Troglitazone"])/6.3867
trog_results = trog_results[["aeid", "Troglitazone", "Trog_NAS",
                             "biological_process_target",
                             "intended_target_family", "intended_target_family_sub",
                             "intended_target_gene_symbol", "technological_target_gene_symbol",
                             "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                             "intended_target_gene_name", "technological_target_gene_name"
                             ]]

# Rosiglitazone Maleate

rosig_results_filter = rosig_trog_results["Rosiglitazone_Maleate"] < 1000000
rosig_nasdili_filter = drug_nasdili["chnm"]=="rosiglitazone maleate"

rosig_results = rosig_trog_results[rosig_results_filter].drop(columns=["Troglitazone"]).copy()
rosig_maleate_nasdili = drug_nasdili[rosig_nasdili_filter].copy()
# No Rosiglitazone Maleate

# Calculate NAS
# rosig_cmax = drug_nasdili.loc[drug_nasdili["chnm"]=="rosiglitazone", "CmaxStand"].values
# CmaxStand: 1.26
rosig_results["Rosig_NAS"] = (1.26-rosig_results["Rosiglitazone_Maleate"])/1.26
rosig_results = rosig_results[["aeid", "Rosiglitazone_Maleate", "Rosig_NAS",
                               "biological_process_target",
                               "intended_target_family", "intended_target_family_sub",
                               "intended_target_gene_symbol", "technological_target_gene_symbol",
                               "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                               "intended_target_gene_name", "technological_target_gene_name"
                             ]]


# Graph 1: Assay Status and Summary
# rosig_results = rosig_trog_results.drop(columns=["Troglitazone"]).copy()
# rosig_results= rosig_results[rosig_results["Rosiglitazone_Maleate"].notnull()]

# trog_results = rosig_trog_results.drop(columns=["Rosiglitazone_Maleate"]).copy()
# trog_results = trog_results[trog_results["Troglitazone"].notnull()]

# Remove rows with NA values in EC50
rosig_trog_results_activated = rosig_trog_results[rosig_trog_results["Rosiglitazone_Maleate"].notnull()].copy()
rosig_trog_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Troglitazone"].notnull()]

# Troglitazone Activated only
trog_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Troglitazone"] < 1000000].copy()
trog_results_activated = trog_results_activated[trog_results_activated["Rosiglitazone_Maleate"] == 1000000]

# Rosiglitazone Maleate Activated only
rosig_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Rosiglitazone_Maleate"] < 1000000].copy()
rosig_results_activated = rosig_results_activated[rosig_results_activated["Troglitazone"] == 1000000]

# Both Activated
both_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Troglitazone"] < 1000000].copy()
both_results_activated = both_results_activated[both_results_activated["Rosiglitazone_Maleate"] < 1000000]

rosig_trog_results.to_csv("rosig_trog_results.csv")
rosig_results.to_csv("rosig_results.csv")
trog_results.to_csv("trog_results.csv")
rosig_results_activated.to_csv("rosig_results_activated.csv")
trog_results_activated.to_csv("trog_results_activated.csv")
both_results_activated.to_csv("both_results_activated.csv")