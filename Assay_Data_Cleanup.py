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

# Extract "zebrafish" data
# zebrafish_results_filter = assay_endpoint_results["organism"] == "zebrafish"
# zebrafish_results = assay_endpoint_results[zebrafish_results_filter].copy()


# Drop "zebrafish" data
# assay_endpoint_results = assay_endpoint_results.drop(assay_endpoint_results[zebrafish_results_filter].index)

# Set 2 - cpd2: Rosiglitazone Maleate (LessDILI)
# Set 2 - cpd3: Troglitazone (MoreDILI)

# Drop rows with both NAN values in endpoint results
rosig_trog_ac50_results = assay_endpoint_results[["aeid", "Set 2 - cpd2", "Set 2 - cpd3"]].copy()
rosig_trog_ac50_results = rosig_trog_ac50_results.dropna(thresh=2)

# Rename the column name
rosig_trog_ac50_results.rename(columns={"Set 2 - cpd2": "Rosig_AC50",
                                        "Set 2 - cpd3": "Trog_AC50"},
                               inplace=True)

# Add assay endpoint component information
component = assay_endpoint_results[["aeid",
                                    "assay_component_endpoint_name",
                                    "assay_name"]].copy()


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

# Test Target: Gene symbol (if not null) / Assay component endpoint name (if Gene null)
rosig_trog_results["test_target"] = rosig_trog_results["technological_target_gene_symbol"]
rosig_trog_results["test_target"] = rosig_trog_results["test_target"].fillna(rosig_trog_results["intended_target_gene_symbol"])
rosig_trog_results["biological_process_target"] = rosig_trog_results["biological_process_target"].fillna("zebrafish whole embryo")
rosig_trog_results["test_target"] = rosig_trog_results["test_target"].fillna(rosig_trog_results["biological_process_target"])

# Test Status
rosig_trog_results["Rosig_Status"] = rosig_trog_results["Rosig_AC50"]
rosig_trog_results.loc[rosig_trog_results["Rosig_Status"] < 1000000, "Rosig_Status"] = "R"
rosig_trog_results["Rosig_Status"] = rosig_trog_results["Rosig_Status"].fillna("N")
rosig_trog_results["Rosig_Status"] = rosig_trog_results["Rosig_Status"].replace(1000000, "nR")

rosig_trog_results["Trog_Status"] = rosig_trog_results["Trog_AC50"]
rosig_trog_results.loc[rosig_trog_results["Trog_Status"] < 1000000, "Trog_Status"] = "T"
rosig_trog_results["Trog_Status"] = rosig_trog_results["Trog_Status"].fillna("N")
rosig_trog_results["Trog_Status"] = rosig_trog_results["Trog_Status"].replace(1000000, "nT")

rosig_trog_results["test_status"] = rosig_trog_results["Rosig_Status"].str.cat(rosig_trog_results["Trog_Status"],                                                                               sep="-")
rosig_trog_results = rosig_trog_results.drop(columns=["Rosig_Status", "Trog_Status"])


# Re-arrange the column position
rosig_trog_results = rosig_trog_results[["aeid", "Rosig_AC50", "Trog_AC50",
                                         "test_status",
                                         "biological_process_target",
                                         "intended_target_family", "intended_target_family_sub",
                                         "test_target",
                                         "intended_target_gene_symbol", "technological_target_gene_symbol",
                                         "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                                         "intended_target_gene_name", "technological_target_gene_name",
                                         "assay_component_endpoint_name",
                                         "assay_name"
                                         ]]
# Build bins for categorizing Activation Level
activation_level = ["Low", "Modest", "High"]
bins = [-10000, -4, 0, 10000]
# (-10000, -4), [-4, 0), [0, 10000)

# Troglitazone

# Activated
trog_results_filter = rosig_trog_results["Trog_AC50"] < 1000000
trog_results = rosig_trog_results[trog_results_filter].drop(columns=["Rosig_AC50"]).copy()

# Calculate NAS
# CmaxStand: 6.3867
trog_results["Trog_NAS"] = (6.3867-trog_results["Trog_AC50"])/6.3867

# Divide Activation Level based on NAS
trog_results["Trog_Activation_Level"] = pd.cut(trog_results["Trog_NAS"], bins=bins,
                                               include_lowest=True,
                                               labels=activation_level, right=False)

trog_results = trog_results[["aeid", "Trog_AC50", "Trog_NAS", "Trog_Activation_Level",
                             "test_status",
                             "biological_process_target",
                             "intended_target_family", "intended_target_family_sub",
                             "test_target",
                             "intended_target_gene_symbol", "technological_target_gene_symbol",
                             "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                             "intended_target_gene_name", "technological_target_gene_name",
                             "assay_component_endpoint_name",
                             "assay_name"
                             ]]


# Rosiglitazone Maleate

# Activated
rosig_results_filter = rosig_trog_results["Rosig_AC50"] < 1000000
rosig_results = rosig_trog_results[rosig_results_filter].drop(columns=["Trog_AC50"]).copy()

# Calculate NAS
# rosig_cmax = drug_nasdili.loc[drug_nasdili["chnm"]=="rosiglitazone", "CmaxStand"].values
# CmaxStand: 1.26
rosig_results["Rosig_NAS"] = (1.26-rosig_results["Rosig_AC50"])/1.26

# Divide Activation Level based on NAS
rosig_results["Rosig_Activation_Level"] = pd.cut(rosig_results["Rosig_NAS"], bins=bins,
                                                 include_lowest=True,
                                                 labels=activation_level, right=False)

rosig_results = rosig_results[["aeid", "Rosig_AC50", "Rosig_NAS","Rosig_Activation_Level",
                               "test_status",
                               "biological_process_target",
                               "intended_target_family", "intended_target_family_sub",
                               "test_target",
                               "intended_target_gene_symbol", "technological_target_gene_symbol",
                               "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                               "intended_target_gene_name", "technological_target_gene_name",
                               "assay_component_endpoint_name",
                               "assay_name"
                               ]]


# Activated Situation Sub-datasets

# Remove rows with NA values in EC50
rosig_trog_results_activated = rosig_trog_results[rosig_trog_results["Rosig_AC50"].notnull()].copy()
rosig_trog_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Trog_AC50"].notnull()]


# Troglitazone Activated only
trog_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Trog_AC50"] < 1000000].copy()
trog_results_activated = trog_results_activated[trog_results_activated["Rosig_AC50"] == 1000000]

trog_results_activated["Trog_NAS"] = (6.3867-trog_results_activated["Trog_AC50"])/6.3867
trog_results_activated["Trog_Activation_Level"] = pd.cut(trog_results_activated["Trog_NAS"], bins=bins,
                                                         include_lowest=True,
                                                         labels=activation_level, right=False)
trog_results_activated = trog_results_activated[["aeid", "Rosig_AC50", "Trog_AC50",
                                                 "Trog_NAS","Trog_Activation_Level",
                                                 "test_status",
                                                 "biological_process_target",
                                                 "intended_target_family", "intended_target_family_sub",
                                                 "test_target",
                                                 "intended_target_gene_symbol", "technological_target_gene_symbol",
                                                 "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                                                 "intended_target_gene_name", "technological_target_gene_name",
                                                 "assay_component_endpoint_name",
                                                 "assay_name"
                                                 ]]

# Rosiglitazone Maleate Activated only
rosig_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Rosig_AC50"] < 1000000].copy()
rosig_results_activated = rosig_results_activated[rosig_results_activated["Trog_AC50"] == 1000000]
rosig_results_activated["Rosig_NAS"] = (1.26-rosig_results_activated["Rosig_AC50"])/1.26
rosig_results_activated["Rosig_Activation_Level"] = pd.cut(rosig_results_activated["Rosig_NAS"], bins=bins,
                                                           include_lowest=True,
                                                           labels=activation_level, right=False)
rosig_results_activated = rosig_results_activated[["aeid", "Rosig_AC50", "Trog_AC50",
                                                   "Rosig_NAS","Rosig_Activation_Level",
                                                   "test_status",
                                                   "biological_process_target",
                                                   "intended_target_family", "intended_target_family_sub",
                                                   "test_target",
                                                   "intended_target_gene_symbol", "technological_target_gene_symbol",
                                                   "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                                                   "intended_target_gene_name", "technological_target_gene_name",
                                                   "assay_component_endpoint_name",
                                                   "assay_name"
                                                   ]]

# Both Activated
both_results_activated = rosig_trog_results_activated[rosig_trog_results_activated["Trog_AC50"] < 1000000].copy()
both_results_activated = both_results_activated[both_results_activated["Rosig_AC50"] < 1000000]
both_results_activated["Rosig_NAS"] = (1.26-both_results_activated["Rosig_AC50"])/1.26
both_results_activated["Trog_NAS"] = (6.3867-both_results_activated["Trog_AC50"])/6.3867
both_results_activated["Rosig_Activation_Level"] = pd.cut(both_results_activated["Rosig_NAS"], bins=bins,
                                                          include_lowest=True,
                                                          labels=activation_level, right=False)
both_results_activated["Trog_Activation_Level"] = pd.cut(both_results_activated["Trog_NAS"], bins=bins,
                                                         include_lowest=True,
                                                         labels=activation_level, right=False)
both_results_activated = both_results_activated[["aeid", "Rosig_AC50", "Trog_AC50",
                                                 "Rosig_NAS","Rosig_Activation_Level",
                                                 "Trog_NAS", "Trog_Activation_Level",
                                                 "test_status",
                                                 "biological_process_target",
                                                 "intended_target_family", "intended_target_family_sub",
                                                 "test_target",
                                                 "intended_target_gene_symbol", "technological_target_gene_symbol",
                                                 "intended_target_entrez_gene_id", "technological_target_entrez_gene_id",
                                                 "intended_target_gene_name", "technological_target_gene_name",
                                                 "assay_component_endpoint_name",
                                                 "assay_name"
                                                 ]]

rosig_trog_results.to_csv("rosig_trog_results.csv")
rosig_results.to_csv("rosig_results.csv")
trog_results.to_csv("trog_results.csv")
rosig_results_activated.to_csv("rosig_results_activated.csv")
trog_results_activated.to_csv("trog_results_activated.csv")
both_results_activated.to_csv("both_results_activated.csv")