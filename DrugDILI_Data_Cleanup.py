import pandas as pd
drug_nasdili = pd.read_csv("medNAS_dili_tt.csv", encoding = "ISO-8859-1")

trog_cmax = drug_nasdili.loc[drug_nasdili["chnm"]=="troglitazone", "CmaxStand"].values

# Drop useless columns
drug_nasdili = drug_nasdili.drop(columns=["casn", "code"])

# Rename column
drug_nasdili = drug_nasdili.rename(columns={"tgt_abbr": "test_target",
                             "testtarget": "test_target_name",
                             "testname": "assay_component_endpoint_name",
                             "CmaxStand": "Cmax",
                             "EC": "AC50"})
drug_nasdili = drug_nasdili[["chnm", "Classification", "Cmax", "AC50", "NAS",
                             "test_target", "test_target_name", "assay_component_endpoint_name"]]

# Troglitazone
trog_nasdili_filter = drug_nasdili["chnm"] == "troglitazone"

trog_nasdili = drug_nasdili[trog_nasdili_filter].copy()

drug_nasdili_test = drug_nasdili.drop(drug_nasdili[trog_nasdili_filter].index).copy()

# Rosiglitazone Maleate
rosig_nasdili_filter = drug_nasdili["chnm"]=="rosiglitazone maleate"
rosig_maleate_nasdili = drug_nasdili[rosig_nasdili_filter].copy()
# No Rosiglitazone Maleate

drug_nasdili.to_csv("drug_nasdili.csv")
trog_nasdili.to_csv("trog_nasdili.csv")
drug_nasdili_test.to_csv("drug_nasdili_test.csv")

