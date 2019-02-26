import pandas as pd
drug_nasdili = pd.read_csv("medNAS_dili_tt.csv", encoding = "ISO-8859-1")


# Drop useless columns
drug_nasdili = drug_nasdili.drop(columns=["casn", "code"])

# Troglitazone
trog_nasdili_filter = drug_nasdili["chnm"] == "troglitazone"

trog_nasdili = drug_nasdili[trog_nasdili_filter].copy()
drug_nasdili_test = drug_nasdili.drop(drug_nasdili[trog_nasdili_filter].index).copy()

trog_nasdili.to_csv("trog_nasdili.csv")
drug_nasdili_test.to_csv("drug_nasdili_test.csv")
