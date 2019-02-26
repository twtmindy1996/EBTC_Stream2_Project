import pandas as pd
rosig_results = pd.read_csv("rosig_results.csv")
trog_results = pd.read_csv("trog_results.csv")
rosig_trog_results = pd.read_csv("rosig_trog_results.csv")
rosig_results["biological_process_target"].unique()
trog_results["biological_process_target"].unique()