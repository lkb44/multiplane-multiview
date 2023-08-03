import medviz as viz
from pathlib import Path
import pandas as pd

df = pd.DataFrame({
  "id": ["Patient-001_ax_ls.mha", "Patient-002_ax_ls.mha", "Patient-003_ax_ls.mha"],
})

path = "/Users/leobao/Documents/MultiPlanePipeline/Data/Axial_Image_NIFTY/Patient-001_ax_ls.nii.gz"

id_fn = lambda x: x.split("-")[-1][:3]

df["id"] = df["id"].apply(id_fn)
print(df)