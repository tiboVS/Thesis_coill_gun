import os
import pandas as pd
import pickle
from scipy.interpolate import LinearNDInterpolator

# Pad naar de map met CSV-bestanden
map_pad = r"D:\MP\Nieuwe data\Hol projectiel\Data spoel n=600"
interpolator_pad = os.path.join(map_pad, "interpolator_n=600.pkl")

# Lees alle CSV-bestanden in de map
alle_data = []

for bestand in os.listdir(map_pad):
    if bestand.endswith(".csv"):
        df = pd.read_csv(os.path.join(map_pad, bestand))
        alle_data.append(df)

# Combineer alle data in één DataFrame
data = pd.concat(alle_data, ignore_index=True)

# **Filter: verwijder alle rijen met Positie > 120**
#data = data[data["Positie"] <= 120].copy()

# Sorteer op stroom en positie
data = data.sort_values(by=["Stroom", "Positie"]).reset_index(drop=True)

# Maak een interpolatiefunctie
interpolator = LinearNDInterpolator(data[["Stroom", "Positie"]], data["Kracht"])

# Sla de interpolator op
with open(interpolator_pad, "wb") as f:
    pickle.dump(interpolator, f)

print("Interpolator opgeslagen.")

# Laad de interpolator
with open(interpolator_pad, "rb") as f:
    interpolator = pickle.load(f)

# Voorbeeld: kracht bij I = 1.5 A en x = 22.75 mm
I_voorbeeld = 1.5
x_voorbeeld = 22.75
kracht_geschat = interpolator(I_voorbeeld, x_voorbeeld)

print(f"Geschatte kracht bij I={I_voorbeeld} A en x={x_voorbeeld} mm: {kracht_geschat}")
