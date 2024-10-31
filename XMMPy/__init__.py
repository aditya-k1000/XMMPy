import pandas as pd
import os

os.environ["SAS_PATH"] = os.environ["SAS_DIR"]

#Regular light curve plotting 
#Naming conventions
names = {"Broadband": "b",
         "Ultrasoft": "us", 
         "Soft": "s", 
         "Medium": "m", 
         "Hard": "h"
         }

#Energy band minimums
bands_min = {"b": 200.,
             "us": 200.,
             "s": 500.,
             "m": 2000.,
             "h": 4500.
             }

#Energy band maximums
bands_max = {"b": 12000.,
             "us": 500.,
             "s": 2000.,
             "m": 4500.,
             "h": 12000.
             }

#Plot colors
colors = {"b": "black",
          "us": "green",
          "s": "red",
          "m": "gold",
          "h": "blue"
             }

#Plot text
text = {"b": "Broadband (0.2 - 12 keV)",
        "us": "Ultrasoft (0.2 - 0.5 keV)",
        "s": "Soft (0.5 - 2 keV)",
        "m": "Medium (2 - 4.5 keV)",
        "h": "Hard (4.5 - 12 keV)"
             }
#Summary info file for plotting
values = pd.DataFrame({
    "Band": list(names.keys()),
    "Identifier": list(names.values()),
    "Energy Min": [bands_min[identifier] for identifier in names.values()],
    "Energy Max": [bands_max[identifier] for identifier in names.values()],
    "Color": [colors[identifier] for identifier in names.values()],
    "Text": [text[identifier] for identifier in names.values()]
}) 

#Hardness ratio plotting
#Hardness Ratios
hr = {"(M - S)/(M + S)": "m-s",
      "(H - S)/(H + S)": "h-s",
      "(S - (M + H)/(S + (M + H))": "s-m-h",
      "((H + M) - (S + U))/((H + M) + (S + U))": "h-m-s-u"
             }

#Hardness Ratio Colors
hr_colors = {"m-s": "gold",
             "h-s": "blue",
             "s-m-h": "green",
             "h-m-s-u": "purple"
             }

hr_values = pd.DataFrame({
    "Formula": list(hr.keys()),
    "Identifier": list(hr.values()),
    "Color": [hr_colors[identifier] for identifier in hr.values()]
})