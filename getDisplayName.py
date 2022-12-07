#! /usr/bin/env python3
# Translates full pango lineages into disGenerates color coding


# Translation table for pangolin and WHO names of variants
pangolin2WHO = {'B.1.1.7': 'Alpha', 'B.1.351': 'Beta', 'P.1': 'Gamma', 'B.1.427': 'Epsilon', 'B.1.429': 'Epsilon',
                'B.1.525': 'Eta', 'B.1.526': 'Iota', 'B.1.617.1': 'Kappa', 'B.1.621': 'Mu', 'B.1.621.1': 'Mu',
                'P.2': 'Zeta', 'B.1.617.3': 'B.1.617.3', 'B.1.617.2': 'Delta', 'AY': 'Delta',
                'B.1.1.529': 'Omicron', 'BA.1': 'BA.1', 'BA.2': 'BA.2', 'BA.2.12': 'BA.2.12',
                'BA.2.75':'BA.2.75', 'BA.4.6':'BA.4.6',
                'BF.7':'BF.7', 'BQ.1':'BQ.1', 'BQ.1.1':'BQ.1.1',
                'BA.3': 'BA.3', 'BA.4':'BA.4', 'BA.5':'BA.5', 'wt': 'wt', 'wt-wuhan': 'wt',
                'BC':'BA.1', 'BG': 'BA.2.12', 'BE':'BA.5',
                'A.21': 'Bat', 'other': 'Other', 'A': 'wt', 'Error':'Error'}



# Convert each variant to a WHO-compatible name, if one exists
def getDisplayName(pangolinName):
    if pangolinName in pangolin2WHO.keys():
        # Exact correspondance to a published name
        return pangolin2WHO[pangolinName]
    elif pangolinName in pangolin2WHO.values():
        # Already an exact match to a published name
        return pangolinName
    else:
        # Check if this is a sublineage of a defined lineage
        for i in range(pangolinName.count('.'), 0, -1):
            superVariant = '.'.join(pangolinName.split('.')[0:i])
            if superVariant in pangolin2WHO:
                return pangolin2WHO[superVariant]

        # No name seems to correspond to it, return itself
        return pangolinName



from matplotlib.colors import to_hex
import matplotlib.pyplot as plt

# Assign a pre-determined color to each display name
rgbColors = plt.get_cmap('tab20b').colors + plt.get_cmap('tab20c').colors[0:16]
colorCycle = [to_hex(color) for color in rgbColors]
def getColor (var_name):
    if var_name.lower() == 'other':
        return '#BBBBBB'
    else:
        color_idx = hash(var_name)%len(colorCycle)
        return colorCycle[color_idx]



# Process the original abundance estimates by Freyja
import pandas as pd
import numpy as np

def import_freyja_demix(filename):
    freyja_raw = pd.read_table(filename, index_col=0)

    # Parse the rows with the detailed subvariant breakdown
    lineages = np.array(freyja_raw.loc['lineages'][0].split(' '))
    abundances = np.array([float(x) for x in freyja_raw.loc['abundances'][0].split(' ')])
    
    # Freyja sometimes throws nan or inf for abundance. Eliminate such entries, if exist.
    valid = np.isfinite(abundances)
    lineages = lineages[valid]
    abundances = abundances[valid]
    
    freyja_names = [ getDisplayName(x) for x in lineages ]
    return (lineages, abundances, freyja_names)



if __name__ == '__main__':
    pass
    
    
