import pandas as pd
import numpy as np

def parse_raw_fe(fe_path):
    """
    Parse an Agilent “FE” dump (.txt) to pull out probe metadata.
    Returns a DataFrame with at least these cols:
      – ProbeName
      – SystematicName   (chr:start–end)
      – ControlType      (0 = oligo, >0 = control)
    """
    cols = None
    rows = []
    with open(fe_path, 'r') as f:
        for line in f:
            # Look for the header line that starts with "FEATURES"
            if line.startswith('FEATURES'):
                cols = line.strip().split()[1:]
            # Once we have cols, capture every subsequent "DATA" row 
            elif cols is not None and line.startswith('DATA'):
                parts = line.strip().split()[1:]
                rows.append(parts)
    raw_df = pd.DataFrame(rows, columns=cols)
    # Convert numeric flags
    for c in ['FeatureNum','Row','Col','SubTypeMask','ControlType']:
        if c in raw_df:
            raw_df[c] = pd.to_numeric(raw_df[c], errors='coerce')
    return raw_df

# — adjust these paths! —
featxt = 'Project\Data\Raw\ChipSeq\Julian2016\FE\GSM1326860_E2f4_wt_replicate3_FE.txt'
proctxt = 'Project\Data\Raw\ChipSeq\Julian2016\Processed\GSM1326860_E2f4_wt_replicate3_Processed.txt'

# 1. Load raw probe metadata
annot = parse_raw_fe(featxt)

# 2. Load the stripped-down intensity table
proc = pd.read_csv(proctxt, sep='\t', 
                   usecols=['ID_REF','IP','INPUT'])

# 3. Merge them
merged = proc.merge(
    annot[['ProbeName','SystematicName','ControlType']],
    left_on='ID_REF', right_on='ProbeName',
    how='left'
)

# 4. Filter out control spots, keep only oligo probes
merged = merged[merged['ControlType']==0].copy()
# 4a. Drop probes with missing intensities (–1 or 0)
merged = merged[(merged['IP'] > 0) & (merged['INPUT'] > 0)].copy()
# 5. Compute log2 enrichment (M-value)
merged['M'] = np.log2(merged['IP'] / merged['INPUT'])

# Quick check
print( merged[['ID_REF','SystematicName','IP','INPUT','M']].head() )

# 6. Save the cleaned, annotated M-values
merged.to_csv(
    'Project\Data\Processed\Julian2016_rep3_Mvalues.csv',
    index=False
)
print(f"Wrote {len(merged)} probes with M-values to CSV.")
