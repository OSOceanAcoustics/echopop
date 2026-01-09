import re
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def extract_sex_from_sheet(df: pd.DataFrame):
    """
    Tries to extract sex from any cell in the first 5 rows. Returns 'male', 'female', 'all', or 
    None.
    """
    text = " ".join(
        str(cell) for row in df.iloc[:5].values for cell in row if isinstance(cell, str)
    )
    match = re.search(r'\b(male|female|all)\b', text, re.IGNORECASE)
    if match:
        return match.group(1).lower()
    return None

def read_all_aged_length_haul_counts(filepath: Path):
    """
    Reads all sheets from the aged length haul counts report Excel file. This dynamically assigns 
    sex to the sheet name and returns a dictionary {sex: DataFrame}. Drops last column and row 
    (subtotal) for each column and row, respectively.
    """
    # Inspect the Excel file
    xls = pd.ExcelFile(filepath)

    # Pre-allocate sheet-specific dictionary
    sheet_dfs = {}
    
    # Iterate through each
    for sheet_name in xls.sheet_names:
        # ---- Read in the file
        df_raw = pd.read_excel(xls, sheet_name=sheet_name, header=None)
        # ---- Extract the sex-assigned sheet
        sex = extract_sex_from_sheet(df_raw)
        # ---- Row 1 (index 1): haul numbers - extract from column 2 to second-to-last 
        # ---- (skip subtotal col)
        haul_numbers = df_raw.iloc[1, 1:-1].dropna().values
        # ---- Row 2 onwards (index 2+): actual data starts here
        lengths = df_raw.iloc[2:, 0].values
        # ---- Find where "Subtotal" row is
        subtotal_idx = None
        for idx, val in enumerate(lengths):
            if isinstance(val, str) and "Subtotal" in str(val):
                subtotal_idx = idx + 2  # +2 because we started at row 2
                break
        # ---- Extract only the data rows (before Subtotal)
        if subtotal_idx:
            df_clean = df_raw.iloc[2:subtotal_idx, 1:-1].copy()
            df_clean.index = df_raw.iloc[2:subtotal_idx, 0].values
        else:
            df_clean = df_raw.iloc[2:, 1:-1].copy()
            df_clean.index = df_raw.iloc[2:, 0].values
        # ---- Set column names to haul numbers
        df_clean.columns = haul_numbers
        # ---- Convert to numeric
        df_clean = df_clean.apply(pd.to_numeric, errors='coerce')
        sheet_dfs[sex] = df_clean
    
    # Return the dictionary of sheets
    return 

def plot_all_aged_length_haul_count_comparisons(echopro_dict, echopop_dict, title_prefix=""):
    """
    For each sex in the intersection of both dicts, plot:
    - EchoPro heatmap
    - EchoPop heatmap
    - Difference heatmap
    """
    # Get all sex keys
    all_sexes = sorted(set(echopro_dict.keys()) | set(echopop_dict.keys()))
    n = len(all_sexes)
    
    # Set up plotting area
    fig, axes = plt.subplots(n, 3, figsize=(18, 8 * n), squeeze=False)
    
    # Iterate through for plots
    for i, sex in enumerate(all_sexes):
        # ---- Extract
        df1 = echopro_dict.get(sex)
        df2 = echopop_dict.get(sex)
        diff = df1 - df2
        if df1 is None or df2 is None:
            # ---- Blank axes if missing
            for j in range(3):
                axes[i, j].axis("off")
            continue  # skip if either is missing
        # ---- Ensure axes are sorted ascending
        df1 = df1.sort_index(axis=0).sort_index(axis=1)
        df2 = df2.sort_index(axis=0).sort_index(axis=1)
        diff = df1 - df2
        diff = diff.sort_index(axis=0).sort_index(axis=1)
        # ---- EchoPro
        sns.heatmap(df1, ax=axes[i, 0], cmap="viridis")
        axes[i, 0].set_title(f"{title_prefix}EchoPro ({sex})")
        axes[i, 0].set_xlabel("Haul number")
        # ---- Echopop
        sns.heatmap(df2, ax=axes[i, 1], cmap="viridis")
        axes[i, 1].set_title(f"{title_prefix}EchoPop ({sex})")
        axes[i, 1].set_ylabel("Length (cm)")        
        # ---- Net comparison
        sns.heatmap(diff, ax=axes[i, 2], cmap="coolwarm", center=0)
        axes[i, 2].set_title(f"Difference ({sex})")
    plt.tight_layout(h_pad=4.0)
    plt.show()