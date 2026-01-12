import re
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, Optional, Tuple
from matplotlib.ticker import MaxNLocator

def extract_sex_from_sheet(df: pd.DataFrame):
    """
    Tries to extract sex from any cell. Returns 'male', 'female', 'all', or None. Allows for 
    optional parentheses and whitespace, e.g., '(male)', ' (ALL) ', etc. Does not match joint 
    expressions like 'male+female', 'female/male', etc.
    """
    
    # Get text
    text = " ".join(
        str(cell) for row in df.values for cell in row if isinstance(cell, str)
    )
    
    # Find all matches that are surrounded by word boundaries or parentheses
    found = []
    for sex in ("all", "male", "female"):
        pattern = rf"(?<![+/])(?<![a-zA-Z0-9])\(?\s*{sex}\s*\)?(?![a-zA-Z0-9])(?!(\s*[+/]))"
        if re.search(pattern, text, re.IGNORECASE):
            found.append(sex)
    
    # Avoid ambiguous sex identifiers
    if len(found) == 1:
        return found[0]
    elif len(found) > 1:
        raise ValueError(f"Multiple sex identifiers found in sheet: {', '.join(found)}.")
    return None

def extract_quantity_type(df: pd.DataFrame):
    """
    Tries to extract the quantity-type (e.g., 'abundance', 'biomass', 'counts'). Returns 'counts', 
    'abundance', 'biomass', or None. Allows for optional parentheses and whitespace, e.g., 
    '(abundance)'.
    """ 
    
    # Get text
    text = " ".join(
        str(cell) for row in df.values for cell in row if isinstance(cell, str)
    )
    
    # Find all matches that are surrounded by word boundaries or parentheses
    found = []
    for quantity in ("counts", "abundance", "biomass"):
        pattern = rf"(?<![+/])(?<![a-zA-Z0-9])\(?\s*{quantity}\s*\)?(?![a-zA-Z0-9])(?!(\s*[+/]))"
        if re.search(pattern, text, re.IGNORECASE):
            found.append(quantity)
    # ---- Get unique in case of duplicates
    found = list(set(found))
            
    # Avoid ambiguous quantity identifiers
    if len(found) == 1:
        return found[0]
    elif len(found) > 1:
        raise ValueError(f"Multiple quantity identifiers found in sheet: {', '.join(found)}.")
    return None
    

def align_dataframes(df1: pd.DataFrame, df2: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align indices of both pandas.DataFrame objects.
    """
    
    # Determine which DataFrame has more unique indices/columns
    n1 = len(df1.index.unique())
    n2 = len(df2.index.unique())
    m1 = len(df1.columns.unique())
    m2 = len(df2.columns.unique())
    
    # Use the index/columns from the DataFrame with the most unique entries
    if (n1 * m1) >= (n2 * m2):
        ref_index = df1.index
        ref_columns = df1.columns
    else:
        ref_index = df2.index
        ref_columns = df2.columns
        
    # Get the union, preserving reference order and appending any extras from the other DataFrame
    all_index = ref_index.append(
        df2.index.difference(ref_index)
    ) if ref_index is df1.index else ref_index.append(df1.index.difference(ref_index))
    all_columns = ref_columns.append(
        df2.columns.difference(ref_columns)
    ) if ref_columns is df1.columns else ref_columns.append(df1.columns.difference(ref_columns))

    # Reindex both DataFrames to the full set, in the chosen order
    df1_aligned = df1.reindex(index=all_index, columns=all_columns).fillna(0.)
    df2_aligned = df2.reindex(index=all_index, columns=all_columns).fillna(0.)
    
    # Return the aligned DataFrames
    return df1_aligned, df2_aligned

def read_pivot_table_report(filepath: Path) -> Dict[str, pd.DataFrame]:
    """
    Reads all sheets from the aged length haul counts report Excel file. This dynamically assigns 
    sex to the sheet name and returns a dictionary {sex: DataFrame}. Drops last column and row 
    (subtotal) for each column and row, respectively.
    """
    
    # Validate file typing
    if not isinstance(filepath, Path):
        raise TypeError(
            f"Argument 'filepath' must be type pathlib.Path. Got type {type(filepath)}."
        )
        
    # Validate file existence
    if not filepath.exists():
        raise FileExistsError(
            f"The input file could not be located at {filepath.as_posix()}."
        )
    
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
        # ---- Row 1 (index 1): column indices - extract from column 2 to second-to-last 
        # ---- (skip subtotal col)
        column_numbers = df_raw.iloc[1, 1:-1].dropna().values
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
        # ---- Set column names to their respective name types
        column_name_mask = np.array(
            [isinstance(vt, float) | isinstance(vt, int) for vt in column_numbers]
        )
        # ---- Coerce the numerics to ints
        column_numbers[column_name_mask] = column_numbers[column_name_mask].astype(int)
        # ---- Capitalize otherwise
        column_numbers[~column_name_mask] = (
            np.char.capitalize(np.array(column_numbers[~column_name_mask], dtype=np.str_))
        )        
        # ---- Set columns
        df_clean.columns = column_numbers
        # ---- Convert to numeric
        df_clean = df_clean.apply(pd.to_numeric, errors="coerce")
        # ---- Get quantity type
        quantity_type = extract_quantity_type(df_raw)
        # ---- Add quantity type
        df_clean.columns = pd.MultiIndex.from_product([[quantity_type], df_clean.columns])
        # ---- Store
        sheet_dfs[sex] = df_clean
    
    # Return the dictionary of sheets
    return sheet_dfs

def plot_haul_count_comparisons(
    echopro: Dict[str, pd.DataFrame], 
    echopop: Dict[str, pd.DataFrame],
    save_filepath: Optional[Path] = None,
    show_plot: bool = True,
):
    """
    For each sex in the intersection of both dicts, plot:
    - EchoPro heatmap
    - EchoPop heatmap
    - Difference heatmap
    
    Parameters
    ----------
    echopro : Dict[str, pd.DataFrame]
    echopop : Dict[str, pd.DataFrame]
    save_filepath : Optional[Path]
        If provided, the plot will be saved to this location.
    show_plot : bool
        If True, the plot will be rendered in the user's session.
    """
    # Get all sex keys
    all_sexes = sorted(set(echopro.keys()) | set(echopop.keys()))
    n = len(all_sexes)
    
    # Set up plotting area
    fig, axes = plt.subplots(n, 3, figsize=(24, 6 * n), squeeze=False)
    
    # Iterate through for plots
    for i, sex in enumerate(all_sexes):
        # ---- Extract
        df1 = echopro.get(sex)["counts"]
        df2 = echopop.get(sex)["counts"]
        diff = df1 - df2
        if df1 is None or df2 is None:
            # ---- Blank axes if missing
            for j in range(3):
                axes[i, j].axis("off")
            continue  # skip if either is missing
        # ---- Ensure axes are sorted ascending
        df1 = df1.sort_index(axis=0).sort_index(axis=1)
        df2 = df2.sort_index(axis=0).sort_index(axis=1)
        diff = diff.sort_index(axis=0).sort_index(axis=1)
        # ---- Compute the shared ranges
        vmin = min(np.nanmin(df1.values), np.nanmin(df2.values))
        vmax = max(np.nanmax(df1.values), np.nanmax(df2.values))
        # ---- EchoPro
        hm1 = sns.heatmap(df1, ax=axes[i, 0], cmap="viridis", linewidths=0.5, linecolor="black",
                          vmin=vmin, vmax=vmax, cbar=True)
        hm1.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Echopop
        hm2 = sns.heatmap(df2, ax=axes[i, 1], cmap="viridis", linewidths=0.5, linecolor="black",
                          vmin=vmin, vmax=vmax, cbar=True)
        hm2.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Net comparison
        if np.all(diff.values == 0):
            hm3 = sns.heatmap(diff, ax=axes[i, 2], cmap="coolwarm", center=0, linewidths=0.5,
                              linecolor="black", cbar=True, vmin=-1, vmax=1)
        else:
            hm3 = sns.heatmap(diff, ax=axes[i, 2], cmap="coolwarm", center=0, linewidths=0.5,
                              linecolor="black", cbar=True)
        hm3.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Remove y-axis ticks from other plots in the row
        axes[i, 0].set_yticklabels(axes[i, 0].get_yticklabels(), rotation=0, ha='right')
        axes[i, 1].set_yticks([])
        axes[i, 2].set_yticks([])
        # ---- Only show x-axis ticks for bottom row
        if i == n - 1:
            for j in range(3):
                axes[i, j].tick_params(axis="x", labelrotation=45, labelsize=10)
        else:
            for j in range(3):
                axes[i, j].set_xticks([])
                axes[i, j].set_xticklabels([])
        # Add sex label as annotation to the left of the row
        axes[i, 0].annotate(
            sex.capitalize(),
            xy=(0, 0.5),
            xycoords="axes fraction",
            fontsize=16,
            ha="right",
            va="center",
            rotation=90,
            xytext=(-axes[i, 0].yaxis.labelpad - 40, 0),
            textcoords="offset points"
        )
                
    # Set column titles at the top
    col_titles = ["EchoPro counts", "Echopop counts", r"Differences ($\mathregular{\Delta}$counts)"]
    for j, title in enumerate(col_titles):
        axes[0, j].set_title(title, fontsize=18)
    
    # Set the outside/shared axis titles
    axes[1, 0].set_ylabel(r"Length ($\mathregular{\ell}$, cm)")
    axes[2, 1].set_xlabel("Haul number")

    plt.tight_layout(pad=3.0)
    plt.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.12)
    
    # Save?
    if save_filepath is not None:
        # ---- Validate directory
        if not save_filepath.parent.exists():
            save_filepath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_filepath, bbox_inches="tight", dpi=300)
        print(
            f"Figure saved at {save_filepath.as_posix()}."
        )
    
    # Show?
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
        
def plot_population_table_comparisons(
    echopro: Dict[str, pd.DataFrame], 
    echopop: Dict[str, pd.DataFrame],
    save_filepath: Optional[Path] = None,
    show_plot: bool = True,
    log_transform: bool = False
):
    """
    For each sex in the intersection of both dicts, plot:
    - EchoPro heatmap
    - EchoPop heatmap
    - Difference heatmap
    
    Parameters
    ----------
    echopro : Dict[str, pd.DataFrame]
    echopop : Dict[str, pd.DataFrame]
    save_filepath : Optional[Path]
        If provided, the plot will be saved to this location.
    show_plot : bool
        If True, the plot will be rendered in the user's session.
    log_transform: bool
        If True, al values wihtin the DataFrames are log-transformed (base-10).
    """
            
    # Get all sex keys
    all_sexes = sorted(set(echopro.keys()) | set(echopop.keys()))
    n = len(all_sexes)

    # Set up plotting area
    fig, axes = plt.subplots(n, 3, figsize=(24, 6 * n), squeeze=False)

    # Iterate through for plots
    for i, sex in enumerate(all_sexes):
        # ---- Extract
        df1 = echopro.get(sex)
        df2 = echopop.get(sex)
        # ---- Get the quantity-types
        df1_type = df1.columns.get_level_values(0).unique()[0]
        df2_type = df2.columns.get_level_values(0).unique()[0]
        # ---- Double-check
        if df1_type != df2_type:
            raise ValueError(
                f"Quantity types between the EchoPro ('{df1_type}' and Echopop ('{df2_type}').)"
            )
        if df1 is None or df2 is None:
            # ---- Blank axes if missing
            for j in range(3):
                axes[i, j].axis("off")
            continue  # skip if either is missing
        # ---- Align
        df1_aligned, df2_aligned = align_dataframes(df1[df1_type], df2[df2_type])
        # ---- Difference in population estimates
        diff = df1_aligned - df2_aligned
        # ---- Ensure axes are sorted ascending
        df1_aligned = df1_aligned.sort_index(axis=0)
        df2_aligned = df2_aligned.sort_index(axis=0)
        diff = diff.sort_index(axis=0)
        # ---- Log-transform?
        if log_transform:
            df1_aligned = np.log10(df1_aligned + 1)
            df2_aligned = np.log10(df2_aligned + 1)
            diff = df1_aligned - df2_aligned
        # ---- Compute the shared ranges
        vmin = min(np.nanmin(df1_aligned.values), np.nanmin(df2_aligned.values))
        vmax = max(np.nanmax(df1_aligned.values), np.nanmax(df2_aligned.values))
        # ---- EchoPro
        hm1 = sns.heatmap(df1_aligned, ax=axes[i, 0], cmap="viridis", linewidths=0.5, 
                          linecolor="black", vmin=vmin, vmax=vmax, cbar=True)
        hm1.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Echopop
        hm2 = sns.heatmap(df2_aligned, ax=axes[i, 1], cmap="viridis", linewidths=0.5, 
                          linecolor="black", vmin=vmin, vmax=vmax, cbar=True)
        hm2.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Net comparison
        if np.all(diff.values == 0):
            hm3 = sns.heatmap(diff, ax=axes[i, 2], cmap="coolwarm", center=0, linewidths=0.5,
                              linecolor="black", cbar=True, vmin=-1, vmax=1)
        else:
            hm3 = sns.heatmap(diff, ax=axes[i, 2], cmap="coolwarm", center=0, linewidths=0.5,
                              linecolor="black", cbar=True)
        hm3.collections[0].colorbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ---- Remove y-axis ticks from other plots in the row
        axes[i, 0].set_yticklabels(axes[i, 0].get_yticklabels(), rotation=0, ha='right')
        axes[i, 1].set_yticks([])
        axes[i, 2].set_yticks([])
        # ---- Only show x-axis ticks for bottom row
        if i == n - 1:
            for j in range(3):
                axes[i, j].tick_params(axis="x", labelrotation=45, labelsize=10)
        else:
            for j in range(3):
                axes[i, j].set_xticks([])
                axes[i, j].set_xticklabels([])
        # Add sex label as annotation to the left of the row
        axes[i, 0].annotate(
            sex.capitalize(),
            xy=(0, 0.5),
            xycoords="axes fraction",
            fontsize=16,
            ha="right",
            va="center",
            rotation=90,
            xytext=(-axes[i, 0].yaxis.labelpad - 40, 0),
            textcoords="offset points"
        )
                
    # Set column titles at the top
    if log_transform:
        col_titles = ["EchoPro abundance (log$_{10}$)", "Echopop abundance (log$_{10}$)", r"Differences ($\mathregular{\Delta}$abundance)"]
    else:
        col_titles = ["EchoPro abundance", "Echopop abundance", r"Differences ($\mathregular{\Delta}$abundance)"]
    for j, title in enumerate(col_titles):
        axes[0, j].set_title(title, fontsize=18)

    # Set the outside/shared axis titles
    axes[1, 0].set_ylabel(r"Length ($\mathregular{\ell}$, cm)")
    axes[2, 1].set_xlabel(r"Age ($\mathregular{\alpha}$, years)")

    plt.tight_layout(pad=3.0)
    plt.subplots_adjust(left=0.08, right=0.97, top=0.93, bottom=0.12)

    # Save?
    if save_filepath is not None:
        # ---- Validate directory
        if not save_filepath.parent.exists():
            save_filepath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_filepath, bbox_inches="tight", dpi=300)
        print(
            f"Figure saved at {save_filepath.as_posix()}."
        )

    # Show?
    if show_plot:
        plt.show()
    else:
        plt.close(fig)