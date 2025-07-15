from typing import List, Tuple
import numpy as np

def transect_mesh_region_2019(
    region : np.number
) -> Tuple[np.number, np.number, List[np.number], List[np.number]]:
    """
    Generate region-specific transect boundaries for the 2019 NWFSC survey.

    This function defines the spatial boundaries for three distinct survey regions used in the 
    2019 Northwest Fisheries Science Center (NWFSC) survey. Each region has specific transect 
    numbering schemes and boundary definitions.

    Parameters
    ----------
    region : np.number
        Region identifier (1, 2, or 3):
        - Region 1: Parallel transects to latitudes from south of SCB to west of Haida Gwaii
        - Region 2: Transects parallel to longitudes north of Haida Gwaii  
        - Region 3: Parallel transects to latitudes west of Haida Gwaii

    Returns
    -------
    Tuple[np.number, np.number, List[np.number], List[np.number]]
        A tuple containing:
        - transect_start (np.number): Starting transect number for the region
        - transect_end (np.number): Ending transect number for the region
        - transect_lower_bound (List[np.number]): Lower boundary values for each transect
        - transect_upper_bound (List[np.number]): Upper boundary values for each transect

    Examples
    --------
    >>> start, end, lower, upper = transect_mesh_region_2019(1)
    >>> print(f"Region 1: transects {start} to {end}")
    Region 1: transects 1 to 119
    >>> print(f"Lower bounds: {lower[:3]}")  # First 3 values
    Lower bounds: [1.1, 2.1, 3.1]

    Notes
    -----
    The boundary values are encoded with decimal places that indicate spatial positions:
    - .1 indicates western boundary (for regions 1 and 3)
    - .4 indicates eastern boundary (for regions 1 and 3)
    - .6 indicates southern boundary (for region 2)
    - .9 indicates northern boundary (for region 2)
    """
    
    # Initialize variables
    transect_start = None
    transect_end = None
    transect_lower_bound = [] # W/S
    transect_upper_bound = [] # E/N

    # Region 1: parallel transects to latitudes from south of SCB to west of Haida Gwaii 
    if region == 1:
        # ---- Southern-most transect
        transect_start = 1
        # ---- Northern-most transect
        transect_end = 119
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)] 
    # Region 2: transects parallel to longitudes north of Haida Gwaii
    elif region == 2:
        # ---- Western-most transect
        transect_start = 121
        # ---- Eastern-most transect
        transect_end = 127
        # ---- Southern boundary
        transect_lower_bound = [i + 0.6 for i in range(transect_start, transect_end + 1)]
        # ---- Northern boundary
        transect_upper_bound = [i + 0.9 for i in range(transect_start, transect_end + 1)] 
    # Region 3: parallel transects to latitudes west of Haida Gwaii
    else:
        # ---- Southern-most transect
        transect_start = 129
        # ---- Northern-most transect
        transect_end = 145
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)] 

    return transect_start, transect_end, transect_lower_bound, transect_upper_bound
