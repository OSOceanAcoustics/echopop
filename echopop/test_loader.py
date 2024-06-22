def parse_transect_numbers()


configuration_dict = self.config

# Get NASC export settings
export_settings = configuration_dict["nasc_exports"]



# Check that directory exists
# ---- Get export file directory
export_file_directory = export_settings["export_file_directory"]
# ---- Drop prepended "/" or "\\" to avoid using an absolute path
if export_file_directory.startswith('/') or export_file_directory.startswith('\\'):
    export_file_directory = export_file_directory[1:]


# ---- Create string
full_path = Path(configuration_dict["data_root_dir"]) /  export_file_directory
# Check if path exists
if not full_path.exists():
    raise FileNotFoundError(
        f"The export file directory {{{full_path}}} not found!"
    )
# Check if any files are in directory
elif not any(full_path.iterdir()):
    raise FileNotFoundError(
        f"The export file directory {{{full_path}}} contains no files!"
    )

# Generate the intervals, cells, and layers dataframes
intervals_df, cells_df, layers_df = batch_read_echoview_exports(full_path, transect_pattern)

cells_mrg = cells_df.merge(layers_df)

# fix strings
cells_df["region_class"] = ( 
    cells_df["region_class"]
    .replace('"','', regex=True)
    .replace(r'^\s*(.*?)\s*$', r'\1', regex=True)            
)

# ---- convert 'region_class' to lowercase
cells_df["region_class"] = cells_df["region_class"].str.lower()

full_data = group_merge(cells_df, [intervals_df, layers_df])


# ---- 
cells_mrg = cells_df.merge(layers_df)



transect_data = full_data
index_variable: Union[str, List[str]] = ["transect_num", "interval"]
region_column: str = "region_class"
region_filter = region_names["all_ages"]
region_filter = region_names["no_age1"]
integrate_nasc(cells_df, region_filter = region_names["all_ages"])["nasc"].sum()
integrate_nasc(cells_df, region_filter = region_names["no_age1"])["nasc"].sum()

gdat = self.input["acoustics"]["nasc_df"]
gdat[gdat["haul_num"] == 17]
mrg_df
np.unique(mrg_df["region_class_y"])
tt = mrg_df.groupby(["group", "region_class_x"])["nasc"].sum()
tt.reset_index().pivot(columns = "group", index = "region_class_x").replace(np.nan, 0.0).sum()

###
configuration_dict["biological"]

# grab CONFIG
CONFIG_DF = pd.DataFrame(CONFIG_MAP)
name = "strata"
dd = CONFIG_DF.loc[name].dropna()
reference = dd[dd.first_valid_index()]

flat_table = pd.json_normalize(configuration_dict)
config_search_lst = flat_table.filter(regex=r'\bstrata\b')
filtered_files = config_search_lst.filter(regex = "filename")
filtered_sheets = config_search_lst.filter(regex = "sheetname")
file = filtered_files.values.flatten()
sheet = filtered_sheets.values.flatten()

# ---- Get path
haul_key_path =  Path(root_dir) / file[0]
# ---- Validate
if not haul_key_path.exists():
    raise FileExistsError(
        f"The haul-to-transect file {{{str(haul_key_path)}}} does not exist!"
    )
# ---- Read in file and sheet
haul_key = pd.read_excel(haul_key_path, sheet_name = sheet[0])
haul_key_filtered = haul_key.filter(reference.keys())
haul_key_filtered = haul_key_filtered.astype(reference, errors='ignore')

#########
# ---- Merge `region_df` with `haul_key_filtered`
region_haul_df = region_df.merge(haul_key_filtered, on = ["haul_num"])
unique_region_classes = region_haul_df["region_class"].unique()
# Create a boolean mask to filter cells_df based on region_class for each group
mask = cells_df['region_class'].isin(unique_region_classes)
filtered_cells_df = cells_df[mask]


mrg_df = filtered_cells_df.merge(region_df, on = ["region_id", "transect_num"])

full_df = group_merge(cells_df, [intervals_df, layers_df])

for i in np.unique(region_)
region_haul_df = region_df.merge(haul_key_filtered, on = ["haul_num"])
group_merge(intervals_df, [region_df, haul_key_filtered])
# ----
full_df = group_merge(cells_df, [intervals_df, layers_df, region_haul_df])


intervals_df
haul_key_filtered
region_df