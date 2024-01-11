from pathlib import Path
from openpyxl import load_workbook

def validate_data_columns( file: Path ,
                           sheet_name: str ,
                           datalayer: str ,
                           validation_settings: dict ):
    """
    Opens a virtual instance of each .xlsx file to validate the presence 
    of require data column/variable names

    Parameters
    ----------
    file: Path
        The file name without the prepended file path
    sheet_name: str
        The Excel sheet name containing the target data
    datalayer: str
        The name of the dataset, such as 'length' or 'strata'
    validation_settings: dict
        The subset CONFIG_MAP settings that contain the target column names
    """
    
    # Open connection with the workbook and specific sheet 
    # This is useful for not calling the workbook into memory and allows for parsing
    # only the necessary rows/column names 
    try:
        workbook = load_workbook(file, read_only=True)
        sheet = workbook[sheet_name]

        # Validate that the expected columns are contained within the parsed 
        # column names of the workbook   
        if datalayer == "vario_krig_para":
            data_columns = [list(row) for row in zip(*sheet.iter_rows(values_only=True))][0]
        else:   
            data_columns = {col.value for col in sheet[1]}

        # Close connection to the work book
        workbook.close()

        # Error evaluation and print message (if applicable)
        if not set(validation_settings.keys()).issubset(set(data_columns)):
            missing_columns = set(validation_settings.keys()) - set(data_columns)
            raise ValueError(f"Missing columns in the Excel file: {missing_columns}")

    except Exception as e:
        print(f"Error reading file '{str(file)}': {e}")