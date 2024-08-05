from sqlalchemy import create_engine, text, Engine, inspect
import sqlalchemy as sqla
import pandas as pd
from typing import Optional, Literal, Union, List
import numpy as np
from pathlib import Path
import re

def sql_create(connection: sqla.Connection, dataframe: pd.DataFrame, table_name: str, 
               primary_keys: Optional[list] = None):
    """
    Generate a SQL command to create a table with dynamic columns, primary keys, and indices.

    Args:
        table_name (str): The name of the table.
        columns (dict): A dictionary where keys are column names and values are data types.
        primary_keys (list, optional): List of column names to be used as primary keys.

    Returns:
        str: The SQL command to create the table.
    """
    # Generate column definitions
    column_definitions = (
        ",\n".join(f"{col} {SQL_DTYPES[type(dataframe[col][0]).__name__]}" 
                   for col in dataframe.columns)
        )
    
    # Generate primary key definition
    primary_key_definition = ""
    if primary_keys:
        primary_key_definition = f",\nPRIMARY KEY ({', '.join(primary_keys)})"
        
    # Combine all parts into the final SQL command
    create_table_command = f"""
    CREATE TABLE IF NOT EXISTS {table_name} (
        {column_definitions}
        {primary_key_definition}
    );
    """
    
    # Execute
    connection.execute(text(create_table_command.strip()))

def sql_map_tables(connection: sqla.Connection):
    """
    """
    inspector = inspect(connection)
    table_names = inspector.get_table_names()
    # result = connection.execute(text("SELECT name FROM sqlite_master WHERE type='table';"))
    # table_names = result.fetch_all()
    # Extract table names from the results
    # table_names = [name[0] for name in table_names]
    return table_names

def sql_validate(connection: sqla.Connection, table_name: str): 
    """
    Check if a table exists in the database.

    Args:
        connection: SQLAlchemy Connection object.
        table_name (str): The name of the table to check.

    Returns:
        bool: True if the table exists, False otherwise.
    """    
    inspector = inspect(connection)
    return table_name in inspector.get_table_names()

def sql_inspect(connection: sqla.Connection, table_name: str, columns: List[str] = None):
    """
    Get a list of all tables present

    Args:
        connection: SQLAlchemy Connection object.

    Returns:
        list: True if the table exists, False otherwise.
    """  

    # Inspect the columns from the table
    if columns is None:
        # ---- Create 'inspector' for the db file
        inspector = inspect(connection)
        # ---- Retrieve column information
        column_info =  inspector.get_columns(table_name)
        # ---- Format as a dictionary and return the output
        return {col['name']: {k: v for k, v in col.items() if k != 'name'} for col in column_info}
    else: 
        # Inspect unique values in specified columns
        # ---- Create SQL command
        sql_command = f"SELECT DISTINCT {', '.join(columns)} FROM {table_name};"
        # ---- Execute
        table = connection.execute(text(sql_command.strip()))
        # ---- Extract unique values
        unique_values = table.fetchall()
        # ---- Format as a dictionary and return the output
        return (
            {col: list(set(row[idx] for row in unique_values)) for idx, col in enumerate(columns)}
        )

def sql_drop(connection: sqla.Connection, table_name: str):
    """
    """
    connection.execute(text(f"DROP TABLE IF EXISTS {table_name}"))
    
def sql_insert(connection: sqla.Connection, table_name: str, columns: list, dataframe: pd.DataFrame,
               id_columns: Optional[list] = None):
    """
    Insert data into a table.

    Args:
        connection (Connection): The SQLAlchemy Connection instance.
        table_name (str): The name of the table.
        columns (list): List of column names.
        data (list of dict): List of dictionaries containing data to insert or update.
        conflict_columns (list): List of column names to use for conflict resolution.
    """
    
    # Prepare the SQL statement for insertion
    # ---- Check whether `columns` is '*'
    if "*" in columns:
        # ---- Create 'inspector' for the db file
        inspector = inspect(connection)
        # ---- Get the column names from the db file
        columns = [col['name'] for col in inspector.get_columns(table_name)]
    # ---- If not a List
    elif not isinstance(columns, list):
        columns = [columns]
    # ---- Prepare the columns as a string of column names
    column_names = ", ".join(columns)

    # Format `id_columns`
    if id_columns is not None and not isinstance(id_columns, list):
        id_columns = [id_columns]
    
    # Convert the DataFrame into a tuple and then into a string
    # ---- Replace NaN with None
    dataframe = dataframe.replace([np.nan], [None])
    # ---- DataFrame to Tuple
    data_tuple = [tuple(row) for row in dataframe.itertuples(index=False)]
    
    def format_value(x):
        if isinstance(x, str):
            return "'{}'".format(x.replace("'", "''"))
        elif isinstance(x, pd.Timestamp):
            return "'{}'".format(x)
        elif x is None:
            return 'NULL'
        else:
            return str(x)
        
    # ---- Tuple to String
    # data_str = ", ".join(
    #     # f"({', '.join(map(lambda x: f'\'{x}\'' if isinstance(x, str) else str(x), row))})"
    #             f"({', '.join(map(lambda x: f'\'{x}\'' 
    #                                 if isinstance(x, str) or isinstance(x, pd.Timestamp) 
    #                                 else 'NULL' if x is None else str(x), row))})"
    #     for row in data_tuple
    # )
    data_str = ", ".join(f"({','.join(map(lambda x: format_value(x), row))})"  for row in data_tuple)
    
    # Construct the "ON CONFLICT, DO UPDATE SET" if needed
    on_conflict_clause = ""
    if id_columns:
        on_conflict_clause = f"""
        ON CONFLICT ({', '.join(id_columns)})
        DO UPDATE SET {', '.join(f'{col}=excluded.{col}' for col in columns)}
        """
    
    # Construct the SQL query
    sql_command = f"""
    INSERT INTO {table_name} ({column_names})
    VALUES {data_str}
    {on_conflict_clause}
    """    
    
    # Execute
    connection.execute(text(sql_command.strip()))
    
    # Commit
    connection.commit()

def sql_update(connection: sqla.Connection, table_name: str, columns: list, 
               dataframe: Optional[pd.DataFrame] = None, operation: Optional[str] = None, 
               condition: Optional[str] = None):
    """
    Insert data into a table.

    Args:
        connection (Connection): The SQLAlchemy Connection instance.
        table_name (str): The name of the table.
        columns (list): List of column names.
        data (list of dict): List of dictionaries containing data to insert or update.
        conflict_columns (list): List of column names to use for conflict resolution.
    """
    
    # Prepare the SQL statement for insertion
    # ---- Check whether `columns` is '*'
    if "*" in columns:
        # ---- Create 'inspector' for the db file
        inspector = inspect(connection)
        # ---- Get the column names from the db file
        columns = [col['name'] for col in inspector.get_columns(table_name)]
    # ---- If not a List
    elif not isinstance(columns, list):
        columns = [columns]

    # Format the SET command
    # ---- Update column by applying arithmetic between table and dataframe
    if operation is not None and dataframe is not None:
        set_list = [f"{column} = {column} {operation} {dataframe[column].values[0]}" 
                    for column in columns]
    # ---- Update column by applying arithmetic within table
    if dataframe is None and operation is not None:
        # ---- Make sure `operation` is a list
        if not isinstance(operation, list):
            operation = [operation]
        # ---- Break up the columns into their components
        set_list = [f"{column} = {calculation}" for column, calculation in zip(columns, operation)]
    # ---- Update column by setting a defined value
    if dataframe is not None and operation is None:
        set_list = [f"{column} = {dataframe[column].values[0]}" for column in columns]
    # ---- Join the list
    set_clause = ', '.join(set_list)

    # Add the WHERE clause if a parsed condition is provided
    if condition is not None:
        # ---- Parse the conditional string
        parsed_condition = parse_condition(condition)
        set_clause += " WHERE " + parsed_condition

    # Complete the full command
    sql_command = f"UPDATE {table_name} SET {set_clause};"

    # Execute
    connection.execute(text(sql_command.strip()))
    
    # Commit
    connection.commit()

def sql_select(connection: sqla.Connection, table_name: str, 
               columns: Optional[Union[list, str]] = None, 
               condition: Optional[str] = None, 
               output_type: type = pd.DataFrame):

    # Columns
    if columns is None:
        column_names = "*"
    elif isinstance(columns, list) or isinstance(columns, pd.Index):
        column_names = ", ".join(columns)
    else:
        column_names = columns

    # Prepare the columns as a string of column names
    # if isinstance(columns, list):
    #     column_names = ", ".join(columns)
    # else:
    #     column_names = columns

    # Format the SQL command
    # sql_command = f"SELECT {column_names} FROM {table_name};"
    sql_command = f"SELECT {column_names} FROM {table_name}"

    # Add the WHERE clause if a parsed condition is provided
    if condition is not None:
        # ---- Parse the conditional string
        parsed_condition = parse_condition(condition)
        sql_command += " WHERE " + parsed_condition

    # Execute the command 
    table = connection.execute(text(sql_command))

    # Fetch the data from the table
    data = table.fetchall()
    
    # Inspect the table to construct a dictionary of expected datatypes for each column
    table_info = sql_inspect(connection, table_name=table_name)
    # ---- Whittle down the information dictionary to isolate just the column datatypes
    table_dtypes = {col: info['type'] for col, info in table_info.items()}

    # Raise error if `output_type` is invalid
    if output_type not in [pd.DataFrame, np.ndarray, str, tuple]:
        raise TypeError(
            f"Argument `output_type` ([{output_type}]) must be either `str`, `tuple`, "
            f"`pandas.DataFrame`, or `numpy.ndarray`."
        )

    # Format the output 
    # ---- DataFrame
    if output_type is pd.DataFrame:
        # ---- Create DataFrame
        output_df = pd.DataFrame(data, columns=table.keys())
        # ---- Format the expected datatypes
        df_dtypes = {col: SQL_DTYPES[type(dtype).__name__] 
                     for col, dtype in table_dtypes.items() if col in columns }
        # ---- Apply the dtypes
        return output_df.astype(df_dtypes)
    else:
        # ---- Get the datatypes that will correspond to each value of the tuples
        tuple_dtypes = [SQL_DTYPES[type(dtype).__name__] for _, dtype in table_dtypes.items()]
        # ---- Convert the `Row` objects to tuples 
        converted_data = [
            tuple(dtype(value) if value is not None else None
                for value, dtype in zip(row, tuple_dtypes))
            for row in data
        ]
        # ---- String
        if output_type is str:
            return [item[0] for item in converted_data]
        # ---- Array
        elif output_type is np.ndarray:
            return np.array([item[0] for item in converted_data])
        # ---- Tuple
        else:
            return converted_data

SQL_COMMANDS = {
    "create": dict(function=sql_create, args=["table_name", "dataframe", "primary_keys"]),
    "drop": dict(function=sql_drop, args=["table_name"]),
    "insert": dict(function=sql_insert, args=["table_name", "columns", "dataframe", "id_columns"]),
    "inspect": dict(function=sql_inspect, args=["table_name", "columns"]),
    "map": dict(function=sql_map_tables, args=[]),
    "select": dict(function=sql_select, args=["table_name", "columns", "output_type", "condition"]),
    "update": dict(function=sql_update, args=["table_name", "columns", "condition", "operation", 
                                              "dataframe"]),
    "validate": dict(function=sql_validate, args=["table_name"]),
}
        
SQL_DTYPES = {
    'int32': 'INTEGER',
    'int64': 'INTEGER',
    'float64': 'FLOAT',
    "float": "FLOAT",
    "int": "INTEGER",
    'bool': 'BOOLEAN',
    "Interval": "TEXT",
    "Timestamp": "DATETIME",
    'object': 'TEXT',
    "str": "TEXT",
    "FLOAT": float,
    "INTEGER": int,
    "DATETIME": str,
    "TEXT": str,
} 

def sql_group_update(db_file: str,
                     dataframe: pd.DataFrame, 
                     table_name: str,
                     columns: List[str],
                     unique_columns: List[str],
                     id_columns: Optional[List[str]] = None):
    
    # Check for unique values contained within the table
    unique_values = SQL(db_file, "inspect", table_name=table_name, columns=unique_columns)

    # Get the unique values in the table
    table_values = {col: dataframe[col].unique().tolist() for col in unique_columns}

    # Find mismatched indices
    new_indices = {col: list(set(table_values[col]) - set(unique_values[col])) 
                   for col in unique_columns}

    # Filter the DataFrame to include only rows with these missing values
    # ---- Create DataFrame copy
    filtered_df = dataframe.copy()
    # ---- Iterate through the extracted dictionary
    for col, missing_vals in new_indices.items():
        if missing_vals:
            filtered_df = filtered_df[filtered_df[col].isin(missing_vals)]
        else:
            # ---- Drop the values that are not contained within the list
            filtered_df = pd.DataFrame(columns=filtered_df.columns)

    # Insert into the table if not otherwise present
    if not filtered_df.empty: 
        SQL(db_file, "insert", table_name=table_name, id_columns=id_columns, dataframe=filtered_df)

    # Update the table
    # ---- Format the conditional string
    case_statements = []
    for col in columns:
        case_stmt = "CASE"
        for _, row in dataframe.iterrows():
            # Construct the filter condition based on unique_columns
            filter_conditions = ' AND '.join([
                f"{col} = '{row[col]}'" if isinstance(row[col], str) else f"{col} = {row[col]}"
                for col in unique_columns
            ])
            # Add the WHEN condition to the CASE statement
            case_stmt += f" WHEN {filter_conditions} THEN {row[col]}"
        case_stmt += " END"
        case_statements.append(f"{col} = {case_stmt}")

    # Construct the full SQL UPDATE statement
    update_clause = ', '.join(case_statements)

    # Format the SQL COMMAND string
    sql_command = f"""        
    UPDATE {table_name}
    SET {update_clause}
    WHERE ({' OR '.join([
        ' AND '.join([
            f"{col} = '{row[col]}'" if isinstance(row[col], str) else f"{col} = {row[col]}"
            for col in unique_columns
        ])
        for _, row in dataframe.iterrows()
    ])});
    """

    # Create engine
    engine = create_engine(f"sqlite:///{db_file}")

    # Execute and commit
    with engine.connect() as connection:
        connection.execute(text(sql_command))
        connection.commit()

    # Dispose engine
    engine.dispose()

def get_table_key_names(db_file: Path, data_dict: dict, table_name: str) -> List[str]:

    # Get the data input column names
    if data_dict[table_name].empty:
        # ---- Inspect the table
        inspected_table = SQL(db_file, "inspect", table_name=table_name)
        # ---- Create a list of the data columns
        table_columns = list(inspected_table.keys())
    else:
        # ---- Get the DataFrame column names
        table_columns = data_dict[table_name].columns

    # Create a list of the primary keys
    key_columns = (
           set(table_columns)
           .intersection(["trawl_partition", "sex", "haul_num", "species_id", "longitude", 
                          "latitude", "stratum"]) 
        )

    # Return a list of the output
    return list(key_columns)

def parse_condition(condition: str):
    # Replace logical operators with SQL equivalents
    condition = condition.replace('&', ' AND ').replace('|', ' OR ')
    
    # Handle "IN" lists and replace square brackets with parentheses
    condition = re.sub(r'(\w+)\s*IN\s*\[(.*?)\]', lambda m: f"{m.group(1)} IN ({m.group(2)})", condition, flags=re.IGNORECASE)
    
    # Handle range conditions for BETWEEN, including floats
    condition = re.sub(r'(\d*\.\d+|\d+)\s*<=\s*(\w+)\s*<=\s*(\d*\.\d+|\d+)', 
                       lambda m: f"{m.group(2)} BETWEEN {m.group(1)} AND {m.group(3)}", condition)
    
    # Handle individual comparisons
    condition = re.sub(r'(\w+)\s*([<>!=]+)\s*(\d*\.\d+|\d+)', lambda m: f"{m.group(1)} {m.group(2)} {m.group(3)}", condition)
    condition = re.sub(r'(\w+)\s*([<>!=]+)\s*(\'[^\']*\')', lambda m: f"{m.group(1)} {m.group(2)} {m.group(3)}", condition)

    # Return the parsed condition
    return condition

def format_sql_select(table_name, column_names, condition_string):
    # Base SQL command to select columns from the table
    sql_command = f"SELECT {column_names} FROM {table_name}"
    
    # Parse the condition string
    parsed_condition = parse_condition(condition_string)
    
    # Add the WHERE clause if a parsed condition is provided
    if parsed_condition:
        sql_command += " WHERE " + parsed_condition
    
    # Add a semicolon at the end of the SQL command
    sql_command += ";"
    
    return sql_command

def format_sql_columns(kwargs: dict):

    # Columns
    if "columns" in kwargs and "condition" not in kwargs:
        if isinstance(kwargs["columns"], list) or isinstance(kwargs["columns"], pd.Index):
            kwargs["columns"] = ", ".join(kwargs["columns"])
    elif "columns" not in kwargs: 
        kwargs["columns"] = "*"

    # ID/Conflict columns
    if "id_columns" in kwargs:
        if isinstance(kwargs["id_columns"], list) or isinstance(kwargs["id_columns"], pd.Index):
            kwargs["id_columns"] = ", ".join(kwargs["id_columns"])      

    # Return the updated `kwargs` dictionary
    return kwargs

# TODO: Documentation
def query_processed_files(root_directory: Path, file_settings: dict, files: List[Path]) -> dict:

    # Get the database name 
    db_name = file_settings["database_name"]

    # Create filepath to the SQL database
    # ---- Create Path to SQL database file
    db_directory = root_directory / "database"
    # ---- Create the directory if it does not already exist
    db_directory.mkdir(parents=True, exist_ok=True)
    # ---- Complete path to the database file
    db_file = db_directory / db_name

    # Create a list of string-formatted Path names
    files_str = [str(file) for file in files]
    # ---- Create DataFrame
    current_files = pd.DataFrame(files_str, columns=["filepath"])

    # Check for the table `files_read`
    files_read_tbl = SQL(db_file, "validate", table_name="files_read")

    # Validate whether the table exists; if not, create the table and then insert
    if not files_read_tbl:
        # ---- Create table
        SQL(db_file, "create", table_name="files_read", dataframe=current_files, 
            primary_keys = ["filepath"])
        # ---- Populate table
        SQL(db_file, "insert", table_name="files_read", dataframe=current_files)
        # ---- Break early
        return files_str, db_file
    
    # Query already existing files
    previous_files = SQL(db_file, "select", table_name="files_read", output_type=str)
    # ---- Insert file list
    SQL(db_file, "insert", table_name="files_read", dataframe=current_files, id_columns=["filepath"])

    # Filter out previously processed files
    # ---- Apply filter by comparing sets and return the output
    return list(set(files_str) - set(previous_files)), db_file

# TODO: Documentation
def sql_data_exchange(database_file: Path, **kwargs):

    # Check whether the `table_name` table exists
    table_exists = SQL(database_file, "validate", **kwargs)

    # If empty and table does not exist
    if kwargs["dataframe"].empty and table_exists:
        return SQL(database_file, "select", **kwargs)

    # Create table if it does not exist and run the initial insertion
    if not table_exists:
        # ---- Create table
        SQL(database_file, "create", **kwargs)
        # ---- Insert into table        
        SQL(database_file, "insert", **kwargs)
        # ---- Return the initial dataframe
        return kwargs.get("dataframe")
    
    # Insert into the table
    SQL(database_file, "insert", **kwargs)
    
    # Select existing data frame the database and return the output
    return SQL(database_file, "select", **kwargs)

def reset_db_files(file_configuration: dict, table_exception: Optional[Union[str, List[str]]] = None):

    # Get all database files
    database_files = file_configuration["database"]

    # Iterate through all keys
    for _, db_file in database_files.items():
        # ---- Map the table names
        table_names = SQL(db_file, "map")
        # ---- Drop any noted exceptions
        if not isinstance(table_exception, list):
            table_exception = [table_exception]
        # ---- Drop exception table name
        if None not in table_exception:
            table_names = list(set(table_names) - set(table_exception))
        # ---- Iterate through
        for table_name in table_names:    
            SQL(db_file, "drop", table_name=table_name)
        # ---- Validate that all tables were removed  
        remaining_tables = SQL(table_names, "map")      
        if set(table_names).intersection(set(remaining_tables)):
            raise ValueError(
                f"Attempted reset of [{str(db_file)}] failed."
            )


# TODO: Documentation
def SQL(db_file: str, command: str, **kwargs):

    # Create engine from `db_file` string
    engine = create_engine(f"sqlite:///{db_file}")
    
    # Format the data columns, if necessary, to fit within the SQL commands
    if command not in ["inspect", "update"]:
        kwargs = format_sql_columns(kwargs)
    
    # Run the command
    try:
        with engine.connect() as connection:
            # ---- Get the function name
            command_function = SQL_COMMANDS[command]["function"]
            # ---- Get the function arguments
            command_args = SQL_COMMANDS[command]["args"]
            # ---- Drop unnecessary keys (update `kwargs`)
            kwargs = {key: value for key, value in kwargs.items() if key in command_args}
            # ---- Return output
            return command_function(connection, **kwargs)
    finally: 
        # ---- Dispose of the engine to release any resources being pooled/used
        engine.dispose()
