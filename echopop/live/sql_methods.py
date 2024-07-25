from sqlalchemy import create_engine, text, Engine, inspect
import sqlalchemy as sqla
import pandas as pd
from typing import Optional

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
    return inspector.get_table_names()

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

def sql_inspect(connection: sqla.Connection, table_name: str):
    """
    Get a list of all tables present

    Args:
        connection: SQLAlchemy Connection object.

    Returns:
        list: True if the table exists, False otherwise.
    """  

    # Create 'inspector' for the db file
    inspector = inspect(connection)

    # Retrieve column information
    column_info =  inspector.get_columns(table_name)

    # Format as a dictionary
    return {col['name']: {k: v for k, v in col.items() if k != 'name'} for col in column_info}


def sql_drop(connection: sqla.Connection, table_name: str):
    """
    """
    connection.execute(text(f"DROP TABLE IF EXISTS {table_name};"))
    
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
    # ---- Identify any possible DATETIME columns
    # datetime_columns = (
    #     {col["name"]: str for col in columns_info 
    #      if isinstance(col["type"], sqla.sql.sqltypes.DATETIME)}
    # )
    # ---- Encapsulate datetimes with quotes by converting to string
    # dataframe = dataframe.astype(datetime_columns)
    # ---- DataFrame to Tuple
    data_tuple = [tuple(row) for row in dataframe.itertuples(index=False)]
    # ---- Tuple to String
    data_str = ", ".join(
        # f"({', '.join(map(lambda x: f'\'{x}\'' if isinstance(x, str) else str(x), row))})"
                f"({', '.join(map(lambda x: f'\'{x}\'' 
                                  if isinstance(x, str) or isinstance(x, pd.Timestamp) 
                                  else 'NULL' if x is None else str(x), row))})"
        for row in data_tuple
    )
    
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

from typing import Literal
import numpy as np
def sql_select(connection: sqla.Connection, table_name: str, columns: list, 
               output_type: type = pd.DataFrame):

    # Prepare the columns as a string of column names
    column_names = ", ".join(columns)

    # Format the SQL command
    sql_command = f"SELECT {column_names} FROM {table_name};"

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
        df_dtypes = {col: SQL_DTYPES[type(dtype).__name__] for col, dtype in table_dtypes.items()}
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
    "inspect": dict(function=sql_inspect, args=["table_name"]),
    "map": dict(function=sql_map_tables, args=[]),
    "select": dict(function=sql_select, args=["table_name", "columns", "output_type"]),
    "validate": dict(function=sql_validate, args=["table_name"]),
}
    
    
SQL_DTYPES = {
    'int32': 'INTEGER',
    'int64': 'INTEGER',
    'float64': 'FLOAT',
    "float": "FLOAT",
    "int": "INTEGER",
    'bool': 'BOOLEAN',
    "Timestamp": "DATETIME",
    'object': 'TEXT',
    "str": "TEXT",
    "FLOAT": float,
    "INTEGER": int,
    "DATETIME": str,
    "TEXT": str,
} 

def format_sql_columns(kwargs: dict):

    # Columns
    if "columns" in kwargs:
        if isinstance(kwargs["columns"], list) or isinstance(kwargs["columns"], pd.Index):
            kwargs["columns"] = ", ".join(kwargs["columns"])
    else: 
        kwargs["columns"] = "*"

    # ID/Conflict columns
    if "id_columns" in kwargs:
        if isinstance(kwargs["id_columns"], list) or isinstance(kwargs["id_columns"], pd.Index):
            kwargs["id_columns"] = ", ".join(kwargs["id_columns"])      

    # Return the updated `kwargs` dictionary
    return kwargs

# TODO: Documentation


# TODO: Documentation
def SQL(db_file: str, command: str, **kwargs):

    # Create engine from `db_file` string
    engine = create_engine(f"sqlite:///{db_file}")
    
    # Format the data columns, if necessary, to fit within the SQL commands
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
            # # ---- SELECT
            # if command == "select":
            #     return pd.read_sql(text(SQL_COMMANDS[command].format(**kwargs)), con=connection)
            # # ---- REPLACE
            # elif command == "replace":
            #     # ---- Extract dataframe
            #     df_to_add = kwargs["dataframe"]
            #     # ---- Replace current
            #     df_to_add.to_sql(name=kwargs["table_name"], 
            #                      con=connection, 
            #                      if_exists="replace", index=False)

            # # ---- INSERT
            # elif command == "insert": 
            #     # ---- Extract dataframe
            #     df_to_add = kwargs["dataframe"]
            #     # ---- Insert into the table
            #     df_to_add.to_sql(name=kwargs["table_name"], con=connection, if_exists="append", 
            #                      index=False)
            # # ---- INSPECT
            # elif command == "inspect":
            #     return inspect(engine).get_table_names()
            # # ---- OTHER COMMAND
            # else: 
            #     connection.execute(text(SQL_COMMANDS[command].format(**kwargs)))
    finally: 
        # ---- Dispose of the engine to release any resources being pooled/used
        engine.dispose()
