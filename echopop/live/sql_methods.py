from sqlalchemy import create_engine, text, Engine, inspect
import sqlalchemy as sqla
import pandas as pd
from typing import Optional

def sql_create(connection: sqla.Connection, df: pd.DataFrame, table_name: str, 
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
        ",\n".join(f"{col} {SQL_DTYPES[type(col).__name__]}" 
                   for col in df.columns)
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

def sql_inspect(connection: sqla.Connection):
    """
    Get a list of all tables present

    Args:
        connection: SQLAlchemy Connection object.

    Returns:
        list: True if the table exists, False otherwise.
    """  
    return inspect(connection).get_table_names()

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
    # ---- If not a List
    if not isinstance(columns, list):
        columns = list(columns)
        
    column_names = ", ".join(columns)
    
    # Convert the DataFrame into a tuple and then into a string
    # ---- DataFrame to Tuple
    data_tuple = [tuple(row) for row in dataframe.itertuples(index=False)]
    # ---- Tuple to String
    if dataframe.columns.size == 1:
        data_str = ", ".join(
            f"{', '.join(map(str, row))}"
            for row in data_tuple
        )
    else:   
        data_str = ", ".join(
            f"({', '.join(map(str, row))})"
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
    VALUES ({data_str})
    {on_conflict_clause}
    """    

    # Execute
    connection.execute(text(sql_command.strip()))
    
    # Commit
    connection.commit()
    

    
SQL_COMMANDS = {
    "create": sql_create,
    "drop": sql_drop,
    "inspect": sql_inspect,
    "validate": sql_validate,
    
    
    
    "check": "SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';",
    "drop": "DROP TABLE IF EXISTS {table_name};",
    "select": "SELECT {columns} FROM {table_name};",
    "index": "CREATE UNIQUE INDEX IF NOT EXISTS {index_name} ON {table_name} ({columns})",
    # "insert": "INSERT INTO {table_name} ({columns});",
    "insert": """
        INSERT INTO {table_name} ({columns}) 
        SELECT {columns} 
        FROM (SELECT VALUES {values} FROM (VALUES {value_placeholder})) AS source ({columns}) 
        {filter_clause};
        """,
    "inspect": None,
}

SQL_DTYPES = {
    'int32': 'INTEGER',
    'int64': 'INTEGER',
    'float64': 'FLOAT',
    'bool': 'BOOLEAN',
    'datetime64[ns]': 'DATETIME',
    'object': 'TEXT',
    "str": "TEXT",
}

def format_sql_columns(kwargs: dict):
    if "columns" in kwargs:
        if isinstance(kwargs["columns"], list):
            kwargs["columns"] = ", ".join(kwargs["columns"])
    else: 
        kwargs["columns"] = "*"
        
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
            # ---- SELECT
            if command == "select":
                return pd.read_sql(text(SQL_COMMANDS[command].format(**kwargs)), con=connection)
            # ---- REPLACE
            elif command == "replace":
                # ---- Extract dataframe
                df_to_add = kwargs["dataframe"]
                # ---- Replace current
                df_to_add.to_sql(name=kwargs["table_name"], 
                                 con=connection, 
                                 if_exists="replace", index=False)

            # ---- INSERT
            elif command == "insert": 
                # ---- Extract dataframe
                df_to_add = kwargs["dataframe"]
                # ---- Insert into the table
                df_to_add.to_sql(name=kwargs["table_name"], con=connection, if_exists="append", 
                                 index=False)
            # ---- INSPECT
            elif command == "inspect":
                return inspect(engine).get_table_names()
            # ---- OTHER COMMAND
            else: 
                connection.execute(text(SQL_COMMANDS[command].format(**kwargs)))
    finally: 
        # ---- Dispose of the engine to release any resources being pooled/used
        engine.dispose()
