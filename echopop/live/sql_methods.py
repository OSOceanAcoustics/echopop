from sqlalchemy import create_engine, text, Engine, inspect

import pandas as pd

SQL_COMMANDS = {
    "create": "CREATE TABLE IF NOT EXISTS {table_name} ({column_definitions});",
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
    'object': 'TEXT'
}

# TODO: Documentation
def SQL(db_file: str, command: str, **kwargs):

    # Create engine from `db_file` string
    engine = create_engine(f"sqlite:///{db_file}")

    # Format `columns`, if there are any and more than 1
    if "columns" in kwargs.keys():
        if isinstance(kwargs["columns"], list):
            kwargs["columns"] = ", ".join(kwargs["columns"])
    else:
        kwargs["columns"] = "*"

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
                # ---- Inser into the table
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
