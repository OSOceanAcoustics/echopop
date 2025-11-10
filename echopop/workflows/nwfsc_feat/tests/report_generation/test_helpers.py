import numpy as np
import pandas as pd

from echopop.workflows.nwfsc_feat import reporter as feat


def test_initialize_workbook(tmp_excel):
    """
    Workbook initialization test
    """
    wb = feat.initialize_workbook(tmp_excel)
    assert wb is not None
    wb.create_sheet("TestSheet")
    wb.save(tmp_excel)
    assert tmp_excel.exists()


def test_format_file_sheet(tmp_excel):
    """
    File formatting test
    """
    wb = feat.initialize_workbook(tmp_excel)
    ws = feat.format_file_sheet("TestSheet", wb)
    assert ws.title == "TestSheet"


def test_format_table_headers(plotting_heatmap_data):
    """
    Table header test
    """
    headers = feat.format_table_headers(plotting_heatmap_data)
    assert isinstance(headers, list)
    assert headers[0] == "Length (cm)"


def test_append_datatable_rows(plotting_heatmap_data, tmp_excel):
    """
    Data appending test
    """
    wb = feat.initialize_workbook(tmp_excel)
    ws = feat.format_file_sheet("Test", wb)
    df = plotting_heatmap_data.copy()
    df.loc["Subtotal"] = df.sum()
    df.columns = [str(c) if hasattr(c, "left") else c for c in df.columns]
    df.index = [str(i) if hasattr(i, "left") else i for i in df.index]
    feat.append_datatable_rows(ws, df)
    assert ws.max_row > 1


def test_append_table_aggregates(plotting_heatmap_data, tmp_excel):
    """
    Aggregate data appending test
    """
    wb = feat.initialize_workbook(tmp_excel)
    ws = feat.format_file_sheet("Test", wb)
    df = plotting_heatmap_data.copy()
    df.loc["Subtotal"] = df.sum()
    tables_dict = {"male": df, "female": df, "all": df}
    df[1] = 999
    feat.append_table_aggregates(ws, df, "all", tables_dict)
    assert ws.max_row > 0


def test_append_sheet_label(tmp_excel):
    """
    Append sheet label test
    """
    wb = feat.initialize_workbook(tmp_excel)
    ws = feat.format_file_sheet("Test", wb)
    feat.append_sheet_label(ws, "Title {SEX}", "male")
    assert ws.max_row > 0


def test_pivot_aged_dataframe(sample_ci_grid_data, proportion_dict):
    """
    Pivot aged dataframe test
    """
    geo = sample_ci_grid_data.head(2).copy()
    prop = proportion_dict["aged"].head(2).copy()
    # Add required columns for the test
    geo["biomass"] = [10, 20]
    prop.index = geo.index
    out = feat.pivot_aged_dataframe(geo, prop, "biomass", "all", ["latitude", "longitude"])
    assert isinstance(out, type(geo))


def test_repivot_table(plotting_heatmap_data):
    """
    Repivot test
    """
    df = plotting_heatmap_data.copy()
    out = feat.repivot_table(df, "val")
    assert isinstance(out, type(df))


def test_pivot_aged_weight_proportions():
    """
    Pivot aged-weight proportions test
    """

    # Columns
    columns = pd.MultiIndex.from_product([["female", "male"], [1, 2]], names=["sex", "stratum"])

    # Index
    index = pd.CategoricalIndex(pd.IntervalIndex.from_tuples([(0, 1), (1, 2)]), name="age_bin")

    # Create DataFrame
    df = pd.DataFrame(
        [[10.5, 15.2, 8.3, 12.7], [11.5, 16.2, 9.3, 13.7]], index=index, columns=columns
    )

    # Run test
    out = feat.pivot_aged_weight_proportions(df)
    assert "all" in out
    assert "male" in out
    assert "female" in out


def test_pivot_haul_tables():
    """
    Pivot haul tables test
    """

    # Columns
    columns = pd.MultiIndex.from_product(
        [
            ["female", "male", "all"],
            pd.CategoricalIndex(
                pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)]), name="length_bin"
            ),
        ],
        names=[None, "length_bin"],
    )

    # Index
    index = pd.Index(np.array([1.0, 2.0, 3.0]), name="haul_num")

    # Create DataFrame
    df = pd.DataFrame(
        [np.random.random(9), np.random.random(9), np.random.random(9)],
        index=index,
        columns=columns,
    )

    # Run test
    out_male = feat.pivot_haul_tables(df, "male")
    assert isinstance(out_male, pd.DataFrame)
    assert all(out_male.index == [3.0, 5.0, 7.0, "Subtotal"])
    assert out_male.index.name is None
    assert all(out_male.columns == [1.0, 2.0, 3.0, "Subtotal"])
    assert out_male.columns.name is None


def test_prepare_aged_biomass_dataframes(sample_ci_grid_data):
    """
    Prepare aged biomass dataframes test
    """
    tables = {
        "all": sample_ci_grid_data.copy(),
        "male": sample_ci_grid_data.copy(),
        "female": sample_ci_grid_data.copy(),
    }
    feat.prepare_aged_biomass_dataframes(tables)
    assert "all" in tables


def test_write_aged_dataframe_report(tmp_excel, plotting_heatmap_data):
    """
    Write aged dataframe report test
    """
    df = plotting_heatmap_data.copy()
    tables = {"all": df, "male": df, "female": df}
    sheetnames = {"male": "Male", "female": "Female", "all": "All"}
    feat.write_aged_dataframe_report(tables, tmp_excel, sheetnames)
    assert tmp_excel.exists()


def test_prepare_age_length_tables():
    """
    Prepare age-length tables test
    """

    # Group 1
    columns_1 = pd.Index(["male", "female", "all"])
    index_1 = pd.Index([2, 4, 6])
    df_1 = pd.DataFrame(
        [np.random.random(3), np.random.random(3), np.random.random(3)],
        index=index_1,
        columns=columns_1,
    )

    # Group 2
    columns_2 = pd.MultiIndex.from_product(
        [["female", "male", "all"], [1, 2]],
    )
    index_2 = pd.Index([2, 4, 6])
    df_2 = pd.DataFrame([np.random.random(6)] * 3, index=index_2, columns=columns_2)

    out = feat.prepare_age_length_tables(df_2, df_1)

    assert isinstance(out, dict)
    assert all([k in out for k in ["male", "female", "all"]])
    assert all([isinstance(v, pd.DataFrame) for v in out.values()])


def test_write_age_length_table_report(tmp_excel):
    """
    Write age-length table report
    """

    # Group 1
    columns_1 = pd.Index(["male", "female", "all"])
    index_1 = pd.Index([2, 4, 6])
    df_1 = pd.DataFrame(
        [np.random.random(3), np.random.random(3), np.random.random(3)],
        index=index_1,
        columns=columns_1,
    )

    # Group 2
    columns_2 = pd.MultiIndex.from_product(
        [["female", "male", "all"], [1, 2]],
    )
    index_2 = pd.Index([2, 4, 6])
    df_2 = pd.DataFrame([np.random.random(6)] * 3, index=index_2, columns=columns_2)

    tables = feat.prepare_age_length_tables(df_2, df_1)

    sheetnames = {"male": "Male", "female": "Female", "all": "All"}
    feat.write_age_length_table_report(tables, tmp_excel, sheetnames, "abundance", "Kriged")
    assert tmp_excel.exists()


def test_write_haul_report(tmp_excel):
    """
    Write haul report test
    """

    # Columns
    columns = pd.MultiIndex.from_product(
        [
            ["female", "male", "all"],
            pd.CategoricalIndex(
                pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)]), name="length_bin"
            ),
        ],
        names=[None, "length_bin"],
    )

    # Index
    index = pd.Index(np.array([1.0, 2.0, 3.0]), name="haul_num")

    # Create DataFrame
    df = pd.DataFrame(
        [np.random.random(9), np.random.random(9), np.random.random(9)],
        index=index,
        columns=columns,
    )

    # Get tables
    tables = {sex: feat.pivot_haul_tables(df, sex) for sex in ["male", "female", "all"]}

    # Sheetname mapping
    sheetnames = {"male": "Male", "female": "Female", "all": "All"}

    # Run test
    feat.write_haul_report(tables, sheetnames, tmp_excel)
    assert tmp_excel.exists()
