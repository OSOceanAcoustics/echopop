"""
Interactive comparison: uncorrected vs. selectivity-corrected length distributions.

Requires: specimen_data_selectivity (DataFrame with 'selectivity_expansion') to be
in scope, along with the standard echopop imports already present in the workflow.

Run after the `count_proportions` block.
"""
import panel as pn
import holoviews as hv
from holoviews import opts

pn.extension()
hv.extension("bokeh")

# ==================================================================================================
# BUILD COMPARISON DATAFRAMES
# ==================================================================================================
def _build_dist(data, weight_col, label, group_cols):
    """
    Aggregate a length distribution and normalise within all group dimensions
    except ``length_bin``.

    Parameters
    ----------
    data : pd.DataFrame
    weight_col : str or None
        Column to sum as the count weight.  ``None`` → row count (unweighted).
    label : str
        ``'Uncorrected'`` or ``'Corrected'``.
    group_cols : list of str
        Grouping columns; must include ``'length_bin'``.
    """
    sub = data.dropna(subset=["length_bin"])
    counts = (
        sub.groupby(group_cols)[weight_col].sum()
        if weight_col
        else sub.groupby(group_cols).size()
    ).rename("count").reset_index()

    norm_keys = [c for c in group_cols if c != "length_bin"]
    counts["proportion"] = counts.groupby(norm_keys)["count"].transform(
        lambda x: x / x.sum()
    )
    counts["length_mid"] = counts["length_bin"].apply(lambda x: float(x.mid))
    counts["type"] = label
    return counts


_base = specimen_data_selectivity.dropna(subset=["length", "weight"])

# ---- By stratum
df_strat = pd.concat([
    _build_dist(_base, None,                    "Uncorrected", ["stratum_ks", "sex", "length_bin"]),
    _build_dist(_base, "selectivity_expansion", "Corrected",   ["stratum_ks", "sex", "length_bin"]),
], ignore_index=True)

# ---- By haul (uid)
df_haul = pd.concat([
    _build_dist(_base, None,                    "Uncorrected", ["uid", "sex", "length_bin"]),
    _build_dist(_base, "selectivity_expansion", "Corrected",   ["uid", "sex", "length_bin"]),
], ignore_index=True)

# ---- Lookup: haul → stratum (for subplot titles)
_haul_to_stratum = (
    _base[["uid", "stratum_ks"]].drop_duplicates()
    .set_index("uid")["stratum_ks"].to_dict()
)

strata = sorted(df_strat["stratum_ks"].unique())
sexes  = sorted(df_strat["sex"].unique())
hauls  = sorted(df_haul["uid"].unique())

# ==================================================================================================
# WIDGETS
# ==================================================================================================
view_sel    = pn.widgets.RadioButtonGroup(
    options=["Stratum", "Haul"], value="Stratum", button_type="default",
)
stratum_sel = pn.widgets.Select(name="Stratum", value=strata[0], options=strata)
haul_sel    = pn.widgets.AutocompleteInput(
    name="Haul (UID)", value=str(hauls[0]),
    options=[str(h) for h in hauls],
    case_sensitive=False, min_characters=0,
)
sex_chk = pn.widgets.CheckBoxGroup(
    name="Sex", value=list(sexes), options=list(sexes), inline=True,
)

# ==================================================================================================
# PLOT BUILDER
# ==================================================================================================
_COLORS = {"Uncorrected": "#6baed6", "Corrected": "#e6550d"}
_DASHES  = {"Uncorrected": "dashed",  "Corrected": "solid"}
_WIDTHS  = {"Uncorrected": 1.5,       "Corrected": 2.2}


def _make_overlay(df, key_col, key_val, selected_sexes, title):
    sub = df[(df[key_col] == key_val) & (df["sex"].isin(selected_sexes))]
    if sub.empty:
        return hv.Curve([], kdims="length_mid", vdims="proportion").opts(
            title=title, width=820, height=420,
        )

    curves = []
    for (label, sex), grp in sub.groupby(["type", "sex"]):
        grp = grp.sort_values("length_mid")
        curves.append(
            hv.Curve(grp, kdims="length_mid", vdims=["proportion", "count"], label=f"{label} – {sex}")
            .opts(
                color=_COLORS[label],
                line_dash=_DASHES[label],
                line_width=_WIDTHS[label],
                alpha=0.9,
                tools=["hover"],
            )
        )

    return hv.Overlay(curves).opts(
        opts.Overlay(
            title=title,
            xlabel="Length (cm)", ylabel="Proportion (within stratum × sex)",
            width=820, height=420,
            legend_position="top_right",
            toolbar="above",
            show_grid=True,
        )
    )


def _plot(view, stratum, haul_str, selected_sexes):
    if not selected_sexes:
        return pn.pane.Markdown("_Select at least one sex._", width=820)

    if view == "Stratum":
        title = f"Length distribution — Stratum {stratum}"
        p = _make_overlay(df_strat, "stratum_ks", stratum, selected_sexes, title)
    else:
        # AutocompleteInput returns a string; cast back to original dtype
        try:
            haul_val = type(hauls[0])(haul_str)
        except (ValueError, TypeError):
            return pn.pane.Markdown(f"_Haul '{haul_str}' not found._", width=820)
        strat_label = _haul_to_stratum.get(haul_val, "?")
        title = f"Length distribution — Haul {haul_val}  (stratum {strat_label})"
        p = _make_overlay(df_haul, "uid", haul_val, selected_sexes, title)

    return pn.pane.HoloViews(p, sizing_mode="stretch_width")


bound_plot = pn.bind(_plot, view_sel, stratum_sel, haul_sel, sex_chk)


@pn.depends(view_sel)
def _controls(view):
    selector = stratum_sel if view == "Stratum" else haul_sel
    return pn.Column(
        pn.pane.Markdown(f"**{'Stratum' if view == 'Stratum' else 'Haul (UID)'}**"),
        selector,
        pn.pane.Markdown("**Sex**"),
        sex_chk,
        width=240,
    )


# ==================================================================================================
# LEGEND / ANNOTATION STRIP
# ==================================================================================================
_legend_md = pn.pane.Markdown(
    "<span style='color:#6baed6; font-weight:bold'>━ ━ &nbsp;Uncorrected</span>"
    "&emsp;&emsp;"
    "<span style='color:#e6550d; font-weight:bold'>───&nbsp;Corrected</span>",
    width=820,
)

# ==================================================================================================
# DASHBOARD
# ==================================================================================================
dashboard = pn.Column(
    pn.pane.Markdown("## Net Selectivity — Length Distribution Comparison"),
    pn.Row(pn.pane.Markdown("**View by:**"), view_sel),
    pn.Row(
        _controls,
        pn.Column(_legend_md, bound_plot),
    ),
    sizing_mode="stretch_width",
)

dashboard.servable()
dashboard
