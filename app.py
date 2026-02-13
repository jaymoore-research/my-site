import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

st.set_page_config(
    page_title="Brain Atlas Viewer",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded",
)

DATA_PATH = (
    "/Users/jaymoore/Desktop/brain-atlas/5cfe2ee0-d62a-487c-b0fe-124f39f4df21.h5ad"
)


@st.cache_resource
def load_data():
    """Load and cache the AnnData object"""
    with st.spinner("Loading data..."):
        adata = sc.read_h5ad(DATA_PATH)
        # Ensure obsm coordinates are DataFrames with cell barcodes as index
        for key in adata.obsm_keys():
            adata.obsm[key] = pd.DataFrame(
                adata.obsm[key],
                index=adata.obs.index,
                columns=[f"{key}_{i}" for i in range(adata.obsm[key].shape[1])],
            )
    return adata


def get_available_genes(adata):
    """Get list of available genes"""
    return sorted(adata.var["feature_name"].tolist())


def get_gene_index(adata, gene_name):
    """Get the index of a gene by its name"""
    mask = adata.var["feature_name"] == gene_name
    if mask.any():
        return np.where(mask)[0][0]
    return None


# Load data
try:
    adata = load_data()
    data_loaded = True
except Exception as e:
    st.error(f"Error loading data: {e}")
    data_loaded = False

if data_loaded:
    # Sidebar
    st.sidebar.header("🔬 Data Overview")
    st.sidebar.markdown(f"**Cells:** {adata.n_obs:,}")
    st.sidebar.markdown(f"**Genes:** {adata.n_vars:,}")

    # Embedding selection
    st.sidebar.header("📊 Visualization")
    embedding_options = [k.replace("X_", "") for k in adata.obsm_keys()]
    selected_embedding = st.sidebar.selectbox(
        "Choose embedding", embedding_options, index=0
    )
    embedding_key = f"X_{selected_embedding}"

    # Color by selection
    st.sidebar.header("🎨 Color By")

    # Categorical metadata options
    categorical_cols = adata.obs.select_dtypes(
        include=["category", "object"]
    ).columns.tolist()
    # Remove ID-like columns that aren't useful for coloring
    exclude_cols = ["observation_joinid", "donor_id", "sample_id"]
    categorical_cols = [c for c in categorical_cols if c not in exclude_cols]

    color_by_options = ["None"] + categorical_cols + ["Gene Expression"]
    color_by = st.sidebar.selectbox("Color by", color_by_options, index=1)

    # Get available genes for later use
    available_genes = get_available_genes(adata)

    # If gene expression selected
    gene_name = None
    if color_by == "Gene Expression":
        gene_name = st.sidebar.selectbox("Search gene", available_genes, index=0)

    # Filters
    st.sidebar.header("🔍 Filters")

    # Cell type filter
    if "cell_type" in adata.obs.columns:
        cell_types = ["All"] + sorted(adata.obs["cell_type"].unique().tolist())
        selected_cell_type = st.sidebar.selectbox("Cell Type", cell_types)
    else:
        selected_cell_type = "All"

    # Disease filter
    if "disease" in adata.obs.columns:
        diseases = ["All"] + sorted(adata.obs["disease"].unique().tolist())
        selected_disease = st.sidebar.selectbox("Disease", diseases)
    else:
        selected_disease = "All"

    # Tissue filter
    if "tissue" in adata.obs.columns:
        tissues = ["All"] + sorted(adata.obs["tissue"].unique().tolist())
        selected_tissue = st.sidebar.selectbox("Tissue", tissues)
    else:
        selected_tissue = "All"

    # Apply filters
    mask = pd.Series([True] * adata.n_obs, index=adata.obs.index)

    if selected_cell_type != "All" and "cell_type" in adata.obs.columns:
        mask &= adata.obs["cell_type"] == selected_cell_type

    if selected_disease != "All" and "disease" in adata.obs.columns:
        mask &= adata.obs["disease"] == selected_disease

    if selected_tissue != "All" and "tissue" in adata.obs.columns:
        mask &= adata.obs["tissue"] == selected_tissue

    # Create filtered view
    filtered_cells = adata.obs.index[mask]
    n_filtered = len(filtered_cells)

    # Main content
    st.title("🧠 Brain Atlas Viewer")
    st.markdown(f"Displaying **{n_filtered:,}** of {adata.n_obs:,} cells")

    if n_filtered == 0:
        st.warning("No cells match the selected filters!")
    else:
        # Get coordinates
        coords = adata.obsm[embedding_key].loc[filtered_cells]

        # Create DataFrame for plotting
        plot_df = pd.DataFrame(
            {"x": coords.iloc[:, 0], "y": coords.iloc[:, 1], "cell_id": filtered_cells}
        )

        # Add metadata columns
        for col in categorical_cols:
            plot_df[col] = adata.obs.loc[filtered_cells, col].values

        # Add continuous columns
        continuous_cols = [
            "fraction_mitochondrial",
            "total_genes",
            "total_UMIs",
            "cell_cycle_score",
        ]
        for col in continuous_cols:
            if col in adata.obs.columns:
                plot_df[col] = adata.obs.loc[filtered_cells, col].values

        # Determine coloring
        if color_by == "Gene Expression" and gene_name:
            gene_idx = get_gene_index(adata, gene_name)
            if gene_idx is not None:
                expression = (
                    adata.X[mask, gene_idx].toarray().flatten()
                    if hasattr(adata.X, "toarray")
                    else adata.X[mask, gene_idx]
                )
                plot_df["expression"] = expression
                color_col = "expression"
                color_scale = "Viridis"
                title_suffix = f" - {gene_name} expression"
            else:
                color_col = None
                color_scale = None
                title_suffix = ""
        elif color_by != "None":
            color_col = color_by
            color_scale = None
            title_suffix = f" - colored by {color_by}"
        else:
            color_col = None
            color_scale = None
            title_suffix = ""

        # Create UMAP plot
        fig = px.scatter(
            plot_df,
            x="x",
            y="y",
            color=color_col,
            color_continuous_scale=color_scale
            if color_by == "Gene Expression"
            else None,
            hover_data=["cell_id", "cell_type"]
            if "cell_type" in plot_df.columns
            else ["cell_id"],
            title=f"{selected_embedding} Embedding{title_suffix}",
            opacity=0.6,
            height=600,
        )

        fig.update_traces(marker=dict(size=3))
        fig.update_layout(
            template="plotly_white",
            xaxis_title=f"{selected_embedding}_1",
            yaxis_title=f"{selected_embedding}_2",
        )

        st.plotly_chart(fig, use_container_width=True)

        # Statistics row
        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("📈 Cell Counts")
            if "cell_type" in adata.obs.columns and n_filtered > 0:
                cell_type_counts = (
                    adata.obs.loc[filtered_cells, "cell_type"].value_counts().head(10)
                )
                st.bar_chart(cell_type_counts)

        with col2:
            st.subheader("🧬 QC Metrics")
            if "total_UMIs" in adata.obs.columns and n_filtered > 0:
                st.metric(
                    "Median UMI/cell",
                    f"{adata.obs.loc[filtered_cells, 'total_UMIs'].median():.0f}",
                )
            if "total_genes" in adata.obs.columns and n_filtered > 0:
                st.metric(
                    "Median genes/cell",
                    f"{adata.obs.loc[filtered_cells, 'total_genes'].median():.0f}",
                )
            if "fraction_mitochondrial" in adata.obs.columns and n_filtered > 0:
                st.metric(
                    "Mean % mitochondrial",
                    f"{adata.obs.loc[filtered_cells, 'fraction_mitochondrial'].mean() * 100:.1f}%",
                )

        with col3:
            st.subheader("🏥 Disease Distribution")
            if "disease" in adata.obs.columns and n_filtered > 0:
                disease_counts = adata.obs.loc[filtered_cells, "disease"].value_counts()
                st.write(disease_counts)

        # Gene expression section
        st.header("🧪 Gene Expression Analysis")

        gene_for_violin = st.selectbox(
            "Select gene for violin plot", available_genes, index=0, key="violin_gene"
        )

        if gene_for_violin:
            gene_idx = get_gene_index(adata, gene_for_violin)
            if gene_idx is not None:
                # Get expression values
                expression_full = (
                    adata.X[:, gene_idx].toarray().flatten()
                    if hasattr(adata.X, "toarray")
                    else adata.X[:, gene_idx]
                )

                # Create expression DataFrame
                expr_df = pd.DataFrame(
                    {
                        "expression": expression_full,
                        "cell_type": adata.obs["cell_type"]
                        if "cell_type" in adata.obs.columns
                        else "Unknown",
                        "disease": adata.obs["disease"]
                        if "disease" in adata.obs.columns
                        else "Unknown",
                    }
                )

                # Filter to selected cells
                expr_df = expr_df.loc[filtered_cells]

                # Group by for violin plot
                if "cell_type" in adata.obs.columns:
                    # Get top cell types by count
                    top_cell_types = (
                        expr_df["cell_type"].value_counts().head(15).index.tolist()
                    )
                    expr_df_filtered = expr_df[
                        expr_df["cell_type"].isin(top_cell_types)
                    ]

                    fig_violin = px.violin(
                        expr_df_filtered,
                        x="cell_type",
                        y="expression",
                        color="cell_type",
                        box=True,
                        title=f"{gene_for_violin} Expression by Cell Type (Top 15)",
                        height=500,
                    )
                    fig_violin.update_layout(xaxis_tickangle=-45, showlegend=False)
                    st.plotly_chart(fig_violin, use_container_width=True)

                # Expression statistics
                st.subheader("Expression Statistics")
                expr_stats = (
                    expr_df.groupby("cell_type")["expression"]
                    .agg(["mean", "median", "std", "count"])
                    .round(2)
                )
                expr_stats = expr_stats.sort_values("mean", ascending=False)
                st.dataframe(expr_stats.head(15))

    # Data download option
    st.sidebar.header("💾 Export")
    if st.sidebar.button("Generate Metadata CSV"):
        csv = adata.obs.to_csv(index=True)
        st.sidebar.download_button(
            label="Download metadata",
            data=csv,
            file_name="brain_atlas_metadata.csv",
            mime="text/csv",
        )

else:
    st.error("Failed to load data. Please check the file path.")
    st.stop()
