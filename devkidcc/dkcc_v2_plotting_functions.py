import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scanpy as sc
import pandas as pd
import anndata as ad
import seaborn as sns

import os

def plot_scpred_scores_distribution(adata, scpred_df, score_cols, split_by=None, fig_size=(12, 6)):
    """Generates a boxplot/stripplot combo for prediction scores, optionally split by sample."""
    
    # Data Integration
    df = scpred_df[score_cols].copy().reset_index(drop=True)
    if split_by:
        df[split_by] = adata.obs[split_by].values

    # Calculate Max Score and Types
    df['Max_Value'] = df[score_cols].max(axis=1)
    df['Highest_Type'] = df[score_cols].idxmax(axis=1).str.replace('scpred_', '') 
    
    ordered_types = [col.replace('scpred_', '') for col in score_cols]

    # Plotting Setup
    fig, ax_score = plt.subplots(figsize=fig_size)
    ax_score.set_title(f"Max Score Distribution {f'(Split by {split_by})' if split_by else ''}", fontsize=16)

    if split_by:
        # Boxplot with legend
        sns.boxplot(
            data=df, x='Highest_Type', y='Max_Value', hue=split_by, 
            order=ordered_types, ax=ax_score, showfliers=False, palette='Set3'
        )
        # Stripplot without redundant legend
        sns.stripplot(
            data=df, x='Highest_Type', y='Max_Value', hue=split_by, 
            order=ordered_types, ax=ax_score, dodge=True, alpha=0.3, 
            s=2, palette='dark:black', legend=False 
        )
        # Position Legend
        ax_score.legend(title=split_by, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        sns.boxplot(data=df, x='Highest_Type', y='Max_Value', order=ordered_types, ax=ax_score, showfliers=False, color='white')
        sns.stripplot(data=df, x='Highest_Type', y='Max_Value', order=ordered_types, ax=ax_score, alpha=0.3, s=2, color='black')

    ax_score.set_xticklabels(ax_score.get_xticklabels(), rotation=45, ha='right')
    ax_score.set_ylim(-0.05, 1.05)
    ax_score.set_ylabel("Max scPred Score")
    ax_score.set_xlabel("Predicted Identity")
    
    plt.tight_layout()
    return fig, ax_score


def create_custom_cmap(high_color, low_color='lightgrey'):
    colors = [low_color, high_color]
    return LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

def plot_scpred_umap_py(adata, score_cols, colors, scpred_df=None, split_by=None, fig_size=(10, 8)):
    """
    Args:
        split_by (str): Column name in adata.obs (e.g., 'sample') to split the boxplots.
    """
    # --- 1. Data Validation and Integration ---
    umap_coords = pd.DataFrame(adata.obsm['X_umap'], columns=['umap1', 'umap2'])
    
    if scpred_df is None:
        scpred_df = adata.obs.copy()
    else:
        scpred_df = scpred_df
    
    # Combine UMAP, scPred scores, and the split_by metadata
    df = pd.concat([umap_coords, scpred_df[score_cols].reset_index(drop=True)], axis=1)
    if split_by:
        df[split_by] = adata.obs[split_by].values

    # --- 2. Calculate Max Score and Types ---
    comparison_data = df[score_cols]
    df['Max_Value'] = comparison_data.max(axis=1)
    df['Highest_Type'] = comparison_data.idxmax(axis=1).str.replace('scpred_', '') 

    color_map_dict = {col.replace('scpred_', ''): color for col, color in zip(score_cols, colors)}

    # --- 3. Visualization Setup ---
    fig = plt.figure(figsize=(fig_size[0] * 2.5, fig_size[1])) 
    ax_umap = fig.add_subplot(1, 2, 1)
    ax_umap.set_title("A. UMAP Projection by Max scPred Score", fontsize=16)
    ax_umap.set_xticks([]); ax_umap.set_yticks([])
    for spine in ax_umap.spines.values(): spine.set_visible(False)
    
    norm = Normalize(vmin=0, vmax=1)
    
    # --- 4. UMAP Plotting ---
    # Define these layout variables inside the function before the loop
    cbar_x_start, cbar_width, cbar_y_start, cbar_gap = 0.48, 0.015, 0.88, 0.02
    cbar_total_height = 0.8
    cbar_height_unit = cbar_total_height / len(score_cols)
    
    for i, (original_col_name, high_color) in enumerate(zip(score_cols, colors)):
        cell_type_clean = original_col_name.replace('scpred_', '')
        
        # 1. Subset the data
        df_subset = df[df['Highest_Type'] == cell_type_clean].copy()
        
        if df_subset.empty: 
            continue
            
        # 2. SORT: Ensure the highest scores in this cell type are plotted last (on top)
        df_subset = df_subset.sort_values(by=original_col_name, ascending=True)
            
        # 3. Create colormap and scatter
        cmap = create_custom_cmap(high_color)
        scatter = ax_umap.scatter(
            df_subset['umap1'], 
            df_subset['umap2'],
            c=df_subset[original_col_name],
            cmap=cmap, 
            norm=norm, 
            s=1.5, 
            alpha=0.8
        )
        
        # 4. Colorbar stack (Using the variables defined above)
        cbar_bottom = cbar_y_start - (i + 1) * cbar_height_unit - (i * cbar_gap)
        cax = fig.add_axes([cbar_x_start, cbar_bottom, cbar_width, cbar_height_unit - cbar_gap])
        fig.colorbar(scatter, cax=cax, orientation='vertical')
        cax.set_ylabel(cell_type_clean, rotation=0, ha='left', va='center', fontsize=10)
        # Adjust layout to make sure the new legend isn't cut off
    plt.subplots_adjust(right=0.85)


## Usage Example:    
 #plot_scpred_umap_py(adata_gfp, 
 #                   score_cols = ["scpred_NPC", "scpred_EN", "scpred_DN", "scpred_PN", "scpred_RC", "scpred_UrEp", "scpred_Endo", "scpred_Stroma"],
 #                   colors = ["green", "black", "orange", "red", "blue", "brown", "pink", "yellow"],
 #                   scpred_df=None, split_by=None, fig_size=(10, 8))   
