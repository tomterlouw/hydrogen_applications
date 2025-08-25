# Define a custom function for autopct to conditionally display percentages and amounts
import numpy as np
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.patches as mpatches  # Import for the legend patches
import matplotlib.lines as mlines  # Import for circle patches in the legend
import matplotlib.patches as patches

import seaborn as sns
import math

def custom_autopct(pct, all_vals, threshold_1=10, threshold_2=6, divider = 1):
    absolute = int(np.round(pct / 100. * np.sum(all_vals)))  # Calculate absolute value in kg
    if pct > threshold_1:  # Show amount in Mt for percentages greater than 5%
        return f'{pct:.0f}%\n{absolute / divider:.0f} Mt'  # Convert kg to Mt
    elif pct > threshold_2:  # Show percentage only if greater than 2%
        return f'{pct:.0f}%'  # Show percentage
    else:
        return ''  # Return empty string for percentages <= 2%

def plot_decarbonization_potential_map(df, end_product_colors, output_filepath, threshold=1, leg_pos = (68.08, -8.8)):
    """
    Plot decarbonization potential maps with pie charts and annotations.

    Parameters:
    - df: pd.DataFrame
        DataFrame with columns ['Country Full', 'lat', 'lon', 'End product', 'color', 'decarbonization_potential_Mt', 'decarbonization_potential_Mt_low']
    - end_product_colors: dict
        Dictionary mapping end product names to colors (must correspond to 'color' column in df)
    - output_filepath: str
        File path to save the resulting figure
    """

    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(15, 10),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )
    plt.subplots_adjust(hspace=0.09)

    for ax, pot_col, title in zip(
        axs,
        ['decarbonization_potential_Mt', 'decarbonization_potential_Mt_low'],
        [r'$\mathbf{a}$ Emission reduction potential: business-as-usual',
         '$\mathbf{b}$ Emission reduction potential: low-carbon']
    ):

        # Add map features
        ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
        ax.add_feature(cfeature.OCEAN, facecolor='#e0e0e0')
        ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='gray')

        # Aggregate data by country
        df_agg = df.groupby('Country Full').agg({
            'lat': 'mean', 'lon': 'mean', pot_col: 'sum'
        }).reset_index()

        # Plot pie charts on map for each country
        for _, row in df_agg.iterrows():
            country_data = df[df['Country Full'] == row['Country Full']]
            if row[pot_col] > threshold:
                grouped_df = country_data.groupby('End product').agg({
                    pot_col: 'sum', 'color': 'first'
                }).reset_index()
                grouped_df = grouped_df[grouped_df[pot_col] > 0]

                pie_ax = inset_axes(
                    ax,
                    width=row[pot_col]/130,
                    height=row[pot_col]/130,
                    loc='center',
                    bbox_to_anchor=(row['lon'], row['lat']),
                    bbox_transform=ax.transData,
                    borderpad=0
                )

                wedges = pie_ax.pie(
                    grouped_df[pot_col],
                    colors=grouped_df['color'],
                    startangle=90,
                    wedgeprops={'linewidth': 0.5, 'edgecolor': 'k'}
                )
                pie_ax.set_aspect('equal')

        # Add country labels with conditions
        for _, row in df_agg.iterrows():
            if (row[pot_col] > 60 if pot_col == 'decarbonization_potential_Mt' else row[pot_col] > 13.5):
                if row["Country Full"] in ["United Kingdom", "Brazil"] and pot_col == 'decarbonization_potential_Mt':
                    ax.text(
                        row['lon'], row['lat'] + 17,
                        f"{row['Country Full']}:\n {int(round(row[pot_col], 0))} Mt/a",
                        fontsize=9, ha='center', va='center', zorder=1000,
                        color='black',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3')
                    )
                elif row["Country Full"] == "China" and pot_col == 'decarbonization_potential_Mt':
                    ax.text(
                        row['lon'] + 6, row['lat'] - (row[pot_col]**0.92 / 3.43),
                        f"{row['Country Full']}: {int(round(row[pot_col], 0))} Mt/a",
                        fontsize=9, ha='center', va='center', zorder=1000,
                        color='black',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3')
                    )
                elif pot_col == 'decarbonization_potential_Mt':
                    ax.text(
                        row['lon'] + 2, row['lat'] - (row[pot_col]**0.90 / 3.43),
                        f"{row['Country Full']}: {int(round(row[pot_col], 0))} Mt/a",
                        fontsize=9, ha='center', va='center', zorder=1000,
                        color='black',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3')
                    )
                else:
                    ax.text(
                        row['lon'] + 2, row['lat'] - 7,
                        f"{row['Country Full']}: {int(round(row[pot_col], 0))} Mt/a",
                        fontsize=9, ha='center', va='center', zorder=1000,
                        color='black',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3')
                    )

        ax.set_title(title, fontsize=14, loc='left')

        # Add inset pie chart summarizing end product share
        grouped_df = df.groupby('End product').agg({pot_col: 'sum', 'color': 'first'}).reset_index()
        grouped_df = grouped_df[grouped_df[pot_col] >= 0]

        ax_inset = inset_axes(
            ax, width=1.2, height=1.2,
            bbox_to_anchor=(0.94, 0.94, 0.058, 0.058),
            bbox_transform=ax.transAxes,
            borderpad=0
        )

        wedges, texts, autotexts = ax_inset.pie(
            grouped_df[pot_col],
            labels=None,
            colors=grouped_df['color'],
            autopct=lambda pct: custom_autopct(pct, grouped_df[pot_col]),
            startangle=90,
            wedgeprops={'linewidth': 0.8, 'edgecolor': 'black'},
            pctdistance=0.56
        )
        plt.setp(autotexts, size=6, color="white", weight='bold')
        ax_inset.set_title(f'Total: {df_agg[pot_col].sum():.0f} Mt', fontsize=8.1, weight='bold', y=0.850)

        # Add a translucent rectangle behind the pie inset
        rect = patches.Rectangle(
                (0.84, 0.69), 0.2, 0.6,
                transform=ax.transAxes,
                facecolor=(0, 0, 0, 0.1),   # black with 10% opacity
                edgecolor=(0, 0, 0, 0.5),     # solid black border
                linewidth=1
            )
        ax.add_patch(rect)
        
        # ========== EUROPEAN ZOOM INSET ==========
        # European extent coordinates
        # Add zoomed-in inset for Europe
        # Add Europe inset map in bottom-left
        europe_ax = inset_axes(
            ax,
            width="30%", height="42%",  # relative to parent
            loc='lower left',
            bbox_to_anchor=(-0.0535, -0.015, 1, 1),
            bbox_transform=ax.transAxes,
            axes_class=cartopy.mpl.geoaxes.GeoAxes,
            axes_kwargs=dict(projection=ccrs.PlateCarree()),
            #borderpad=1.5
        )
    
        zoom_extent = [-10, 28, 32, 70] # [lon_min, lon_max, lat_min, lat_max]
        europe_ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
        europe_ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
        europe_ax.add_feature(cfeature.OCEAN, facecolor='#e0e0e0')
        europe_ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.3)
        europe_ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='gray')
        europe_ax.set_title("Europe", fontsize=8, loc='left', weight='bold')
    
        mark_inset(ax, europe_ax, loc1=2, loc2=4, fc="none", ec='darkgray', linestyle='')
        # Then access the connecting lines and hide them
        for line in ax.patches[-1].get_children():
            if hasattr(line, 'set_visible'):
                line.set_visible(False)

        # Aggregate data by country
        df_agg = df.groupby('Country Full').agg({
            'lat': 'mean', 'lon': 'mean', pot_col: 'sum'
        }).reset_index()
        
        # Plot pies in the Europe inset
        for _, row in df_agg.iterrows():
            if zoom_extent[0] <= row['lon'] <= zoom_extent[1] and zoom_extent[2] <= row['lat'] <= zoom_extent[3]:  # Filter only Europe
                country_data = df[df['Country Full'] == row['Country Full']]
                if row[pot_col] > threshold:
                    grouped_df = country_data.groupby('End product').agg({
                        pot_col: 'sum', 'color': 'first'
                    }).reset_index()
                    grouped_df = grouped_df[grouped_df[pot_col] > 0]
    
                    pie_ax = inset_axes(
                        europe_ax,
                        width=row[pot_col]/130,
                        height=row[pot_col]/130,
                        loc='center',
                        bbox_to_anchor=(row['lon'], row['lat']),
                        bbox_transform=europe_ax.transData,
                        borderpad=0
                    )
    
                    wedges = pie_ax.pie(
                        grouped_df[pot_col],
                        colors=grouped_df['color'],
                        startangle=90,
                        wedgeprops={'linewidth': 0.5, 'edgecolor': 'k'}
                    )
                    pie_ax.set_aspect('equal')
                    
    # Legend creation
    legend_patches = [
        mpatches.Patch(color=color,
                       label=label.replace('H2', 'H$_2$').replace('CH4', 'CH$_4$'))
        for label, color in end_product_colors.items()
    ]

    size_legend = [mlines.Line2D([], [], marker='o', color='white',
                                markerfacecolor='gray', markersize=size,
                                label=f'{int(size * 2)} Mt') for size in [6.25, 12.5, 25]]

    size_legend = [mlines.Line2D([], [], marker='', color='white', markersize=0, label='')] + size_legend

    all_legend_patches = legend_patches + size_legend
    plt.legend(handles=all_legend_patches, bbox_to_anchor=leg_pos, ncol=6, fontsize=9, frameon=False)

    # Set visible frames on the maps
    for ax in axs.flatten():
        ax.spines['geo'].set_visible(True)
        ax.spines['geo'].set_linewidth(1.5)
        ax.spines['geo'].set_color('gray')

    # Save and show
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    if '.png' in str(output_filepath):
        plt.savefig(output_filepath.replace('figs/', 'figs/pdf/').replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_hydrogen_pie_map(df, end_product_colors, column_plot = 'kg_h2', font_size_leg=15, leg_cols=6,
                          groupeby = 'End product', output_path="figs/h2_prod_map.png", div_factor = 1e9, threshold=0.5):
    df_2 = df.copy()
    df_2[column_plot] = df_2[column_plot] / div_factor
    
    fig, ax = plt.subplots(figsize=(20, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
    ax.add_feature(cfeature.OCEAN, facecolor='#e0e0e0')
    ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='gray')

    df_filtered = df_2.sort_values(by=column_plot, ascending=True)
    df_agg = df_filtered.groupby('Country Full').agg({'lat': 'mean', 'lon': 'mean', column_plot : 'sum'}).reset_index()

    for _, row in df_agg.iterrows():
        country_data = df[df['Country Full'] == row['Country Full']]
        if row[column_plot ] > threshold:
            grouped_df = country_data.groupby(groupeby).agg({column_plot : 'sum', 'color': 'first'}).reset_index()
            grouped_df = grouped_df[grouped_df[column_plot ] > 0]

            pie_ax = inset_axes(ax, width=row[column_plot ]/8.5, height=row[column_plot ]/8.5, loc='center',
                                bbox_to_anchor=(row['lon'], row['lat']),
                                bbox_transform=ax.transData, borderpad=0)
            pie_ax.pie(
                grouped_df[column_plot ],
                colors=grouped_df['color'],
                startangle=90,
                wedgeprops={'linewidth': 0.5, 'edgecolor': 'k'}
            )
            pie_ax.set_aspect('equal')

    for _, row in df_agg.iterrows():
        if abs(row[column_plot ]) > 4:
            x_offset, y_offset = 0, - (row[column_plot]**0.64 * 3.9)
            if row['Country Full'] == "United Kingdom":
                x_offset, y_offset = -27, row[column_plot] * 1.7
            elif row['Country Full'] == "Spain":
                x_offset, y_offset = -26, row[column_plot] * 1.5
            elif row['Country Full'] == "China":
                x_offset, y_offset = 2, - (row[column_plot]**0.85 * 2.4)

            ax.text(row['lon'] + x_offset, row['lat'] + y_offset, 
                    f"{row['Country Full']}: {round(row[column_plot], 1)} Mt/a", 
                    fontsize=12.5, ha='center', va='center',
                    color='black', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))

    ax.set_title(r'$\mathbf{a}$ Planned low-carbon hydrogen production and associated use (2043)', fontsize=22, loc='left')
    ax.spines['geo'].set_visible(True)
    ax.spines['geo'].set_linewidth(1.5)
    ax.spines['geo'].set_color('gray')

    # Add inset pie chart (global hydrogen use by application)
    grouped_df = df.groupby([groupeby]).agg({column_plot: 'sum', 'color': 'first'}).reset_index()
    ax_inset = inset_axes(ax, width=2.2, height=2.2, 
                          bbox_to_anchor=(0.934, 0.934, 0.065, 0.065),
                          bbox_transform=ax.transAxes, borderpad=0)
    wedges, texts, autotexts = ax_inset.pie(
        grouped_df[column_plot],
        labels=None,
        colors=grouped_df['color'],
        autopct=lambda pct: custom_autopct(pct, grouped_df[column_plot], threshold_1=10, threshold_2=5, divider = 1e9),
        startangle=90,
        wedgeprops={'linewidth': 1, 'edgecolor': 'black'},
        pctdistance=0.65
    )
    plt.setp(autotexts, size=12, weight="bold", color="white")
    ax_inset.set_title('H$_2$ use by application', fontsize=13.2, weight='bold', y=0.892)

    # Legend setup
    legend_patches = [mpatches.Patch(color=color, label=label.replace('H2', 'H$_2$').replace('CH4', 'CH$_4$')) 
                      for label, color in end_product_colors.items()]
    size_legend = [mlines.Line2D([], [], marker='', color='white', markersize=0, label='')] + \
                  [mlines.Line2D([], [], marker='o', color='white', markerfacecolor='gray', markersize=size, 
                                 label=f'{int(size/6)} Mt') for size in [6.25, 12.5, 25]]
    
    # Add a translucent rectangle behind the pie inset
    rect = patches.Rectangle(
            (0.818, 0.675), 0.20, 0.6,
            transform=ax.transAxes,
            facecolor=(0, 0, 0, 0.1),   # black with 10% opacity
            edgecolor=(0, 0, 0, 0.5),     # solid black border
            linewidth=1
        )
    ax.add_patch(rect)

    plt.legend(handles=legend_patches + size_legend, bbox_to_anchor=(0.94, -1.8), 
               ncol=leg_cols, fontsize=font_size_leg, frameon=False)
    
    # ========== EUROPEAN ZOOM INSET ==========
    # European extent coordinates
    # Add zoomed-in inset for Europe
    # Add Europe inset map in bottom-left
    europe_ax = inset_axes(
        ax,
        width="30%", height="42%",  # relative to parent
        loc='lower left',
        bbox_to_anchor=(-0.05, -0.01, 1, 1),
        bbox_transform=ax.transAxes,
        axes_class=cartopy.mpl.geoaxes.GeoAxes,
        axes_kwargs=dict(projection=ccrs.PlateCarree()),
        #borderpad=1.5
    )

    zoom_extent = [-10, 28, 32, 70] # [lon_min, lon_max, lat_min, lat_max]
    europe_ax.set_extent(zoom_extent, crs=ccrs.PlateCarree())
    europe_ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')
    europe_ax.add_feature(cfeature.OCEAN, facecolor='#e0e0e0')
    europe_ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.3)
    europe_ax.add_feature(cfeature.COASTLINE, linewidth=0.5, color='gray')
    europe_ax.set_title("Europe", fontsize=12, loc='left', weight='bold')

    mark_inset(ax, europe_ax, loc1=2, loc2=4, fc="none", ec='darkgray', linestyle='')
    # Then access the connecting lines and hide them
    for line in ax.patches[-1].get_children():
        if hasattr(line, 'set_visible'):
            line.set_visible(False)
    
    # Plot pies in the Europe inset
    for _, row in df_agg.iterrows():
        if zoom_extent[0] <= row['lon'] <= zoom_extent[1] and zoom_extent[2] <= row['lat'] <= zoom_extent[3]:  # Filter only Europe
            country_data = df[df['Country Full'] == row['Country Full']]
            if row[column_plot] > threshold:
                grouped_df = country_data.groupby(groupeby).agg({column_plot : 'sum', 'color': 'first'}).reset_index()
                grouped_df = grouped_df[grouped_df[column_plot] > 0]
    
                pie_ax = inset_axes(europe_ax, width=row[column_plot]/8.5, height=row[column_plot]/8.5, loc='center',
                                    bbox_to_anchor=(row['lon'], row['lat']),
                                    bbox_transform=europe_ax.transData, borderpad=0)
                pie_ax.pie(
                    grouped_df[column_plot],
                    colors=grouped_df['color'],
                    startangle=90,
                    wedgeprops={'linewidth': 0.5, 'edgecolor': 'k'}
                )
                pie_ax.set_aspect('equal')

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    if '.png' in str(output_path):
        plt.savefig(output_path.replace('figs/', 'figs/pdf/').replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_decarbonization_potential_boxplots(results_df_wo_ff, agg_df, replace_low_carbon_product, sec_y_axis=False,
                                            output_path='figs/boxplot_decarbonization_potential.png'):
    font_size = 23

    # Melt data
    melted_df = pd.melt(
        results_df_wo_ff,
        id_vars=['End product', 'Unit'],
        value_vars=['decarbonization_potential_tCO2_tH2', 'decarbonization_potential_tCO2_tH2_low'],
        var_name='Scenario',
        value_name='Decarbonization Potential'
    )

    scenario_map = {
        'decarbonization_potential_tCO2_tH2': 'Business\n-as-usual',
        'decarbonization_potential_tCO2_tH2_low': 'Low-carbon'
    }
    melted_df['Scenario'] = melted_df['Scenario'].map(scenario_map)

    end_products = sorted(results_df_wo_ff['End product'].unique())
    n_cols, max_y = 7, 50
    n_rows = math.ceil(len(end_products) / n_cols)
    scaling_factor = results_df_wo_ff.kg_h2.sum() / 1e12  # Gt

    colors = {
        'Business\n-as-usual': '#CC0000',
        'Low-carbon': '#2da095'
    }

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 14))
    axes = axes.flatten()
    sns.set(style="whitegrid")

    for i, product in enumerate(end_products):
        df_prod = melted_df[melted_df['End product'] == product]

        n_val = (
                        df_prod.groupby(['End product', 'Scenario'])['Decarbonization Potential']
                        .count()
                        .reset_index(name='n')
                    )['n'][0]
        ax = axes[i]

        if product == "Biofuels":
            divider = make_axes_locatable(ax)
            ax_top = divider.new_vertical(size="50%", pad=0.2)
            fig.add_axes(ax_top)

            ax.set_ylim(-max_y, 13)
            ax_top.set_ylim(40, 90)
            ax_top.set_xticks([])

            for sp in ['top', 'bottom']:
                ax.spines[sp].set_visible(True)
                ax_top.spines[sp].set_visible(True)
                
            d = .02
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
            ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
            kwargs['transform'] = ax_top.transAxes
            ax_top.plot((-d, +d), (-d, +d), **kwargs)
            ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

            for a in [ax, ax_top]:
                sns.boxplot(x='Scenario', y='Decarbonization Potential', data=df_prod,
                            ax=a, palette=colors, width=0.5)
                a.axhline(0, color='black', linestyle='--')
                a.axvline(0.5, color='black', linestyle='--')

            ax.set_yticks(np.arange(-max_y, 14, 20))
            ax.set_yticklabels(np.arange(-max_y, 14, 20), fontsize=font_size - 6)
            ax_top.set_yticks(np.arange(40, 91, 20))
            ax_top.set_yticklabels(np.arange(40, 91, 20), fontsize=font_size - 6)
            ax.set_xticklabels(list(scenario_map.values()), rotation=90, ha='center', fontsize=font_size - 3)
            ax_top.set_title(f'$\mathbf{{{chr(97 + i)}}}$ {product}', fontsize=font_size, y=1.06)
            
            ax_top.annotate(replace_low_carbon_product[product], xy=(0.55, 82), fontsize=font_size - 8)

            for scen_code, scen_label in scenario_map.items():
                mean_val = agg_df.loc[product][scen_code].item()
                if mean_val <= 13:
                    ax.annotate(f'{mean_val:.1f}', xy=(scen_label, mean_val), xytext=(scen_label, mean_val + 2),
                                ha='center', fontsize=font_size - 2, weight='bold')
                elif mean_val >= 40:
                    ax_top.annotate(f'{mean_val:.1f}', xy=(scen_label, mean_val), xytext=(scen_label, mean_val + 2),
                                    ha='center', fontsize=font_size - 2, weight='bold')

            ax_top.set_ylabel('', fontsize=font_size)

            if sec_y_axis:
                sec_ax = ax_top.twinx()
                sec_ax.set_ylim(ax_top.get_ylim()[0] * scaling_factor, ax_top.get_ylim()[1] * scaling_factor)
                sec_ax.tick_params(axis='y', labelsize=font_size - 6)
                sec_ax.grid(False)
                ax_top.set_ylabel('', fontsize=font_size)

        else:
            sns.boxplot(x='Scenario', y='Decarbonization Potential', data=df_prod,
                        ax=ax, palette=colors, width=0.6)
            ax.set_ylim(-max_y, max_y)
            ax.set_yticks(np.arange(-max_y, max_y + 1, 20))
            ax.set_yticklabels(np.arange(-max_y, max_y + 1, 20), fontsize=font_size - 6)
            ax.set_xticklabels(list(scenario_map.values()), rotation=90, ha='center', fontsize=font_size - 3)
    
            ax.axhline(0, color='black', linestyle='--')
            ax.axvline(0.5, color='black', linestyle='--')

            title = product.replace("H2", "Hydrogen").replace("CH4", "CH$_4$").replace("Other Ind", "High-temperature heat")
            ax.set_title(f'$\mathbf{{{chr(97 + i)}}}$ {title}', fontsize=font_size, y=1.02)

            for scen_code, scen_label in scenario_map.items():
                mean_val = agg_df.loc[product][scen_code].item()
                ax.annotate(f'{mean_val:.1f}', xy=(scen_label, mean_val), xytext=(scen_label, mean_val + 2),
                            ha='center', fontsize=font_size - 2, weight='bold')

            offset = max_y * 0.80 if product == 'Power' else max_y * 0.9
            h_x=0.52 if product == 'Power' else 0.55
            ax.annotate(replace_low_carbon_product[product], xy=(h_x, offset), 
                        fontsize=font_size - 10 if product == 'Power' else font_size - 8)

            if sec_y_axis:
                if i != 20:
                    sec_ax = ax.twinx()
                    sec_ax.set_ylim(ax.get_ylim()[0] * scaling_factor, ax.get_ylim()[1] * scaling_factor)
                    sec_ax.tick_params(axis='y', labelsize=font_size - 6)
                    sec_ax.grid(False)

        ax.set_ylabel('', fontsize=font_size)
        ax.set_xlabel('')

        ax.grid(False) 

        # Add n annotations for each box
        for j, scen_label in enumerate(df_prod['Scenario'].unique()):
            if j == 0:
                ax.text(j-0.2, ax.get_ylim()[1] - (max_y * 0.01), f"n={n_val}",
                        ha='center', va='top', fontsize=font_size - 10)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.text(0.01, 0.73, 'Emission reduction potential (tCO$_2$-eq./tH$_2$)', 
             ha='center', va='center', fontsize=font_size-2, rotation=90)
    fig.text(0.01, 0.28, 'Emission reduction potential (tCO$_2$-eq./tH$_2$)', 
             ha='center', va='center', fontsize=font_size-2, rotation=90)

    if sec_y_axis:
        fig.text(0.98, 0.73, 'Max. emission reduction potential (GtCO$_2$-eq./a)', 
                 ha='center', va='center', fontsize=font_size - 4, rotation=90)
        fig.text(0.98, 0.28, 'Max. emission reduction potential (GtCO$_2$-eq./a)', 
                 ha='center', va='center', fontsize=font_size - 4, rotation=90)

    plt.tight_layout(pad=3)
    plt.subplots_adjust(top=0.9, wspace=0.40)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    if '.png' in str(output_path):
        plt.savefig(output_path.replace('figs/', 'figs/pdf/').replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_violin_end_products(
    results_df,
    agg_df,
    end_products,
    scenario_labels,
    end_product_colors,
    calculate_lca_impact,
    MY_METHODS,
    new_db_name_low_acts,
    hydrogen_impact=False,
    height_plot = 12,
    font_size = 20,
    output_path='figs/violinplot_end_products_cc.png'
):
    """
    Creates violin plots of climate impact for a list of hydrogen end products, comparing application-specific impacts 
    with hydrogen supply-specific impacts and reference scenario averages.

    Parameters
    ----------
    results_df : pd.DataFrame
        DataFrame containing LCA results with columns such as 'End product', 'lca_impact_climate change_unit', 
        'lca_impact_h2_spec', 'Activity_Name_prod', and 'Unit'.
    
    agg_df : pd.DataFrame
        Aggregated DataFrame containing mean reference scenario values for each product with columns 
        'lca_impact_climate change_unit_ref' and 'lca_impact_climate change_unit_ref_low'.
    
    end_products : list of str
        List of end product names to include in the violin plots.
    
    scenario_labels : dict
        Dictionary mapping scenario column names to human-readable labels.
    
    end_product_colors : dict
        Dictionary mapping end product names to plot colors.
    
    calculate_lca_impact : callable
        Function that takes (activity_name, reference_name, country, impact_method, db) and returns a float impact.
    
    MY_METHODS : list
        List of impact methods; the first method is used for calculating minimal hydrogen impact.
    
    new_db_name_low_acts : str
        Name of the database for calculating low-carbon hydrogen impacts.
    
    hydrogen_impact : bool, optional (default=True)
        If True, includes hydrogen supply-specific violin plots and vertical divider.
    
    output_path : str, optional
        File path for saving the output plot.
    """
    n_products = len(end_products)
    n_cols = 7
    n_rows = math.ceil(n_products / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, height_plot))
    sns.set(style="whitegrid")
    axes = axes.flatten()

    linestyles = ["--", "-."]
    colors_lines = ['darkred', 'darkgreen']

    end_products.sort()

    for i, product in enumerate(end_products):
        ax = axes[i]
        sel_product_df = results_df[results_df['End product'] == product]
        n_val = len(sel_product_df)
        unit = sel_product_df['Unit'].unique().item()
        act_name_min = sel_product_df.Activity_Name_prod.unique()[0]

        lca_impact_min = calculate_lca_impact(
            act_name_min,
            'hydrogen, gaseous, 30 bar' if act_name_min == 'nan' else act_name_min,
            "CH",
            impact_method=MY_METHODS[0],
            db=new_db_name_low_acts
        )

        # Application-specific impact violin
        sns.violinplot(
            y='lca_impact_climate change_unit',
            data=sel_product_df,
            ax=ax,
            x=-0.1 if hydrogen_impact else 0,
            color=end_product_colors[product],
            cut=0,
            bw_adjust=0.25,
            width=0.8 if hydrogen_impact else 0.5
        )

        if hydrogen_impact:
            sns.violinplot(
                y='lca_impact_h2_spec',
                data=sel_product_df,
                ax=ax,
                x=0.1,
                color='lightblue',
                cut=0,
                bw_adjust=0.25,
                legend=False
            )

        for j, scenario in enumerate(['lca_impact_climate change_unit_ref', 'lca_impact_climate change_unit_ref_low']):
            scenario_label = scenario_labels[scenario]
            #if hydrogen_impact:
            mean_value = agg_df.loc[product][scenario].item()
            ax.plot([-.4, 0.4], [mean_value, mean_value], color=colors_lines[j],
                    linestyle=linestyles[j], linewidth=2, label=f'{scenario_label}')

            if j == 1:
                ax.plot(
                    lca_impact_min,
                    marker='D',
                    linestyle='',
                    markersize=11,
                    markeredgecolor='red',
                    linewidth=3,
                    markerfacecolor='green',
                    label='Electrolysis using (very) low-carbon power (2040)'
                )

            # to show the n values
            if j==0: 
                if hydrogen_impact==False:
                    ax.text(0.215, ax.get_ylim()[1] * 0.99 if product!='Refining' else ax.get_ylim()[1] * 0.05,
                                        f"n={n_val}", ha='center', va='top', zorder=100,
                                        fontsize=font_size - 6, transform=ax.transData)
                else:
                    ax.text(1.12, ax.get_ylim()[1] - (abs(ax.get_ylim()[1]) * 0.015),
                                        f"n={n_val}", ha='center', va='top', zorder=100,
                                        fontsize=font_size - 6, transform=ax.transData)

        ax.set_title(
            f'$\mathbf{{{chr(97 + i)}}}$ {product.replace("CH4", "CH$_4$").replace("H2", "Hydrogen").replace("Other Ind", "High-temperature heat")}',
            fontsize=font_size - 1
        )
        ax.set_ylabel(f'Climate impact (kgCO$_2$-eq./{unit})', fontsize=font_size - 3)
        ax.set_xlabel('')
        ax.set_ylim(
            min(
                lca_impact_min * 0.8 if lca_impact_min > 0 else lca_impact_min * 1.2,
                sel_product_df['lca_impact_climate change_unit'].min() * 1.3,
                0
            ),
        )

        ax.grid(False) 

        ax.yaxis.set_tick_params(labelsize=font_size - 4)
        ax.xaxis.set_visible(False)

        if hydrogen_impact:
            ax.axvline(x=0.5, color='black', linestyle='--', linewidth=1)
            ax.text(0.75, -0.16, 'Hydrogen', rotation=90, transform=ax.transAxes,
                    ha='center', va='center', fontsize=font_size - 2)
            ax.text(0.25, -0.18, 'Application', rotation=90, transform=ax.transAxes,
                    ha='center', va='center', fontsize=font_size - 2)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    handles, labels = ax.get_legend_handles_labels()
    legend = fig.legend(
        handles, labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.04),
        ncol=3,
        fontsize=font_size,
        title_fontsize=font_size,
        frameon=False
    )
    legend.set_title('Averages of scenarios')
    legend.get_title().set_fontweight('bold')

    plt.tight_layout(pad=2)
    plt.subplots_adjust(top=0.8, wspace=0.6 if hydrogen_impact else 0.9)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    if '.png' in str(output_path):
        plt.savefig(output_path.replace('figs/', 'figs/pdf/').replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
    plt.show()