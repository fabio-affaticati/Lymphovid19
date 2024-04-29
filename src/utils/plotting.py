# Copyright © 2024 fabio-affaticati

import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import seaborn as sns
import matplotlib.pyplot as plt


FONT = 'JetBrains Mono'
FONTSNS =  'DejaVu Sans'


def volcano_plot(total, plotname, file_name, plotsdir):
    
    minchange = min(total['log2_fold_change']) -1
    maxchange = max(total['log2_fold_change']) +1
    logsig = abs(np.log10(0.05))

    fig = px.scatter(total.drop_duplicates(subset=['cluster']), x='log2_fold_change', y='log10_corrected_p_value_fdr_bh',
                    color='chain',
                    #text = 'cluster',
                    hover_data = ['motif', 'negative', 'positive'],
                    size = 'dot_size',
                    symbol='chain',
                    marginal_x='rug', marginal_y='rug',
                    title=plotname,
                    template='plotly_white',  # Choose a template (e.g., 'plotly', 'plotly_white', 'plotly_dark')
                    color_discrete_map={'A': 'blue', 'B': 'orange'},  # Define colors for categories
                    opacity=0.3,  # Adjust marker opacity
                    #symbol='circle',  # Set marker symbol
                    #size_max=10,  # Set maximum marker size
                    width=1500, height=800
                    )

    fig.update_layout(shapes=[dict(type='line', x0=minchange, x1=maxchange,
                                y0=logsig, y1=logsig,
                                line=dict(color='red', width=1, dash='dash'))])
    
    
    fig.update_traces(marker=dict(line=dict(color='black', width=.5)))
    
    fig.update_layout(font=dict(
            size=14,
            family=FONT,
        ),
        legend_title= "Chain",
        xaxis_title='Log2FoldChange',
        yaxis_title='Log10CorrectedPvalue',
        xaxis_range=[minchange, maxchange],
    )
    #fig.show()
    pio.write_image(fig, plotsdir + "volcano_" + file_name + ".png", engine='kaleido')



def calculate_dot_size(df):
    if df['negative'] == 0:
        return df['positive']
    elif df['positive'] == 0:
        return df['negative']
    else:
        return np.round( abs( max(df['positive'], df['negative']) / min(df['positive'],df['negative']) ) )
    
    
def barplot_unique_sequences(data, metadata, plotsdir):
    
    uniques = pd.DataFrame({'countUniqueCDR3': data.groupby(['sample_id', 'TCR_Chain'])['junction_aa'].nunique()}).reset_index()
    uniques = pd.merge(uniques, metadata, on = 'sample_id', how='left')

    ggplot2_colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00']
    
    fig = px.bar(uniques, x="sample_id", y="countUniqueCDR3", color="TCR_Chain", log_y = True, #text="countUniqueCDR3", 
                facet_row="CONDITION", 
                #pattern_shape="condition", pattern_shape_sequence=[".", "x", "+"],
                width=1500, height=800,
                template='plotly_white', 
                color_discrete_map=dict(zip(uniques['TCR_Chain'].unique(), ggplot2_colors)),
                )
    fig.update_layout(
        title='Unique CDR3 counts',
        xaxis_title='Samples',
        font=dict(family=FONT, size=12),
        #paper_bgcolor='rgb(243, 243, 243)',
        #plot_bgcolor='rgb(243, 243, 243)',
        barmode='group',
        bargap=0.1,
    )

    fig.update_traces(
        textfont_size=12, textangle=0, textposition="outside", cliponaxis=False, marker_line_width=2,
        
    )
    fig.update_yaxes(tickangle=45, tickmode='array', showgrid = False)
    fig.update_xaxes(tickangle=45, showgrid=True)
    #fig.show()
    pio.write_image(fig, plotsdir+"barplot_counts.png", engine='kaleido')
    
    
def scatter_clone_fractions_old(clone_fractions_pivoted, plotsdir):
    
    # clone_fractions_pivoted is a multiiindex dataframe with the following structure:
    # SAMPLE | TIMEPOINTS | clonotype1 | clonotype2 | ... | clonotypeN
    # sample1 | baseline | 0.1 | 0.2 | ... | 0.3
    # sample1 | V1 | 0.2 | 0.3 | ... | 0.4
    # ...
    # sampleN | baseline | 0.3 | 0.4 | ... | 0.5
    # sampleN | V1 | 0.4 | 0.5 | ... | 0.6
    
    # plot a scatter plot of the clone fractions at baseline and V1 (axis x and y, respectively)
    # one trace for every sample (first level of the multiindex)
    
    fig = go.Figure()
    
    for sample in clone_fractions_pivoted.index.get_level_values(0).unique():
        sliced_pivoted_clones = clone_fractions_pivoted.loc[(sample, ['baseline', 'V1']), :]
        try:
            # Add scatter to go.Figure()
            fig.add_trace(go.Scatter(x=sliced_pivoted_clones.loc[sample, 'baseline'],
                                    y=sliced_pivoted_clones.loc[sample, 'V1'],
                                    mode='markers', opacity=0.5,
                                    name=sample))
        except KeyError:
            pass
        
    # Set axis labels
    fig.update_layout(xaxis_title='baseline',
                    yaxis_title='V1')
    # change scale of axis to log
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    # Set title
    fig.update_layout(title='Clone fractions at baseline and V1')
    # Change background color to white
    fig.update_layout(plot_bgcolor='white', font=dict(family=FONT, size=12))
    # change width and height of plot 
    fig.update_layout(width=1000, height=1000)
    #pio.write_image(fig, plotsdir+"clonefractionscatter.png", engine='kaleido')
    fig.write_image(plotsdir+"clonefractionscatter.png", scale = 4)
    
    
    
def scatter_clone_fractions(clone_fractions_pivoted, plotsdir):
    
    num_samples = len(clone_fractions_pivoted.index.get_level_values(0).unique())
    fig, axes = plt.subplots(num_samples, 1, figsize=(8, 6 * num_samples))

    for i, sample in enumerate(clone_fractions_pivoted.index.get_level_values(0).unique()):
        if sample in ['s_10', 's_4', 's_8', 's_9']:
            try:
                sliced_pivoted_clones = clone_fractions_pivoted.loc[(sample, ['baseline', 'V1']), :]
                
                # Add scatter plot to each subplot
                sns.scatterplot(data=sliced_pivoted_clones,
                                x=sliced_pivoted_clones.loc[sample, 'baseline'],
                                y=sliced_pivoted_clones.loc[sample, 'V1'],
                                ax=axes[i])
                axes[i].set_title(sample)
                axes[i].set_xlabel('baseline')
                axes[i].set_ylabel('V1')
                axes[i].set_xscale('log')
                axes[i].set_yscale('log')
                # add a diagonal line
                axes[i].plot([min(sliced_pivoted_clones.loc[sample, 'baseline']), 1],
                             [min(sliced_pivoted_clones.loc[sample, 'V1']), 1], color='black')
                
                # change font of sns plot
                for item in ([axes[i].title, axes[i].xaxis.label, axes[i].yaxis.label] +
                        axes[i].get_xticklabels() + axes[i].get_yticklabels()):
                    item.set_fontsize(12)
                    item.set_fontname(FONTSNS)
                
            except KeyError:
                pass
            
            
    plt.tight_layout()
    plt.savefig(plotsdir + "clonefractionscatter.png", dpi=300)
    plt.close()
    
    
def alpha_diversity_tcr(abundances, pairs, res_dir):
    
    patients, condition = [el[0] for el in list(abundances.index)], [el[1] for el in list(abundances.index)]
    
    metrics = ['Shannon_index']
    shannon = alpha_diversity('shannon', abundances, patients)

    stat = pd.DataFrame({"Shannon_index": shannon,
                         "Condition": condition})
    
    # do a Mann-Whitney U test for each pair of conditions
    for pair in pairs:
        pos = pair[0]
        neg = pair[1]
        for metric in metrics:
            res = [mannwhitneyu(stat[stat['Condition'] == pos][metric],
                                stat[stat['Condition'] == neg][metric])[1]]
    
    fig, axes = plt.subplots(1, 1, figsize=(4, 4))
    
    for i, metric in enumerate(metrics):

        order = sorted(stat['Condition'].unique())

        sns.swarmplot(data=stat, y = metric, x = 'Condition',hue = 'Condition', linewidth=.5, size=5, order=order, ax=axes)
        sns.boxplot(data=stat, y = metric, x = 'Condition', hue = 'Condition', dodge= False, saturation = 0, width=0.4, showfliers=False, order=order, ax=axes)
        plt.tight_layout()
        
    # annotate the plot with the p-value
    axes.text(0.1, max(shannon)-0.01, 'MW-pvalue={:.5f}'.format(res[0]), ha='center', va='center')

    # change font of sns plot
    for item in ([axes.title, axes.xaxis.label, axes.yaxis.label] +
             axes.get_xticklabels() + axes.get_yticklabels()):
        item.set_fontsize(12)
        item.set_fontname(FONTSNS)
        
    plt.legend([],[], frameon=False)
    plt.title('Baseline Shannon entropies')
    fig.savefig(res_dir + 'alpha_diversity.png', dpi=600, bbox_inches='tight')
    plt.close()