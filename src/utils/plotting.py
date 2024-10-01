# 2024 fabio-affaticati
import pandas as pd
import numpy as np


import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import seaborn as sns

from skbio.diversity import alpha_diversity
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from src.utils import tools

FONT = 'JetBrains Mono'
FONTSNS =  'DejaVu Sans'


def plot_nuniques_clonecounts(data, plotsdir):
    
    fig = px.scatter(data.query('TCR_Chain == "TRB"'), x='countUniqueCDR3', y='cloneCount',
                color='CONDITION', symbol='TIMEPOINTS', opacity=0.7,
                color_discrete_sequence=['red', 'blue'])

    fig.add_vline(x = 0, line_width=2, line_color="black")
    fig.update_layout(
        title='Number of unique CDR3s vs Summed Clone Counts',
        xaxis_title='Number of Unique CDR3s',
        yaxis_title='Summed Clone Counts',
        width=800,
        plot_bgcolor='white',
        font=dict(
            family=FONT,
            size=12,
            color='black'
        ),
        legend=dict(
            title='Condition',
            title_font_family=FONT,
            title_font_size=12,
            font=dict(
                family=FONTSNS,
                size=10,
                color='black'
            )
        )
    )

    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        showgrid=False
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        showgrid=False
    )
    fig.show()
    fig.write_image(plotsdir + 'summed_clonecounts_vs_unique_cdr3s.png', scale = 4)
    


def test_global_repertoire_differences(data, plotsdir):
    
    repertoire_comparison = data[['sample_id', 'SAMPLE', 'countUniqueCDR3', 'cloneCount', 'normalizedUniqueCDR3', 'CONDITION', 'TIMEPOINTS', 'TCR_Chain']]
    repertoire_comparison = repertoire_comparison.sort_values(by=['SAMPLE'])
    repertoire_comparison = repertoire_comparison.query('TCR_Chain == "TRB"')

    ref_dict = {}
    for condition in repertoire_comparison['CONDITION'].unique():
            ref_dict[condition] = {}
            for timepoint in repertoire_comparison['TIMEPOINTS'].unique():
                ref_dict[condition][timepoint] = repertoire_comparison.query('CONDITION == @condition and TIMEPOINTS == @timepoint')['normalizedUniqueCDR3'].values
                ref_dict[condition]['Mean_'+timepoint] = np.mean(repertoire_comparison.query('CONDITION == @condition and TIMEPOINTS == @timepoint')['normalizedUniqueCDR3'].values)
                ref_dict[condition]['SD_'+timepoint] = np.std(repertoire_comparison.query('CONDITION == @condition and TIMEPOINTS == @timepoint')['normalizedUniqueCDR3'].values)
                ref_dict[condition]['Labels_'+timepoint] = repertoire_comparison.query('CONDITION == @condition and TIMEPOINTS == @timepoint')['SAMPLE'].values
                
    ### Mann Whitney U test
    mw_tests = {}
    for timepoint in repertoire_comparison['TIMEPOINTS'].unique():
        stat, p_value = mannwhitneyu(ref_dict['Healthy'][timepoint], ref_dict['Lymphomas'][timepoint], alternative='two-sided')
        mw_tests[timepoint] = {'stat': stat, 'p_value': p_value}
        
    max_x = max([max(ref_dict[condition][timepoint]) for condition in repertoire_comparison['CONDITION'].unique() for timepoint in repertoire_comparison['TIMEPOINTS'].unique()])

    fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharex=True)
    for i, condition in enumerate(repertoire_comparison['CONDITION'].unique()):
        for j, timepoint in enumerate(repertoire_comparison['TIMEPOINTS'].unique()):
            ax[i][j].errorbar(x=ref_dict[condition][timepoint], 
                            y=ref_dict[condition]['Labels_'+timepoint],
                            xerr=ref_dict[condition]['SD_'+timepoint], 
                            fmt='o', color='k')
            
            ax[i][j].set_title(condition+' '+timepoint)
            ax[i][j].axvline(ref_dict[condition]['Mean_'+timepoint], ls='--')
            
            # get pvalue from mw_tests
            p_value = mw_tests[timepoint]['p_value']
            stat = mw_tests[timepoint]['stat']
            if p_value > 0.05:
                ax[i][j].set_facecolor('lightgrey')
            else:
                if stat > 0:
                    ax[0][j].set_facecolor('lightcoral')
                    ax[0][j].annotate('MW two-sided\npval = {:.4f}'.format(p_value), xy=(0.6, .8), xycoords='axes fraction', fontsize=8, color='black')
                else:
                    ax[1][j].set_facecolor('lightcoral')
                    ax[1][j].annotate('MW two-sided\npval = {:.4f}'.format(p_value), xy=(0.6, .8), xycoords='axes fraction', fontsize=8, color='black')
            
            if j == 0:
                ax[i][j].set_ylabel('Sample')
            if (j == 1) and (i == 1):
                ax[i][j].set_xlabel('Normalized Unique CDR3s')
            if i == 1:
                plt.setp(ax[1][j].get_xticklabels(), rotation=45)
            

    plt.xlim(0, max_x+0.05)
    fig.show()
    fig.savefig(plotsdir + 'errorbar_plot.png', dpi=600, bbox_inches='tight')   




def plot_repertoire_sizes(data, plotsdir):
    
    fig = make_subplots(rows=2, cols=2, shared_yaxes=True, subplot_titles=("Healthy TRA", "Lymphomas TRA", "Healthy TRB", "Lymphomas TRB"))

    timepoint_colors = {
        "baseline": '#004D40',
        "V1": '#1E88E5',
        "V3": '#FFC107'
    }

    healthy_data = data.query('CONDITION == "Healthy" and TCR_Chain == "TRA"')

    for timepoint in healthy_data["TIMEPOINTS"].unique():
        filtered_data = healthy_data[healthy_data["TIMEPOINTS"] == timepoint]
        fig.add_trace(
            go.Bar(x=filtered_data["SAMPLE"], y=filtered_data["normalizedUniqueCDR3"],
                name=timepoint, #text=filtered_data["normalizedUniqueCDR3"], 
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=True),
            row=1, col=1
        )
        
    lymphomas_data = data.query('CONDITION == "Lymphomas" and TCR_Chain == "TRA"')
    for timepoint in lymphomas_data["TIMEPOINTS"].unique():
        filtered_data = lymphomas_data[lymphomas_data["TIMEPOINTS"] == timepoint]
        fig.add_trace(
            go.Bar(x=filtered_data["SAMPLE"], y=filtered_data["normalizedUniqueCDR3"],
                name=timepoint, #text=filtered_data["normalizedUniqueCDR3"],
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=False),
            row=1, col=2
        )
        
    healthy_data = data.query('CONDITION == "Healthy" and TCR_Chain == "TRB"')
    for timepoint in healthy_data["TIMEPOINTS"].unique():
        filtered_data = healthy_data[healthy_data["TIMEPOINTS"] == timepoint]
        fig.add_trace(
            go.Bar(x=filtered_data["SAMPLE"], y=filtered_data["normalizedUniqueCDR3"],
                name=timepoint, #text=filtered_data["normalizedUniqueCDR3"], 
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=False),
            row=2, col=1
        )
        
    lymphomas_data = data.query('CONDITION == "Lymphomas" and TCR_Chain == "TRB"')
    for timepoint in lymphomas_data["TIMEPOINTS"].unique():
        filtered_data = lymphomas_data[lymphomas_data["TIMEPOINTS"] == timepoint]
        fig.add_trace(
            go.Bar(x=filtered_data["SAMPLE"], y=filtered_data["normalizedUniqueCDR3"],
                name=timepoint, #text=filtered_data["normalizedUniqueCDR3"],
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=False),
            row=2, col=2
        )
        
        
        
    # Update axis title for all subplots
    fig.update_xaxes(title_text="Sample", row=1, col=1)
    fig.update_xaxes(title_text="Sample", row=1, col=2)
    fig.update_xaxes(title_text="Sample", row=2, col=1)
    fig.update_xaxes(title_text="Sample", row=2, col=2)

    fig.update_yaxes(title_text="Normalised unique CDR3", row=1, col=1)
    fig.update_yaxes(title_text="Normalised unique CDR3", row=2, col=1)
        
    fig.update_traces(
        textfont_size=12, textangle=0, textposition="outside", cliponaxis=False, marker_line_width=2,
        
    )
    fig.update_layout(
        font=dict(family=FONT, size=12),
        #paper_bgcolor='rgb(243, 243, 243)',
        #plot_bgcolor='rgb(243, 243, 243)',
        barmode='group',
        bargap=0.1,
        title_text="Comparison of Unique CDR3 Ratios",
        xaxis_title="Sample",
        yaxis_title="Normalised unique CDR3",
        legend_title="Timepoints",
        height = 1000,
        width = 1800,
        template="plotly_white"
    )
    fig.update_yaxes(tickangle=45, tickmode='array', showgrid = False)
    fig.update_xaxes(tickangle=45, showgrid=True)
    fig.show()
    fig.write_image(plotsdir+"barplot.png", scale = 4)

    
    
    
    
    
def alpha_diversity_tcr(abundances, pairs, plotsdir, type_test):
    
    patients, condition = [el[0] for el in list(abundances.index)], [el[1] for el in list(abundances.index)]
    
    metrics = ['Normalised Shannon entropy']
    shannon = alpha_diversity('shannon', abundances, patients)

    stat = pd.DataFrame({"Shannon_index": shannon,
                         "Condition": condition})
    norm_factor = np.log2(len(abundances.columns))
    stat['Normalised Shannon entropy'] = stat['Shannon_index']/norm_factor
    
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

        # assign to the Conditions specific colors
        colors = {'Healthy': 'blue', 'Lymphomas': 'red'}
        sns.swarmplot(data=stat, y=metric, x='Condition', hue='Condition', linewidth=.5, size=5, palette=colors, order=order, ax=axes)
        sns.boxplot(data=stat, y = metric, x = 'Condition', hue = 'Condition', dodge= False, saturation = 0, width=0.4, showfliers=False, palette=colors, order=order, ax=axes)
        # select the hues of the boxplot

            
        plt.tight_layout()
        
    # annotate the plot with the p-value
    axes.text(.5, .9, 'MW-pvalue={:.3f}'.format(res[0]), ha='center', va='center')

    # change font of sns plot
    for item in ([axes.title, axes.xaxis.label, axes.yaxis.label] +
             axes.get_xticklabels() + axes.get_yticklabels()):
        item.set_fontsize(12)
        item.set_fontname(FONTSNS)
        
    plt.legend([],[], frameon=False)
    # set limits of y axis
    axes.set_ylim(0,1)
    plt.title(type_test + ' Shannon entropies')
    fig.savefig(plotsdir + 'alpha_diversity' + type_test + '.png', dpi=600, bbox_inches='tight')
    plt.close()





def plot_gene_entropy_comparison(stat, highlight_samples, ylimits, plotname, title, plotsdir):
    
    
    timepoint_colors = {
    "baseline": '#004D40',
    "V1": '#1E88E5',
    "V3": '#FFC107'
    }


    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, subplot_titles=("Healthy", "Lymphomas"))

    healthy_data = stat.query('Condition == "Healthy"')
    for timepoint in healthy_data["Timepoint"].unique():
        filtered_data = healthy_data[healthy_data["Timepoint"] == timepoint]
        fig.add_trace(
            go.Box(y=filtered_data["Shannon_index"],
                name=timepoint, text=filtered_data["sample_id"], boxpoints = 'all',
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=True),
            row=1, col=1
        )
        
        
    # Filter data to get the highlighted samples
    highlight_data = healthy_data[healthy_data['sample_id'].isin(highlight_samples)]

    # Add highlighted samples as a scatter plot trace
    fig.add_trace(
        go.Scatter(
            x=highlight_data['Timepoint'],
            y=highlight_data['Shannon_index'],
            mode='markers',
            marker=dict(color='red', size=10, symbol='star'),
            name='Already positive individuals',
            text=highlight_data['sample_id'],  # Hover text for highlighted samples
            showlegend=True),
        row=1, col=1
    )
        
    lymphoma_data = stat.query('Condition == "Lymphomas"')
    for timepoint in lymphoma_data["Timepoint"].unique():
        filtered_data = lymphoma_data[lymphoma_data["Timepoint"] == timepoint]
        fig.add_trace(
            go.Box(y=filtered_data["Shannon_index"],
                name=timepoint, text=filtered_data["sample_id"], boxpoints = 'all',
                marker=dict(color=timepoint_colors[timepoint], line=dict(width=0.5)),
                showlegend=False),
            row=1, col=2
        )


    # Filter data to get the highlighted samples
    highlight_data = lymphoma_data[lymphoma_data['sample_id'].isin(highlight_samples)]

    # Add highlighted samples as a scatter plot trace
    fig.add_trace(
        go.Scatter(
            x=highlight_data['Timepoint'],
            y=highlight_data['Shannon_index'],
            mode='markers',
            marker=dict(color='red', size=10, symbol='star'),
            name='Already positive individuals',
            text=highlight_data['sample_id'],  # Hover text for highlighted samples
            showlegend=False),
        row=1, col=2
    )

    fig.update_layout(
        font=dict(family=FONT, size=12),
        #paper_bgcolor='rgb(243, 243, 243)',
        #plot_bgcolor='rgb(243, 243, 243)',
        barmode='group',
        bargap=0.1,
        title_text=title,
        xaxis_title="Condition",
        yaxis_title="Normalised Shannon index",
        legend_title="Timepoints",
        height = 600,
        width = 1000,
        template="plotly_white"
    )


    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        showgrid=False
    )
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        showgrid=True,
        range= ylimits
    )
    fig.write_image(plotsdir + plotname + '.png', scale = 4)
    fig.show()
    
    
def testing_and_plot_singlegene_usage(gene_usage, variables_to_test, plotsdir):
    
    for variable in variables_to_test:
        print(f'Tested varialbe: {variable}')
        for timepoint in ['baseline', 'V1', 'V3']:
            pvalues = []
            v_genes = []
            filtered_data = gene_usage.query('TIMEPOINTS == @timepoint')
            for unique_v_call in filtered_data['v_call'].unique():
                #print(f'V gene: {unique_v_call}')
                res = mannwhitneyu(filtered_data.query('CONDITION == "Lymphomas" and v_call == @unique_v_call')[variable],
                                    filtered_data.query('CONDITION == "Healthy" and v_call == @unique_v_call')[variable])[1]
                pvalues.append(res)
                v_genes.append(unique_v_call)
            
            
            corrected_pvalues = multipletests(pvalues, method ='fdr_bh')[1]
            
            fig = go.Figure()
            
            for adj_pvalue, v_gene in zip(corrected_pvalues, v_genes):
                if adj_pvalue <= 0.05:
                    adj_pvalue
                    print(f'{v_gene} at {timepoint} has a significant adjusted p-value of {adj_pvalue}')
                    fig.add_trace(go.Box(
                        y=filtered_data.query('CONDITION == "Lymphomas" and v_call == @v_gene')[variable],
                        name=v_gene,
                        line_color='red',
                        showlegend=False,
                        width=0.3 
                    ))

                    # Create boxplot for Healthy condition
                    fig.add_trace(go.Box(
                        y=filtered_data.query('CONDITION == "Healthy" and v_call == @v_gene')[variable],
                        name=v_gene,
                        line_color='blue',
                        showlegend=False,
                        width=0.3 
                    ))
                    
                    # Add p-value annotation
                    offset = 0
                    
                    if variable == 'Normalised_Shannon_entropy':
                        offset = max(filtered_data.query('v_call == @v_gene')[variable]) + 0.1
                    elif variable == 'unique_counts':
                        offset = max(filtered_data.query('v_call == @v_gene')[variable]) + 4
                    elif variable == 'normalizedUniqueCounts':
                        offset = max(filtered_data.query('v_call == @v_gene')[variable]) + 0.000005
                        
                    fig.add_annotation(
                        x=v_gene,  # Place the annotation near the 'Healthy' group for this v_gene
                        # y position in the scale of the y-axis
                        y=offset,  # Position the annotation above the max y-value for the v_gene
                        # asterisk noation for significance
                        #text=f'adj pvalue = {adj_pvalue:.4f}',  # Format the p-value with three decimal places
                        text=f'{tools.convert_pvalue_to_asterisks(adj_pvalue)}',
                        showarrow=False,
                        font=dict(color="black", size=14),
                        xanchor='center'
                    )
                        
                        
                else:
                    pass
                    print(f'{v_gene} at {timepoint} has a non-significant adjusted p-value of {adj_pvalue}')
                    
            # if Figure empty, skip saving
            if len(fig.data) == 0:
                continue
            
            
            else:
                
                if variable == 'normalizedUniqueCounts':
                    
                    fig.update_yaxes(title_text="Normalised unique CDR3 counts", range=[0, 0.000025])
                    
                elif variable == 'unique_counts':
                    
                    fig.update_yaxes(title_text="Unique CDR3 counts",)
                
                elif variable == 'Normalised_Shannon_entropy':
                    
                    fig.update_yaxes(title_text="Normalised Shannon entropy")
                    
                # Update traces for all boxplots
                fig.update_traces(
                    jitter=0.05,  # add some jitter on points for better visibility
                    marker=dict(opacity=0.6),  # make markers slightly transparent
                    boxpoints='all',
                )
                
                # Update layout for the plot
                fig.update_layout(
                    title_text=f"Timepoint {timepoint} (only sig results reported ; MW BH-corrected pvalue <= 0.05)",
                    #yaxis=dict(range=[-1, 12]),  # Set the desired y-axis limits
                    xaxis_title="V gene",  # Custom x-axis title
                    #yaxis_title=f"{variable}",  # Custom y-axis title
                    height=500,
                    width=900,
                    font=dict(
                        family=FONT,
                    ),
                    boxmode='group',
                    template="plotly_white",
                    bargap=0.2
                )
            
                # Adding custom legend
                fig.add_trace(go.Scatter(
                    x=[None], y=[None],
                    mode='markers',
                    marker=dict(color='blue'),
                    showlegend=True,
                    name='Healthy'
                ))
                fig.add_trace(go.Scatter(
                    x=[None], y=[None],
                    mode='markers',
                    marker=dict(color='red'),
                    showlegend=True,
                    name='Lymphomas'
                ))
                
                fig.update_xaxes(
                    mirror=True,
                    ticks='outside',
                    showline=True,
                    linecolor='black',
                    showgrid=False
                )
                fig.update_yaxes(
                    mirror=True,
                    ticks='outside',
                    showline=True,
                    linecolor='black',
                    showgrid=True
                )

                #fig.show()
                fig.write_image(f'{plotsdir}{variable}_{timepoint}_covid.png', width = 1200, scale = 4)