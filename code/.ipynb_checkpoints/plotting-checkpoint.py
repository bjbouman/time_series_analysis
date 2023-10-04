### import libraries
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap #for own cmap
import matplotlib.colors as plt_colors # for converting rgb to hex
from matplotlib import colors
from matplotlib import cm
import seaborn as sns
plt.rcParams['pdf.fonttype'] = 42 #for saving PDF with changeable text
plt.rcParams['ps.fonttype'] = 42 #for saving PDF with changeable text
from utils import *

### functions
def plot_genesets(adata, geneset, save_path=None, save=False, plot="UMAP"):
    '''
    Plot a selected set of gene markers in a UMAP or dotplot format.
    
    Parameters
    ----------
    adata
        Annotated data matrix
    geneset
        Dataframe with in each column marker genes assigned to one
        specific cell type. The column names should include the 
        corresponding cell type names.
    key: `str`
        String for how to name the plot when saving
    plot: `str` (default: UMAP)
        String determining whether to return a set of UMAPs, a 
        dotplot or both.
    save: `bool` (default:False)
        Save the plot to folder figures.
    '''
    ### calculate scores for the genes in the set
    genes_set_scores = calculate_geneset_scores(adata, geneset)
    
    if plot=="UMAP" or plot=="both": 

        # plot scores for selected genes set
        no_of_rows = int(len(genes_set_scores.columns)/5)+1
        fig, axs = plt.subplots(no_of_rows, 5, 
                                figsize=(20,4*no_of_rows), 
                                gridspec_kw={'wspace':0.3, 'hspace':0.3})
        axs = axs.ravel()

        for i in range(len(genes_set_scores.columns)):

            column = genes_set_scores.columns[i]

            # add score for plotting to the adata object for plotting
            adata.obs['marker_gene_score'] = genes_set_scores[column]

            # UMAP per cell type in the set
            sc.pl.umap(adata, color='marker_gene_score', title=column,
                       frameon=False, ncols=5, show=False, ax=axs[i], size=30)

        for ax in axs:
            ax.axis('off')
            if ax.get_legend(): ax.legend().set_visible(False)

        # save figure
        if save==True:
            save_path = "../figures/1.2.3.UMAPs_celltype_marker_"+key+".png"
            fig.savefig(save_path, bbox_inches='tight', format='png', dpi=300)

        del adata.obs['marker_gene_score']
    
    
    if plot=="dotplot" or plot=="both":
        
        # add score for plotting to the adata object for plotting
        adata.obs[genes_set_scores.columns] = genes_set_scores

        # plot scores for selected genes set
        fig, axs = plt.subplots(1, 1, figsize=(genes_set_scores.shape[1]*0.4,len(adata.obs["clusters"].cat.categories)*0.4),
                                gridspec_kw={'wspace':0.3, 'hspace':0.3})

        sc.pl.dotplot(adata, var_names=genes_set_scores.columns, 
                      groupby="clusters", standard_scale="var", 
                      show=False, ax=axs)

        for column in genes_set_scores.columns:
            del adata.obs[column] 

        # save figure
        if save==True:
            fig.savefig(save_path, bbox_inches='tight', format='png', dpi=300)

def plot_UMAP_per_timepoint(
    adata, 
    variable, 
    min_value, 
    max_value, 
    var_label=None, 
    save=False, 
    time_key="time", 
    save_path=None
):
    '''
    Plot the expression of a selected variable (e.g. a gene or a gene set score) in
    UMAPs for each timepoint
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    variable: `str` 
        Variable in adata object to be plotted.
    time_key: `str` (default:"time")
        The name of the column in adata.obs that you want to use for
        splitting the data.    
    '''
    # create own cmap
    colors_cmap = ["lightgrey", "cornsilk", "orange",
              "red", "purple", "midnightblue"]
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors_cmap)
    
    # create plot
    timepoints = adata.obs[time_key].cat.categories

    fig, axs = plt.subplots(1, len(timepoints), figsize=(2.25*len(timepoints),2), gridspec_kw={'wspace':0, 'hspace':0.1})

    axs = axs.ravel()

    for idx, timepoint in enumerate(timepoints):
        sc.pl.umap(adata[adata.obs[time_key]==timepoint], 
                   color=variable, ax=axs[idx], show=False, size=15, 
                   title="", vmin=min_value, vmax=max_value, 
                   cmap=custom_cmap, add_outline=True)

        plt.gcf().axes[-1].remove()

    for ax in axs.flat:
        ax.set_axis_off()

    axs[0].get_yaxis().set_visible(True)  

    cb_ax = fig.add_axes([0.72, 0.0, 0.16, 0.04])
    norm = colors.Normalize(vmin=min_value, vmax=max_value)
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), orientation="horizontal", cax=cb_ax)

    # add timepoint label to each UMAP
    for i, timepoint in enumerate(timepoints):
        plt.text(0.4, 0.1, timepoint,ha='center',va='center',transform = axs[i].transAxes)

    if var_label==None:
        fig.text(0.72, 0.08, variable, va='center')   
    else:
        fig.text(0.72, 0.08, var_label, va='center')


    ### save figure
    if save==True:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
    
    plt.show()            
        
def plot_discrete_expression(
    adata, 
    gene, 
    split_key="clusters", 
    time_key="time", 
    save=False, 
    save_path=None,
):
    '''
    Plot the expression of the selected gene with experimental time 
    on the x-axis and expression on the y-axis. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    gene: `str` 
        Gene(s) to be plotted.
    split_key: `str` (default:"clusters")
        The name of the column in adata.obs that you want to use to
        split the data in the plot.
    time_key: `str` (default:"time")
        The name of the column in adata.obs that you want to use for
        splitting the data. 
    '''
    
    # create plot
    fig, axs = plt.subplots(1, 1,figsize = (7,2.5))
    
    # get colors matching the clusters
    labels = adata.obs[split_key].cat.categories
    colors = adata.uns[split_key+"_colors"]
    color_dict = dict(zip(labels,colors))

    # calculate size of each cluster
    cluster_sizes = {}
    for cluster in adata.obs[split_key].unique():
        cluster_sizes[cluster] = adata[adata.obs[split_key]==cluster].obs[time_key].value_counts().values

    for cluster in adata.obs[split_key].cat.categories:
        timepoints = adata.obs[time_key].cat.categories
        std_per_timepoint = []
        mean_per_timepoint = []

        for time in timepoints:
            std = np.std([adata[(adata.obs[split_key]==cluster)&(adata.obs[time_key]==time)&(adata.obs["hashtags"]==i)][:,gene].X.mean() for i in ('tag2', 'tag3', 'tag4', 'tag1')])
            std_per_timepoint.append(std)
            mean = adata[(adata.obs[split_key]==cluster)&(adata.obs[time_key]==time)][:,gene].X.mean()
            mean_per_timepoint.append(mean)

        axs.errorbar(timepoints, mean_per_timepoint, yerr=std_per_timepoint,
                     fmt='o', c=color_dict[cluster], capsize=6)
        axs.plot(timepoints, mean_per_timepoint, c=color_dict[cluster])

    axs.set_title(gene)
    axs.grid(False)
    axs.yaxis.tick_right()
        
    #remove spines
    axs.spines['left'].set_visible(False)
    axs.spines['top'].set_visible(False)
            
    fig.text(1, 0.5, 'gene expression', va='center', rotation='vertical')
    axs.set(xlabel='experimental time')
    
    ### save figure
    if save==True:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
    
    plt.show()
        
### function to create heatmap with response genes
def plot_heatmap_scores(
    adata,
    scores, 
    groups, 
    row_linkage, 
    name="", 
    save=False, 
    save_path=None
):

    # # create own cmap
    # colors = ["deepskyblue","lightskyblue","aliceblue","white",
    #           "mistyrose", "lightsalmon", "orangered"]
    # custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    
    # add row colors so ax_row_colors is created
    row_colors = ["#ffffff"]*len(groups)

    # plot heatmap
    fig = sns.clustermap(
        scores.reset_index(drop=True), 
        row_linkage=row_linkage, 
        col_cluster=False,
        method="ward", 
        figsize=(15, 30), 
        cmap='Reds',
        yticklabels=False,
        xticklabels=False, 
        row_colors=row_colors,
        cbar_kws={"orientation": "horizontal"})

    # change position colorbar
    fig.cax.set_position([0.06, .81, .15, .01])

    # remove row dendogram
    fig.ax_row_dendrogram.set_visible(False)

    # add cluster labels 
    fig.ax_col_dendrogram.scatter(adata.obs["clusters"].cat.categories,[0]*len(adata.obs["clusters"].cat.categories), color=adata.uns['clusters_colors'], s=600)
    fig.ax_col_dendrogram.set_position([0.23, 0.8, 0.76, 0.03])
    for i in range(len(adata.obs["clusters"].cat.categories)):
        fig.ax_col_dendrogram.text(x=i-0.15,
                                   y=0.1, rotation=90,
                                   s=adata.obs["clusters"].cat.categories[i],
                                   fontdict=dict(color="black",size=20))

    # add black lines
    group_sizes = groups.value_counts().sort_index()
    size=0
    for i in group_sizes.index[0:len(group_sizes)-1]:
        size+=group_sizes[i]
        fig.ax_heatmap.axhline(size,color='black',ls='-', xmin=-0.18, xmax=1, clip_on = False, linewidth=3)

    # add group numbers and ISG numbers for each group
    net_ind = (group_sizes.cumsum()-(group_sizes/2))
    net_names = ["group " + str(i) for i in range(1, len(group_sizes)+1)]
    #net_names[7] = ""
    net_names[9] = ""
    fig.ax_row_colors.set_yticks(net_ind)
    fig.ax_row_colors.set_yticklabels(net_names, fontdict=dict(color="black",size=20))
    fig.ax_row_colors.yaxis.set_tick_params(size=0)

    # save figure
    if save==True:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
    
    plt.show()
        
def plot_in_pseudotime(
    adata,
    items,
    clusters=None,
    colors=None,
    labels=None,
    plot_scores=False,
    smoothing_sampling_rounds=20,
    smoothing_sample_size=0.5,
    smoothing_window_size=300,
    plot_all_clusters=False,
    plot_all_clusters_separate=False,
    plot_legend=True,
    legend_location="right",
    add_CI=True,
    scale=False,
    ymin=None,
    ymax=None,
    save=False,
    save_path=None,
):
    '''
    Plot the expression of the selected gene(s) in pseudotime.   
    '''
    
    # preset set of clusters if specified
    if plot_all_clusters == True:
        clusters = [[
            'HSCs #1', 'HSCs #2', 'LMPPs #1', 'LMPPs #2', 'myel. prog. #1',
            'myel. prog. #2', 'myel. prog. #3', 'ery. prog. #1', 'ery. prog. #2',
            'ery. prog. #3', 'MK prog.', 'eosinophil prog.',
        ]]
    
    elif plot_all_clusters_separate == True:
        clusters = [
            ['HSCs #1'], ['HSCs #2'], ['LMPPs #1'], ['LMPPs #2'], ['myel. prog. #1'],
            ['myel. prog. #2'], ['myel. prog. #3'], ['ery. prog. #1'], ['ery. prog. #2'],
            ['ery. prog. #3'], ['MK prog.'], ['eosinophil prog.'],
        ]

    # abort function if given colors are more/less than no. of clusters to be plotted
    if colors is not None:
        if len(colors) != len(clusters):
            sys.exit("Number of colors do not match the number of clusters.")
            
    # abort function if given labels are more/less than no. of clusters to be plotted
    if labels is not None:
        if len(labels) != len(clusters):
            sys.exit("Number of labels do not match the number of clusters.")
            
    # abort if multiple items and multiple clusters are specified
    if (len(items) > 1 & len(clusters) > 1):
        sys.exit("Enter either multiple items or multiple clusters. Never both.")
    
    # get colors matching the cell types
    cell_type_labels = adata.obs["clusters"].cat.categories
    cell_type_colors = adata.uns["clusters_colors"]
    color_dict = dict(zip(cell_type_labels, cell_type_colors))
    
    # set colors used for plotting
    if colors is None:
        if len(items) > 1:
            colors = [plt.cm.Pastel2(idx) for idx in range(len(items))]
        elif len(clusters) > 1:
            colors = [color_dict[cluster[0]] for cluster in clusters]
        elif len(items) == 1:
            colors = [plt.cm.Pastel2(0)]
            
    # set labels used for legend
    if labels is None:
        if len(clusters) > 1:
            labels = [', '.join(cluster) for cluster in clusters]
        elif len(items) >= 1:
            labels = items
    
    # intialize figure
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))

    # intialize first color and first label
    color_idx = 0
    label_idx = 0
    
    # loop over items or clusters to be plotted in figure
    for idx, item in enumerate(items):

        for cluster in clusters:

            # select clusters
            mask = adata.obs['clusters'].isin(cluster)

            # select x (time) and y (expression) data
            x = adata[mask].obs["pt_ordering"].values
            
            
            # select y (expression or score) data
            if item in adata.var_names.values:
                y = adata[:,item][mask].X.flatten()
                smoothing_method = "bin"
                if(scale==True):
                    ylabel = "scaled gene expression"
                else:
                    ylabel = "gene expression"
                
            elif item in adata.obs.columns.values:
                y = adata.obs[item][mask].values
                smoothing_method = "loess"
                if(scale==True):
                    ylabel = "scaled gene set score"
                else:
                    ylabel = "gene set score"
            
            # apply smoothing 
            xgrid, ygrids, mean, stderr = bootstrap_smoothing(
                x, y,
                smoothing_method,
                smoothing_sampling_rounds,
                smoothing_sample_size,
                smoothing_window_size,
            )
            
            upper_margin = mean+1.96*stderr
            lower_margin = mean-1.96*stderr

            # scale data
            if(scale==True):
                lower_margin = (lower_margin-np.nanmin(mean))/(np.nanmax(mean)-np.nanmin(mean))
                upper_margin = (upper_margin-np.nanmin(mean))/(np.nanmax(mean)-np.nanmin(mean))
                mean = (mean-np.nanmin(mean))/(np.nanmax(mean)-np.nanmin(mean))

            # plot mean and error margin
            if add_CI==True:
                ax.fill_between(xgrid, lower_margin, upper_margin, alpha=0.1, color=colors[color_idx])
            ax.plot(xgrid, mean, color=colors[color_idx], label=labels[label_idx])
            
            # set new color and new label for next plotted line
            color_idx += 1
            label_idx += 1
    
    # add legend
    if plot_legend == True:        
        if legend_location == 'right':
            plt.legend(loc="lower left", bbox_to_anchor=(1, 0), ncol=1, frameon=False)
        elif legend_location == 'top':
            plt.legend(loc="lower left", bbox_to_anchor=(0, 1), ncol=3, frameon=False)

    # add labels axes
    plt.ylabel(ylabel)
    plt.xlabel("pseudotime")
    
    # remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # change y axis
    ax.set_ylim(ymin, ymax)
    
    # save figure
    if (save==True) & (save_path is not None):
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
   
    plt.show()
    
def plot_groups_in_pattern(
    adata, 
    patterns,
    all_response_genes, 
    groups_scores,
    text="percentage",
    save=False,
    save_path=None,
):
    '''
    Plot the distribution of selected patterns. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    all_response_genes: XXX  
    patterns: `str` 
        Pattern(s) to be plotted.
    text: `str`
        Whether to write the percentage ("percentage") or the no. of genes ("genes").
    save: `bool` (default: False)
        Whether to save the figure or not.
    '''
    
    ### create plot
    fig, axs = plt.subplots(len(patterns), 2,
                            figsize = (22,1.5*len(patterns)),
                            gridspec_kw={'wspace':0.01, 'hspace':2,
                                         'width_ratios': [1, 12]})

    # create plot with selected patterns
    for idx, pattern in enumerate(patterns):
        
        # select response genes in pattern 
        response_genes = all_response_genes[all_response_genes.pattern==("pattern"+pattern)].index

        # get expression of response genes
        expressions = adata[:,response_genes].X.todense()
        expressions = pd.DataFrame(expressions, columns=response_genes, index=adata.obs["pseudotime"])
        expressions = expressions.sort_index()
        expressions = expressions.T

        # take mean expression per 100 cells
        expressions = expressions.groupby(np.arange(len(expressions.columns))//100, axis=1).mean()

        # scale expression between 0 and 1
        expressions = expressions.sub(expressions.min(axis=1), axis=0).divide((expressions.max(axis=1) - expressions.min(axis=1)), axis=0)

        # create plot with pattern
        axs[idx,0].plot(expressions.mean(axis=0).rolling(10).sum(), 
                        linewidth=2, color="black")
        axs[idx,0].spines['top'].set_visible(False)
        axs[idx,0].spines['right'].set_visible(False)
        axs[idx,0].spines['left'].set_visible(False)
        axs[idx,0].spines['bottom'].set_visible(False)
        axs[idx,0].get_xaxis().set_ticks([])
        axs[idx,0].get_yaxis().set_ticks([])
        axs[idx,0].set_ylabel("pattern "+pattern)
        
        ###################################################
        ### create bar plot
        groups_in_pattern = all_response_genes[all_response_genes.pattern==("pattern"+pattern)].group.value_counts()
        groups_percentages = groups_in_pattern/groups_in_pattern.sum()*100
        #groups_percentages = groups_percentages.sort_index()

        total = 0

        # add label for each percentage
        for percentage in groups_percentages:
            axs[idx,1].barh(("pattern"+pattern), percentage, left=total, 
                            color=plt.cm.Pastel2(0), edgecolor='black')
            total += percentage

        # add percentages as text to barplot
        for i, rect in enumerate(axs[idx,1].patches):

            # get coordinates percentage label
            width = rect.get_width()
            x = rect.get_x()
            label_x = x + width / 2

            # write label (percentage or no. of genes)
            if text=="pattern":
                percentage_text = f'{width:.1f}'+'%'
                min_width = 4.0
            elif text == "genes":
                percentage_text = f'{groups_in_pattern.iloc[i]}'+" gene(s)"
                min_width = 5.0

            # add percentage label (only when height is greater than value)
            if width > min_width:
                axs[idx,1].text(label_x, 0, percentage_text, ha='center', va='center')
         
        # remove spines and ticks
        axs[idx,1].spines['right'].set_visible(False)
        axs[idx,1].spines['left'].set_visible(False)
        axs[idx,1].spines['top'].set_visible(False)
        axs[idx,1].spines['bottom'].set_visible(False)
        axs[idx,1].set_xticks([])
        axs[idx,1].set_yticks([])
        
        ###################################################
        ### add umap with specificity
        umap_x = 0
        
        for i, rect in enumerate(axs[idx,1].patches):
            
            width = rect.get_width()
            x = rect.get_x()
            umap_x += width/2
            
            # add only when height is greater than value
            if width > 6.0:
                
                axin1 = axs[idx,1].inset_axes([(umap_x-2.5), 0.5, 5, 1], 
                                              transform=axs[idx,1].transData)
                # score_colors = np.array([plt_colors.rgb2hex(plt.cm.Reds(score)) for score in groups_scores[str(groups_percentages.index[i])]])
                # adata.uns['clusters_colors'] = score_colors
                # sc.pl.umap(adata, color="clusters", ax=axin1, show=False, size=3, title="")
                # axin1.get_legend().remove()
                # axin1.spines['right'].set_visible(False)
                # axin1.spines['top'].set_visible(False)
                # axin1.spines['left'].set_visible(False)
                # axin1.spines['bottom'].set_visible(False)
                # axin1.set_xlabel("")
                # axin1.set_ylabel("")
                axin1.set_axis_off()
                
                axin1.text(
                    0, 0.4, 
                    groups_percentages.index[i],
                    horizontalalignment='center', 
                    verticalalignment='center', 
                    transform=axin1.transAxes,
                    bbox=dict(boxstyle="circle",pad=0.3, fc="lightgrey", alpha=0.5, lw=2))
            
            umap_x += width/2

    fig.patch.set_facecolor('xkcd:white')
    
    ### save figure
    if save==True:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
    
    plt.show()
    
def plot_patterns_in_group(adata, 
    groups,
    all_response_genes, 
    groups_scores,
    text="pattern",
    save=False,
    save_path=None,
):
    '''
    Plot the distribution of  of selected pattern groups. 
    
    Parameters
    ----------
    adata
        Annotated data matrix.
    groups: `str` 
        Cell type specificity group(s) to be plotted.
    text: `str`
        Whether to write the percentage ("percentage") or the no. of genes ("genes").
    save: `bool` (default: False)
        Whether to save the figure or not.
    '''
    
    ### create plot
    fig, axs = plt.subplots(len(groups), 2,
                            figsize = (22,1*len(groups)),
                            gridspec_kw={'wspace':0.01, 'hspace':1.3,
                                         'width_ratios': [1, 25]})
    
    for idx, group in enumerate(groups):
    
        ###################################################
        ### create plot with pattern of group
        
        score_colors = np.array([plt_colors.rgb2hex(plt.cm.Reds(score)) for score in groups_scores[group]])        
        adata.uns['cell_types_colors'] = score_colors
        #sc.pl.umap(adata, color="cell_types", ax=axs[idx,0], show=False, size=3, title="")

        #axs[idx,0].get_legend().remove()
        axs[idx,0].spines['right'].set_visible(False)
        axs[idx,0].spines['top'].set_visible(False)
        axs[idx,0].spines['left'].set_visible(False)
        axs[idx,0].spines['bottom'].set_visible(False)
        axs[idx,0].set_xticks([])
        axs[idx,0].set_yticks([])
        
        axs[idx,0].set_ylabel("group "+group)
        axs[idx,0].set_xlabel("")
    
        ###################################################
        ### create bar plot
        patterns_in_group = all_response_genes[all_response_genes.group==int(group)].pattern.value_counts()
        patterns_percentages = patterns_in_group/patterns_in_group.sum()*100
        #patterns_percentages = patterns_percentages.sort_index()

        total = 0

        for percentage in patterns_percentages:
            axs[idx,1].barh(("group"+group), percentage, left=total, 
                            color=plt.cm.Pastel2(2), edgecolor='black')
            total += percentage

        # add percentages as text to barplot
        for i, rect in enumerate(axs[idx,1].patches):
           
            # get coordinates percentage label
            width = rect.get_width()
            x = rect.get_x()
            label_x = x + width / 2
            
            # write label (percentage or no. of genes)
            if text=="pattern":
                percentage_text = f'{width:.1f}'+'%'
                min_width = 4.0
            elif text == "genes":
                percentage_text = f'{patterns_in_group.iloc[i]}'+" gene(s)"
                min_width = 5.0

            # add only when height is greater than value
            if width > min_width:
                axs[idx,1].text(label_x, 0, percentage_text, ha='center', va='center')
                
        # remove spines and ticks
        axs[idx,1].spines['right'].set_visible(False)
        axs[idx,1].spines['left'].set_visible(False)
        axs[idx,1].spines['top'].set_visible(False)
        axs[idx,1].spines['bottom'].set_visible(False)
        axs[idx,1].set_xticks([])
        axs[idx,1].set_yticks([])
        
        ###################################################
        ### add patterns above barplot
        umap_x = 0
        
        for i, rect in enumerate(axs[idx,1].patches):
            
            width = rect.get_width()
            x = rect.get_x()
            umap_x += width/2
            
            # add only when height is greater than value
            if width > 6.0:
                
                axin1 = axs[idx,1].inset_axes([(umap_x-4), 0.5, 8, 0.5], 
                                              transform=axs[idx,1].transData)
                # select response genes in pattern 
                response_genes = all_response_genes[all_response_genes.pattern==("pattern"+patterns_percentages.index[i].replace("pattern",""))].index

                # get expression of response genes
                expressions = adata[:,response_genes].X.todense()
                expressions = pd.DataFrame(expressions, columns=response_genes, index=adata.obs["pseudotime"])
                expressions = expressions.sort_index()
                expressions = expressions.T

                # take mean expression per 100 cells
                expressions = expressions.groupby(np.arange(len(expressions.columns))//100, axis=1).mean()

                # scale expression between 0 and 1
                expressions = expressions.sub(expressions.min(axis=1), axis=0).divide((expressions.max(axis=1) - expressions.min(axis=1)), axis=0)

                # create plot with pattern
                axin1.plot(expressions.mean(axis=0).rolling(10).sum(), 
                                linewidth=2, color="black")
                axin1.spines['top'].set_visible(False)
                axin1.spines['right'].set_visible(False)
                axin1.spines['left'].set_visible(False)
                axin1.spines['bottom'].set_visible(False)
                axin1.get_xaxis().set_ticks([])
                axin1.get_yaxis().set_ticks([])
                
                
                axin1.text(0.5, 1.2, patterns_percentages.index[i].replace("pattern",""), 
                           horizontalalignment='center', verticalalignment='center', transform=axin1.transAxes,
                           bbox=dict(boxstyle="circle",pad=0.3, fc="lightgrey", alpha=0.5, lw=2))
            
            umap_x += width/2

    fig.patch.set_facecolor('xkcd:white')
    
    ### save figure
    if save==True:
        fig.savefig(save_path, bbox_inches='tight', format='pdf', dpi=300)
    
    plt.show()