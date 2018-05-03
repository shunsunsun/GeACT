# Credit: Wenting
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot
import umi_utils
import collections
import cluster
import plot_utils
import reference
import cached_pipeline
import star
import enrichr
import  mpl_toolkits.axes_grid1 as axes_grid
import numpy_utils as npu

import data.exp210 as u2os_data
import data.exp204 as hek_data

def references():
    return ["GRCh38.p7_gencode_comprehensive", "ERCC"]

def map(sample):
    return cached_pipeline.map(sample, star.get_index(references()), require_spacer_match = False, trim1 = 50, trim2 = 0)

def get_bams(samples):
    return collections.OrderedDict([(sample.name, map(sample)) for sample in samples])
    
def load_adjusted_counts(samples, references):
    """Returns UMI counts using either 3 or 7 max_distance depending on expression level"""
    
    bams    = get_bams(samples)
    counts3 = umi_utils.load_sample_counts(bams, references, max_distance = 3)
    counts7 = umi_utils.load_sample_counts(bams, references, max_distance = 7)
    return counts3.where(counts3 > 30, counts7)

def load_mean_adjusted_counts(samples, references):
    """Returns UMI counts using either 3 or 7 max_distance depending on mean expression level"""

    bams    = get_bams(samples)
    counts3 = umi_utils.load_sample_counts(get_bams(u2os_data.samples), references, max_distance = 3)
    counts7 = umi_utils.load_sample_counts(get_bams(u2os_data.samples), references, max_distance = 7)

    means3  = np.tile(counts3.mean(0)[None, :], [counts3.shape[0], 1])
    return counts7.where(means3 < 10, counts3)

def merge_diag(upper, lower, index = None, columns = None):
    if index is None:
        index = upper.index
    if columns is None:
        columns = upper.columns
    return pd.DataFrame(np.triu(upper.loc[index, columns], k = 0) + np.tril(lower.loc[index, columns], k = -1), index = index, columns = columns)


def format_expression_axes(data, axes = None, mappable = None, cax = None):
    if axes is None:
        axes = pyplot.gca()
    transposed = axes.get_xlim()[1] == data.shape[0]
    
    u2os_start = next(i for (i, x) in enumerate(data.index) if x.startswith("exp210"))

    set_yticks      = axes.set_xticks if transposed else axes.set_yticks
    set_yticklabels = axes.set_xticklabels if transposed else axes.set_yticklabels
    set_xticks      = axes.set_yticks if transposed else axes.set_xticks
    yaxis           = axes.xaxis if transposed else axes.yaxis
    
    set_yticks([u2os_start/2, u2os_start, u2os_start + (data.shape[0] - u2os_start)/2.])
    ticks = yaxis.get_major_ticks()
    ticks[0].tick1On = False
    ticks[2].tick1On = False
    set_yticklabels(["HEK293T", "", "U2OS"])
    set_xticks([])

    pyplot.colorbar(label = "Expression level (Z-Score)", ax = axes, orientation = "vertical" if transposed else "horizontal", mappable = mappable, cax = cax)

    
def plot_expr_corr(expression, correlation, plot_genes = None, ax = None, expr_kwargs = dict(), corr_kwargs = dict()):
    expr_kwargs.setdefault("cmap", pyplot.get_cmap("cool"))
    expr_kwargs.setdefault("aspect", "auto")
    expr_kwargs.setdefault("vmin", -3)
    expr_kwargs.setdefault("vmax", 3)

    corr_kwargs.setdefault("cmap", pyplot.get_cmap("bwr"))
    corr_kwargs.setdefault("aspect", "equal")
    corr_kwargs.setdefault("vmin", -.25)
    corr_kwargs.setdefault("vmax", .25)

    if plot_genes is None:
        plot_genes = correlation.columns
    expression  = expression.loc[:, plot_genes]
    correlation = correlation.loc[plot_genes, plot_genes]

    if ax is None:
        fig = pyplot.figure(figsize = (15, 10))
        ax2 = fig.add_subplot(1, 1, 1)
    else:
        ax2 = ax
        
    div  = axes_grid.make_axes_locatable(ax2)
    cax2 = div.append_axes("right", size = .5, pad = .2)
    cax1 = div.append_axes("left", size = .5, pad = 1)
    ax1  = div.append_axes("left", size = 3, pad = .2, sharey = ax2)

    im1  = plot_utils.heatmap(expression.T, axes = ax1, **expr_kwargs)
    im2  = plot_utils.heatmap(correlation, axes = ax2, **corr_kwargs)

    format_expression_axes(expression, mappable = im1, cax = cax1)
    pyplot.colorbar(mappable = im2, cax = cax2, label = "Correlation coefficient")
    ax2.set_xticks([])
    
    
    return (ax1, cax1, ax2, cax2)


def compare_counting_methods(references):
    def cluster_counts(counts):
        normed = umi_utils.normalize_counts(umi_utils.filter_counts(counts))
        return cluster.correlation_cluster(normed, abs_corr = True)
    
    bams    = get_bams(u2os_data.samples)

    counts3 = umi_utils.load_sample_counts(get_bams(u2os_data.samples), references, max_distance = 3)
    counts7 = umi_utils.load_sample_counts(get_bams(u2os_data.samples), references, max_distance = 7)

    countsA = load_adjusted_counts(u2os_data.samples)
    countsB = load_mean_adjusted_counts(u2os_data.samples)

    clustered3 = cluster_counts(counts3)
    clustered7 = cluster_counts(counts7)
    clusteredA = cluster_counts(countsA)
    clusteredB = cluster_counts(countsB)

    # plot correlation using different counting methods
    for clustered in [clustered3, clustered7, clusteredA, clusteredB]:
        pyplot.figure(figsize = (10, 10))
        pyplot.subplots_adjust(0, 0, 1, 1, wspace = 0, hspace = 0)
        plot_utils.heatmap(clustered, vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))

    # plot comparison of correlations
    pyplot.figure(figsize = (10, 10))
    pyplot.subplots_adjust(0, 0, 1, 1, wspace = 0, hspace = 0)
    plot_utils.heatmap(merge_diag(clusteredA, clustered3), vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))

    pyplot.figure(figsize = (10, 10))
    pyplot.subplots_adjust(0, 0, 1, 1, wspace = 0, hspace = 0)
    plot_utils.heatmap(merge_diag(clusteredA, clustered7), vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))

def clusters_overview(u2os_counts):
    normed    = umi_utils.normalize_counts(umi_utils.filter_counts(u2os_counts))
    clustered = cluster.correlation_cluster(normed, abs_corr = True)
    
    # identify cluster boundaries
    boundaries       = cluster.boundaries2(clustered, .8, .1)
    large_boundaries = [(a, b) for (a, b) in boundaries if b - a >= 10]


    # check clusters against Enrichr
    enrichments = []
    databases   = ["Reactome_2016", "WikiPathways_2016", "KEGG_2016", "ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]
    for (a, b) in large_boundaries:
        genes = clustered.columns[a:b]
        enrichments.append([enrichr.enrichment(genes, database, retries = 3) for database in databases])

    # These regions seem to be out of date
    # zoom_regions = {
    #     "misc":        [2550, 2950],
    #     "p53":         [7258, 7483],
    #     "bone_ecm":    [8923, 8949],
    #     "protein":     [9015, 9082],
    #     "cholesterol": [10959, 10972],
    # }

    zoom_regions = {
        "bone_ecm":    [7951, 7981],
    }
        
    heatmap_args = dict(vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))
    pyplot.figure(figsize = (10, 10))
    pyplot.subplots_adjust(0, 0, 1, 1, wspace = 0, hspace = 0)
    plot_utils.heatmap(clustered, **heatmap_args)
    plot_utils.outline_clusters(known_clusters.values())

    pyplot.figure()
    for (i, (name, (a, b))) in enumerate(zoom_regions.iteritems()):
        padding = max(int((b - a)*.2), 10)
        a -= padding
        b += padding
        pyplot.subplot(1, len(zoom_regions), i + 1)
        plot_utils.heatmap(clustered.iloc[a:b, a:b], **heatmap_args)
        pyplot.title("%s (%d-%d)" % (name, a, b))
        
        
    # plot correlations with numbered clusters
    heatmap_args = dict(vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))
    pyplot.figure(figsize = (10, 10))
    pyplot.subplots_adjust(0, 0, 1, 1, wspace = 0, hspace = 0)
    plot_utils.heatmap(clustered, **heatmap_args)
    plot_utils.outline_clusters(large_boundaries)
    for (i, (a, b)) in enumerate(large_boundaries):
        pyplot.gca().annotate(str(i), xy = ((a + b)/2., (a + b)/2.), xytext = (b + 1, a - 1))

    # print table of enrichment results for each cluster
    with open("temp/results/u2os_cluster_enrichment_condensed.txt", "w") as file:
        file.write("Cluster ID\tNumber of Genes\tGene Names\tEnrichment\n")
        for (i, (a, b)) in enumerate(large_boundaries):
            genes    = clustered.columns[a:b]
            enriched = [(db, results[0]) for (db, results) in zip(databases, enrichments[i]) if len(results) > 0 and results[0].q < .05 and len(results[0].genes) > 1]

            file.write("\t".join([
                str(i),
                str(len(genes)),
                ", ".join(genes),
                ", ".join(["%s (%s q=%f)" % (result.name, db, result.q) for (db, result) in enriched])
            ]))
            file.write("\n")

def compare_cell_types(hek_counts, u2os_counts):
    (u2os_filtered, hek_filtered) = umi_utils.joint_filter_counts([u2os_counts, hek_counts])
    (u2os_normed, hek_normed)     = umi_utils.joint_normalize_counts([u2os_filtered, hek_filtered])

    u2os_clustered  = cluster.correlation_cluster(u2os_normed, abs_corr = True)
    hek_clustered   = cluster.correlation_cluster(hek_normed, abs_corr = True)
    joint_clustered = cluster.joint_correlation_cluster([u2os_normed, hek_normed], abs_corr = True)

    merged_clustered = merge_diag(u2os_clustered, hek_clustered, index = joint_clustered.index, columns = joint_clustered.columns)

    joint_normed      = pd.concat([hek_normed, u2os_normed], 0)
    joint_transformed = (joint_normed - joint_normed.mean(0))/joint_normed.std(0)
    joint_ward        = cluster.ward_cluster(joint_transformed)
    # Reorder joint_ward so that all hek cells come before u2os cells
    # The cells are already completely separated by hek vs u2os, so this
    # just ensures the order is consitent with later correlation plots
    joint_ward = joint_ward.iloc[
        [i for (i, name) in enumerate(joint_ward.index) if name.startswith("exp204")] + [i for (i, name) in enumerate(joint_ward.index) if name.startswith("exp210")]
    ]

    de_clustered = cluster.correlation_cluster(pd.concat([hek_normed, u2os_normed], 0), abs_corr = False)

    expr_kwargs = dict(vmin = -1, vmax = 1, cmap = pyplot.get_cmap("cool"), aspect = "auto")
    corr_kwargs = dict(vmin = -.25, vmax = .25, cmap = pyplot.get_cmap("bwr"))

    # downsample cells for display
    # with npu.RandomSeedContext(0):
    #     hek_sampling  = np.random.choice(hek_filtered.index, size = 100, replace = False)
    #     u2os_sampling = np.random.choice(u2os_filtered.index, size = 100, replace = False)
    # downsampled_ward = joint_ward.loc[list(hek_sampling) + list(u2os_sampling), :]


    plot_expr_corr(joint_ward, de_clustered, expr_kwargs = expr_kwargs, corr_kwargs = corr_kwargs)
    
    # differential expression heatmaps
    pyplot.figure(figsize = (10, 3))
    plot_utils.heatmap(joint_ward, **expr_kwargs)
    format_expression_axes(joint_ward)
    pyplot.tight_layout()
    pyplot.savefig("temp/results/differential_expression_1.png")

    pyplot.figure(figsize = (3, 10))
    plot_utils.heatmap(joint_ward.T, **expr_kwargs)
    format_expression_axes(joint_ward)
    pyplot.tight_layout()
    pyplot.savefig("temp/results/differential_expression_2.png")

    # plot using correlation to order genes
    (ax1, cax1, ax2, cax2) = plot_expr_corr(joint_ward, merged_clustered, expr_kwargs = expr_kwargs, corr_kwargs = corr_kwargs)
    pyplot.savefig("temp/results/diff_corr.png")

    # zoom in
    region = range(9675, 10360)
    (ax1, cax1, ax2, cax2) = plot_expr_corr(joint_ward, merged_clustered.iloc[:, 9675:10360], expr_kwargs = expr_kwargs, corr_kwargs = corr_kwargs)
    pyplot.savefig("temp/results/diff_corr2.png")

    
    zoom_clusters = [
        # differential clustering, differential expression
        range(9211, 9221),                                                        # bone ECM
        #range(9247, 9253),
        range(9253, 9266),
        # common clustering, differential expression
        # range(9870, 9882),                                                      # not sure what this module does
        set(range(9310, 9336)) - set([9316, 9317, 9318, 9325, 9333, 9335]),
        # differential clustering, no differential expression5
        range(8143, 8173),
        # common clustering, no differential expression
        range(10222, 10252),
        range(10308, 10340),
    ]
    zoom_genes = [merged_clustered.columns[i] for r in zoom_clusters for i in r]
    for (i, indices) in enumerate(zoom_clusters):
        print "\nCluster %d:" % (i + 1)
        print "\n".join(merged_clustered.columns[list(indices)])
    plot_expr_corr(joint_ward, merged_clustered, plot_genes = zoom_genes, expr_kwargs = expr_kwargs, corr_kwargs = corr_kwargs)
    pyplot.savefig("temp/results/diff_corr3.png")

    
def main():
    u2os_counts = load_adjusted_counts(u2os_data.samples, references())
    hek_counts  = load_adjusted_counts(hek_data.hiseq_filtered_cell_samples, references())


