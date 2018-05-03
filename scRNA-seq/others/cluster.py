# Credit: Wenting
import numpy as np
import numpy_utils as npu
import pandas as pd
import scipy.cluster.hierarchy as hcluster
import scipy.spatial.distance as scidist

def masked_corr(data, num_outliers = 1):
    mask = npu.outlier_mask(data.values, num_outliers = num_outliers, axis = 1)
    return npu.masked_corrcoef(data.T, mask.T)

def correlation_cluster(data, abs_corr = False, mask_outliers = 0):
    corr          = np.nan_to_num(masked_corr(data, mask_outliers))
    np.fill_diagonal(corr, 1.)
    corr_dist     = scidist.squareform(corr, checks = False)
    corr_dist     = 1 - np.abs(corr_dist) if abs_corr else 1 - corr_dist
    gene_clusters = hcluster.linkage(corr_dist, method = "average")
    gene_leaves   = hcluster.leaves_list(gene_clusters)

    clustered_corr  = corr[gene_leaves, :][:, gene_leaves]
    clustered_genes = [data.columns[i] for i in gene_leaves]
    rv              = pd.DataFrame(clustered_corr, columns = clustered_genes, index = clustered_genes)
    
    # transform the linkage matrix for the reordered columns
    new_indices          = np.argsort(gene_leaves).astype(np.int)
    new_linkage          = gene_clusters.copy()
    old0                 = np.nonzero(new_linkage[:, 0] < len(gene_leaves))[0]
    old1                 = np.nonzero(new_linkage[:, 1] < len(gene_leaves))[0]
    new_linkage[old0, 0] = new_indices[new_linkage[old0, 0].astype(np.int)]
    new_linkage[old1, 1] = new_indices[new_linkage[old1, 1].astype(np.int)]
    rv.linkage           = new_linkage
    
    return rv

def joint_correlation_cluster(data_list, abs_corr = False, mask_outliers = 0):
    corrs     = [np.nan_to_num(masked_corr(data, mask_outliers)) for data in data_list]
    for corr in corrs:
        np.fill_diagonal(corr, 1.)
    corr_dists = []
    for corr in corrs:
        corr_dist = scidist.squareform(corr, checks = False)
        corr_dist = 1 - np.abs(corr_dist) if abs_corr else 1 - corr_dist
        corr_dists.append(corr_dist)
    corr_dist = np.min(corr_dists, 0)
    gene_clusters = hcluster.linkage(corr_dist, method = "average")
    gene_leaves   = hcluster.leaves_list(gene_clusters)

    clustered_corr  = corr[gene_leaves, :][:, gene_leaves]
    clustered_genes = [data_list[0].columns[i] for i in gene_leaves]

    return pd.DataFrame(clustered_corr, columns = clustered_genes, index = clustered_genes)

def ward_cluster(data):
    linkage1 = hcluster.ward(data)
    linkage2 = hcluster.ward(data.T)
    return data.iloc[hcluster.leaves_list(linkage1), hcluster.leaves_list(linkage2)]

def boundaries(linkage, correlation_threshold = .1):
    cluster_ids = hcluster.fcluster(linkage, 1 - correlation_threshold, criterion = "distance")
    endpoints   = np.nonzero(np.ediff1d(cluster_ids, to_begin = 1, to_end = 1))[0]
    return zip(endpoints[:-1], endpoints[1:])

def linkage_subset(linkage, leaf_ids):
    num_orig_leaves = len(linkage) + 1
    id_map          = {leaf: i for (i, leaf) in enumerate(leaf_ids)}
    new_linkage     = np.zeros((len(leaf_ids) - 1, 4), dtype = linkage.dtype)

    next_id = len(leaf_ids)
    for (i, old_link) in enumerate(linkage):
        id       = i + num_orig_leaves
        left_id  = id_map.get(old_link[0], -1)
        right_id = id_map.get(old_link[1], -1)
        if left_id == -1 or right_id == -1:
            id_map[id] = max(left_id, right_id)
        else:
            new_link    = new_linkage[next_id - len(leaf_ids), :]
            new_link[0] = left_id
            new_link[1] = right_id
            new_link[2] = old_link[2]
            
            left_count  = new_linkage[left_id - len(leaf_ids), 3] if left_id >= len(leaf_ids) else 1
            right_count = new_linkage[right_id - len(leaf_ids), 3] if right_id >= len(leaf_ids) else 1
            new_link[3] = left_count + right_count
            id_map[id]  = next_id
            next_id    += 1
    return new_linkage

def boundaries2(data, percentile, threshold):

    def find_clusters(tree):
        if tree.is_leaf():
            return ([tree.id, tree.id + 1], False)
        
        (left_clusters, left_frozen)   = find_clusters(tree.left)
        (right_clusters, right_frozen) = find_clusters(tree.right)

        if left_frozen or right_frozen:
            return (left_clusters[:-1] + right_clusters, True)
        else:
            joint_values = np.sort(data.values[left_clusters[0]:left_clusters[-1], right_clusters[0]:right_clusters[-1]].flatten())
            if joint_values[int(len(joint_values)*percentile)] >= threshold:
                return ([left_clusters[0], right_clusters[-1]], False)
            else:
                return (left_clusters[:-1] + right_clusters, False)

    tree = hcluster.to_tree(data.linkage)
    (clusters, frozen) = find_clusters(tree)

    return zip(clusters[:-1], clusters[1:])

def coclusters(clustered1, boundaries1, clustered2, differences = False, size_threshold = 10, corr_threshold = .5):
    common = []
    for (a, b) in boundaries1:
        if b - a < size_threshold:
            continue
        genes  = clustered1.columns[a:b]
        mean1  = np.mean(np.abs(np.nan_to_num(scidist.squareform(clustered1.loc[genes, genes], checks = False))))
        mean2  = np.mean(np.abs(np.nan_to_num(scidist.squareform(clustered2.loc[genes, genes], checks = False))))
        if mean2 > mean1*corr_threshold:
            if not differences:
                common.append((a, b))
        elif differences:
            common.append((a, b))

    return common

def average_clusters(data, boundaries1, boundaries2, scale_clusters = True):
    filtered = np.array(data)[[i for (a, b) in boundaries1 for i in range(a, b)], :][:, [i for (a, b) in boundaries2 for i in range(a, b)]]

    cluster_sizes1 = [b - a for (a, b) in boundaries1]
    cluster_sizes2 = [b - a for (a, b) in boundaries2]
    split_points1  = np.cumsum(cluster_sizes1[:-1])
    split_points2  = np.cumsum(cluster_sizes2[:-1])
    averaged       = np.array([np.mean(x, 1) for x in np.hsplit(filtered, split_points2)]).T
    averaged       = np.array([np.mean(x, 0) for x in np.vsplit(averaged, split_points1)])

    if scale_clusters:
        scaled = np.repeat(averaged, cluster_sizes1, 0)
        scaled = np.repeat(scaled, cluster_sizes2, 1)
        if isinstance(data, pd.DataFrame):
            return pd.DataFrame(
                scaled,
                index   = [data.index[i] for (a, b) in boundaries1 for i in range(a, b)],
                columns = [data.columns[i] for (a, b) in boundaries2 for i in range(a, b)],
            )
        else:
            return scaled
    else:
        return pd.DataFrame(averaged) if isinstance(data, pd.DataFrame) else averaged

def random_cluster_averages(data, cluster_sizes, num_samplings):
    data            = np.array(data)
    boundary_points = [0] + list(np.cumsum(cluster_sizes))
    boundaries      = zip(boundary_points[:-1], boundary_points[1:])
    averages        = []
    for i in range(num_samplings):
        permutation = np.random.permutation(len(data))
        averages.append(average_clusters(data[permutation, :][:, permutation], boundaries, boundaries, False))
    return np.array(averages)
    
