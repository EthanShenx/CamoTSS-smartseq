#!/usr/bin/env python
"""
Deep dive into CamoTSS algorithm to understand the core functionality.
This script tests individual components of the algorithm.
"""

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.optimize import linear_sum_assignment
import statistics
import editdistance
import pickle
import os
from pathlib import Path

# Test the clustering algorithm
def test_clustering():
    print("Testing clustering algorithm...")
    
    # Simulate read positions (TSS positions)
    np.random.seed(42)
    positions = np.concatenate([
        np.random.normal(1000, 10, 60),  # Cluster 1: ~60 reads around pos 1000
        np.random.normal(1300, 10, 40),  # Cluster 2: ~40 reads around pos 1300
        np.random.normal(2000, 10, 30),  # Cluster 3: ~30 reads around pos 2000
    ])
    
    # Add cell barcodes for each position
    cell_barcodes = ['cell_' + str(i % 50) for i in range(len(positions))]
    
    # Add cigar strings (simulated)
    cigar_strings = ['14S49M' if i < 100 else '15S48M' for i in range(len(positions))]
    
    # Combine into tuples as CamoTSS does
    read_info = [(pos, cb, cs) for pos, cb, cs in zip(positions, cell_barcodes, cigar_strings)]
    
    print(f"Generated {len(read_info)} simulated reads")
    
    # Apply CamoTSS clustering algorithm
    InnerDistance = 100  # Default parameter
    minCount = 50  # Default parameter
    
    # Prepare data for clustering
    posiarray = np.array([t[0] for t in read_info]).reshape(-1,1)
    CBarray = np.array([t[1] for t in read_info]).reshape(-1,1)
    cigartuplearray = np.array([t[2] for t in read_info]).reshape(-1,1)
    
    # Perform hierarchical clustering
    clusterModel = AgglomerativeClustering(
        n_clusters=None,
        linkage='average',
        distance_threshold=InnerDistance
    )
    
    clusterModel = clusterModel.fit(posiarray)
    labels = clusterModel.labels_
    
    # Count reads per cluster
    label, count = np.unique(labels, return_counts=True)
    selectlabel = label[count >= minCount]  # Only clusters with sufficient reads
    selectcount = count[count >= minCount]
    
    print(f"Found {len(selectlabel)} clusters with at least {minCount} reads:")
    for i, (lbl, cnt) in enumerate(zip(selectlabel, selectcount)):
        cluster_positions = posiarray[labels == lbl].flatten()
        print(f"  Cluster {i}: {cnt} reads, positions {cluster_positions.min():.0f}-{cluster_positions.max():.0f}")
    
    # Prepare clustered results
    altTSSls = []
    for i in range(0, len(selectlabel)):
        mask = labels == selectlabel[i]
        altTSSls.append([posiarray[mask], CBarray[mask], cigartuplearray[mask]])
    
    print(f"Final clusters after filtering: {len(altTSSls)}")
    return altTSSls

# Test the feature calculation for ML filtering
def test_feature_calculation():
    print("\nTesting feature calculation for ML filtering...")
    
    # Simulate a cluster result from clustering
    np.random.seed(42)
    cluster_pos = np.random.normal(1000, 10, 60)  # 60 reads in this cluster
    cluster_cb = ['cell_' + str(i % 30) for i in range(60)]  # 30 unique cells
    cluster_cigar = ['14S49M' if i < 30 else '15S48M' for i in range(60)]  # Some have unencoded G
    
    # Calculate the 4 features used by CamoTSS
    count = len(cluster_pos)  # Total UMI count
    std = statistics.stdev(cluster_pos)  # Standard deviation of positions
    summit_count = np.max(np.unique(cluster_pos, return_counts=True)[1])  # Max reads at single position
    unencoded_G_percent = sum([('14S' in ele) or ('15S' in ele) or ('16S' in ele) 
                              for ele in cluster_cigar]) / count  # Fraction with unencoded G
    
    print(f"Cluster features:")
    print(f"  UMI_count: {count}")
    print(f"  SD: {std:.2f}")
    print(f"  summit_UMI_count: {summit_count}")
    print(f"  unencoded_G_percent: {unencoded_G_percent:.2f}")
    
    # These 4 features are fed into the ML model for filtering
    features = [count, std, summit_count, unencoded_G_percent]
    print(f"  Feature vector: {features}")
    
    return features

# Test the Hungarian algorithm for TSS annotation
def test_hungarian_assignment():
    print("\nTesting Hungarian algorithm for TSS annotation...")
    
    # Simulate known TSS positions from reference
    ref_tss_positions = [995, 1005, 1302, 1308]  # Known TSS from reference
    gene_id = "ENSG00000123456"
    
    # Simulate detected clusters (from clustering)
    detected_clusters = [
        {'positions': np.random.normal(1000, 5, 55)},  # Should match ref pos 995 or 1005
        {'positions': np.random.normal(1305, 5, 45)}   # Should match ref pos 1302 or 1308
    ]
    
    print(f"Reference TSS positions: {ref_tss_positions}")
    print(f"Detected clusters: {len(detected_clusters)}")
    
    # Create cost matrix (distance between detected cluster and reference TSS)
    cost_mtx = np.zeros((len(detected_clusters), len(ref_tss_positions)))
    for i in range(len(detected_clusters)):
        for j in range(len(ref_tss_positions)):
            # Calculate cost as distance to modal position in cluster
            cluster_vals = detected_clusters[i]['positions']
            position, count = np.unique(cluster_vals, return_counts=True)
            mode_position = position[np.argmax(count)]
            cost_mtx[i, j] = abs(mode_position - ref_tss_positions[j])
    
    print(f"Cost matrix:\n{cost_mtx}")
    
    # Apply Hungarian algorithm
    row_ind, col_ind = linear_sum_assignment(cost_mtx)
    print(f"Optimal assignment - detected[{row_ind}] -> ref[{col_ind}]")
    
    # Map to transcript IDs
    ref_transcripts = [f"ENST00000{i:06d}" for i in range(len(ref_tss_positions))]
    assigned_transcripts = [ref_transcripts[j] for j in col_ind]
    print(f"Assigned transcripts: {assigned_transcripts}")
    
    return assigned_transcripts

# Test the sliding window algorithm for CTSS detection
def test_sliding_window():
    print("\nTesting sliding window algorithm for CTSS detection...")
    
    # Simulate reads within a TSS cluster
    np.random.seed(42)
    # Simulate positions with some having higher counts (potential CTSS)
    positions = list(range(950, 1050))  # Region from 950 to 1050
    # Create some peaks in the counts to simulate real data
    counts = [5] * len(positions)
    counts[25] = 45  # Peak at position 975
    counts[45] = 38  # Peak at position 995
    counts[65] = 52  # Peak at position 1015
    
    # Add unencoded G filter simulation
    cigars = ['14S49M' if i in [25, 45, 65] else '50M' for i in range(len(positions))]
    
    # Filter for reads with unencoded G (14S, 15S, 16S)
    filtered_data = [(pos, cnt, cigar) for pos, cnt, cigar in 
                     zip(positions, counts, cigars) 
                     if ('14S' in cigar) or ('15S' in cigar) or ('16S' in cigar)]
    
    print(f"Filtered reads with unencoded G: {len(filtered_data)}")
    
    # Apply sliding window
    window_size = 15  # Default parameter
    minCTSSCount = 10  # Minimum count threshold
    minFC = 2.0  # Minimum fold change threshold
    
    results = []
    for i in range(len(filtered_data) - window_size + 1):
        window_data = filtered_data[i:i + window_size]
        if len(window_data) < window_size:
            continue
            
        # Calculate fold change for middle position vs average of window
        middle_idx = 0  # In CamoTSS implementation, it's leftIndex=0
        middle_count = window_data[middle_idx][1]
        avg_count = sum(item[1] for item in window_data) / len(window_data)
        fold_change = (middle_count + 1) / (avg_count + 1)  # Add pseudocount
        
        results.append({
            'position': window_data[middle_idx][0],
            'count': middle_count,
            'fold_change': fold_change,
            'avg_window': avg_count
        })
    
    # Sort by fold change (descending)
    results.sort(key=lambda x: x['fold_change'], reverse=True)
    
    print(f"Top 5 CTSS candidates (sliding window results):")
    for i, result in enumerate(results[:5]):
        print(f"  {i+1}. Pos {result['position']}: count={result['count']:.0f}, "
              f"FC={result['fold_change']:.2f}, avg={result['avg_window']:.1f}")
    
    # Apply filters
    filtered_results = [r for r in results 
                       if r['count'] >= minCTSSCount and r['fold_change'] >= minFC]
    
    print(f"\nAfter applying filters (minCount≥{minCTSSCount}, minFC≥{minFC}): "
          f"{len(filtered_results)} CTSS retained")
    
    return filtered_results

def main():
    print("Deep Dive Analysis of CamoTSS Algorithm Components\n")
    print("="*60)
    
    # Test clustering
    clusters = test_clustering()
    
    # Test feature calculation
    features = test_feature_calculation()
    
    # Test Hungarian assignment
    transcripts = test_hungarian_assignment()
    
    # Test sliding window for CTSS
    ctss_candidates = test_sliding_window()
    
    print("\n" + "="*60)
    print("SUMMARY OF ALGORITHM COMPONENTS:")
    print("1. Clustering: Hierarchical clustering with distance threshold")
    print("2. Filtering: ML model using 4 features (count, SD, summit, unencoded-G%)")
    print("3. Annotation: Hungarian algorithm to match clusters to known TSS")
    print("4. CTSS Detection: Sliding window with fold-change calculation")
    
    print("\nThe algorithm is designed to:")
    print("- Detect TSS clusters from scRNA-seq 5' data")
    print("- Filter false positives using ML model")
    print("- Annotate clusters to known transcripts or identify novel TSS")
    print("- Detect cleavage sites (CTSS) within TSS clusters")
    print("- Output cell-by-TSS count matrices in AnnData format")

if __name__ == "__main__":
    main()