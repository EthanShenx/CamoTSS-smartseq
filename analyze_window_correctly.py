#!/usr/bin/env python
"""
Correct analysis of the sliding window algorithm
"""

import numpy as np
import pandas as pd

def analyze_sliding_window_correctly():
    print("CORRECT ANALYSIS OF SLIDING WINDOW ALGORITHM")
    print("="*60)
    
    # Simulate the real scenario: we have reads at specific positions
    # Each read has position, cell barcode, cigar string
    reads_data = [
        # High count at position 975 (with unencoded G)
        (975, 'cell_A', '14S49M'),
        (975, 'cell_B', '14S49M'),
        (975, 'cell_C', '14S49M'),
        (975, 'cell_D', '14S49M'),
        (975, 'cell_E', '14S49M'),
        (975, 'cell_F', '14S49M'),
        (975, 'cell_G', '14S49M'),
        (975, 'cell_H', '14S49M'),
        (975, 'cell_I', '14S49M'),
        (975, 'cell_J', '14S49M'),
        (975, 'cell_K', '14S49M'),
        (975, 'cell_L', '14S49M'),
        # More reads at 975
        (975, 'cell_M', '14S49M'),
        (975, 'cell_N', '14S49M'),
        (975, 'cell_O', '14S49M'),
        
        # Some reads at nearby positions (lower count)
        (976, 'cell_P', '50M'),
        (977, 'cell_Q', '50M'),
        (978, 'cell_R', '50M'),
        
        # High count at position 995 (with unencoded G)
        (995, 'cell_S', '15S48M'),
        (995, 'cell_T', '15S48M'),
        (995, 'cell_U', '15S48M'),
        (995, 'cell_V', '15S48M'),
        (995, 'cell_W', '15S48M'),
        (995, 'cell_X', '15S48M'),
        (995, 'cell_Y', '15S48M'),
        (995, 'cell_Z', '15S48M'),
    ]
    
    print(f"Simulated {len(reads_data)} reads")
    print(f"Reads with unencoded G (14S/15S/16S): {[r for r in reads_data if '14S' in r[2] or '15S' in r[2] or '16S' in r[2]]}")
    
    # Filter for unencoded G as CamoTSS does
    filterls = [i for i in reads_data if ('14S' in i[2]) or ('15S' in i[2]) or ('16S' in i[2])]
    print(f"After unencoded G filter: {len(filterls)} reads")
    
    # Calculate the TSS position and corresponding counts
    promoterTSS = [read[0] for read in filterls]
    TSS, count = np.unique(promoterTSS, return_counts=True)
    
    print(f"Unique positions with unencoded G: {TSS}")
    print(f"Counts at these positions: {count}")
    
    # Create desired dataframe with all positions in range
    TSS_start, TSS_end = 970, 1000
    originaldf = pd.DataFrame({'count': count}, index=TSS)
    desireddf = pd.DataFrame({'count': 0}, index=[i for i in range(int(TSS_start), int(TSS_end))])
    desireddf['count'] = desireddf.index.map(originaldf['count']).fillna(0)
    desireddf.reset_index(inplace=True)
    nonzeroarray = np.array(desireddf)
    
    print(f"All positions in range [{TSS_start}, {TSS_end}): {len(desireddf)}")
    print(f"Non-zero counts: {desireddf[desireddf['count'] > 0]}")
    
    # Sort by position for positive strand
    strand = '+'
    if strand == '+':
        sortfinalarray = nonzeroarray[nonzeroarray[:, 0].argsort()]
    elif strand == '-':
        sortfinalarray = nonzeroarray[nonzeroarray[:, 0].argsort()[::-1]]

    TSS_full = sortfinalarray.T[0]
    count_full = sortfinalarray.T[1]
    
    print(f"Sorted positions: {TSS_full[TSS_full > 0][:10]}...")  # First 10
    print(f"Sorted counts: {count_full[count_full > 0][:10]}...")  # First 10
    
    # Now apply sliding window
    windowSize = 15
    leftIndex = 0
    minCTSSCount = 10
    minFC = 2.0
    
    print(f"\nApplying sliding window (size={windowSize}, leftIndex={leftIndex}):")
    
    storels = []
    for i in range(len(TSS_full) - windowSize + 1):
        onewindow = TSS_full[i: i + windowSize]
        correspondingcount = count_full[i: i + windowSize]
        
        # The key insight: it takes the leftmost position [leftIndex=0], not the middle!
        middlecount = correspondingcount[leftIndex]
        avgcount = sum(correspondingcount) / len(correspondingcount)
        foldchange = (middlecount + 1) / (avgcount + 1)
        
        if middlecount > 0:  # Only show windows with non-zero counts
            print(f"  Window starting at pos {int(onewindow[leftIndex])}: "
                  f"pos range [{int(onewindow[0])}, {int(onewindow[-1])}], "
                  f"selected count={int(middlecount)}, avg={avgcount:.2f}, FC={foldchange:.3f}")
        
        storels.append([onewindow[leftIndex], middlecount, foldchange])

    # Sort by fold change descending
    foldchangels = [i[2] for i in storels]
    sortindex = sorted(range(len(foldchangels)), key=lambda k: foldchangels[k], reverse=True)
    allsortls = [storels[i] for i in sortindex]
    
    print(f"\nTop 5 candidates by fold-change:")
    for i in range(min(5, len(allsortls))):
        pos, cnt, fc = allsortls[i]
        if cnt > 0:  # Only show non-zero
            print(f"  Pos {int(pos)}: count={int(cnt)}, FC={fc:.3f}")
    
    # Apply filters
    keepCTSS = [ele for ele in allsortls if (ele[1] > minCTSSCount) and (ele[2] > minFC)]
    
    print(f"\nApplying filters: count > {minCTSSCount} AND FC > {minFC}")
    print(f"CTSS passing filters: {len(keepCTSS)}")
    
    if keepCTSS:
        print("Passed CTSS:", [(int(item[0]), int(item[1]), item[2]) for item in keepCTSS])
    else:
        print("No CTSS passed filters.")
        print("\nREASON FOR FAILURE:")
        print("1. The algorithm uses 'leftIndex=0', meaning it evaluates the LEFTMOST position in each window")
        print("2. For a position to be evaluated, it must be at the start of a window")
        print("3. Our high-count positions (975, 995) might not align with window starts")
        print("4. The fold-change is calculated as (leftmost_count + 1) / (window_avg + 1)")
        print("5. Even with high counts, if the window contains many zero-count positions, FC will be low")

def demonstrate_working_example():
    print("\n" + "="*60)
    print("DEMONSTRATING A WORKING EXAMPLE")
    print("="*60)
    
    # Create a scenario where we know it should work
    # Have high count at position that will be the start of a window
    import pandas as pd
    
    # Simulate a region where position 1000 has high count
    TSS_start, TSS_end = 1000, 1020
    all_positions = list(range(int(TSS_start), int(TSS_end)))
    
    # Create counts: position 1000 has high count, others have low
    all_counts = [0] * len(all_positions)
    all_counts[0] = 40  # Position 1000 has 40 reads (high)
    all_counts[5] = 5   # Position 1005 has 5 reads (low)
    
    # Create the dataframe as CamoTSS does
    desireddf = pd.DataFrame({'count': all_counts}, index=all_positions)
    nonzeroarray = np.array(desireddf.reset_index())
    
    # Sort (for positive strand, it's already sorted)
    TSS_full = nonzeroarray.T[0]  # positions
    count_full = nonzeroarray.T[1]  # counts
    
    print(f"Positions: {TSS_full}")
    print(f"Counts: {count_full}")
    
    # Apply sliding window
    windowSize = 10
    leftIndex = 0
    minCTSSCount = 10
    minFC = 2.0
    
    storels = []
    for i in range(len(TSS_full) - windowSize + 1):
        onewindow = TSS_full[i: i + windowSize]
        correspondingcount = count_full[i: i + windowSize]
        
        middlecount = correspondingcount[leftIndex]
        avgcount = sum(correspondingcount) / len(correspondingcount)
        foldchange = (middlecount + 1) / (avgcount + 1)
        
        print(f"Window starting at {int(onewindow[leftIndex])}: count={int(middlecount)}, avg={avgcount:.2f}, FC={foldchange:.3f}")
        
        storels.append([onewindow[leftIndex], middlecount, foldchange])

    # Apply filters
    keepCTSS = [ele for ele in storels if (ele[1] > minCTSSCount) and (ele[2] > minFC)]
    print(f"\nCTSS passing filters: {len(keepCTSS)}")
    if keepCTSS:
        print("Passed:", [(int(item[0]), int(item[1]), item[2]) for item in keepCTSS])

if __name__ == "__main__":
    analyze_sliding_window_correctly()
    demonstrate_working_example()