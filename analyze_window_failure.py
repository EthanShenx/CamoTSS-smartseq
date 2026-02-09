#!/usr/bin/env python
"""
Analysis of why the sliding window test failed
"""

import numpy as np

def analyze_sliding_window_failure():
    print("ANALYZING SLIDING WINDOW FAILURE")
    print("="*50)
    
    # Recreate the scenario from the failed test
    np.random.seed(42)
    positions = list(range(950, 1051))  # 101 positions from 950 to 1050
    counts = [3] * len(positions)  # Base level of 3
    
    # Add peaks at specific positions
    counts[25] = 45  # Position 975: high count
    counts[45] = 38  # Position 995: high count  
    counts[65] = 52  # Position 1015: high count
    
    # Simulate cigar strings (some with unencoded G)
    cigars = ['50M'] * len(positions)
    cigars[25] = '14S49M'  # Position 975 has unencoded G
    cigars[45] = '15S48M'  # Position 995 has unencoded G
    cigars[65] = '16S47M'  # Position 1015 has unencoded G
    
    # Filter for unencoded G as CamoTSS does
    genereads = [(pos, f'cell_{i}', cigar) 
                 for i, (pos, count, cigar) in enumerate(zip(positions, counts, cigars))
                 if ('14S' in cigar) or ('15S' in cigar) or ('16S' in cigar)]
    
    print(f"Total positions: {len(positions)}")
    print(f"Positions with unencoded G: {len(genereads)}")
    print(f"Unencoded G positions: {[r[0] for r in genereads]}")
    print(f"Their counts: {[counts[pos-950] for pos in [r[0] for r in genereads]]}")
    
    TSS_start, TSS_end = 950, 1050
    strand = '+'
    windowSize = 15
    minCTSSCount = 10
    minFC = 2.0
    
    # Calculate TSS positions and counts within region
    promoterTSS = []
    for read in genereads:
        tss = read[0]
        if (tss >= TSS_start) and (tss <= TSS_end):
            promoterTSS.append(tss)
    
    TSS, count = np.unique(promoterTSS, return_counts=True)
    nonzeroarray = np.asarray((TSS, count)).T
    
    print(f"TSS positions after filtering: {TSS}")
    print(f"Counts at these positions: {count}")
    
    # Sort by position for positive strand
    if strand == '+':
        sortfinalarray = nonzeroarray[nonzeroarray[:, 0].argsort()]
        TSS = sortfinalarray.T[0]
        count = sortfinalarray.T[1]
    
    print(f"Sorted TSS positions: {TSS}")
    print(f"Sorted counts: {count}")
    
    print(f"\nApplying sliding window (size={windowSize}):")
    
    # Apply sliding window exactly as CamoTSS does
    storels = []
    leftIndex = 0  # This is crucial - it's always 0 in CamoTSS
    
    for i in range(len(TSS) - windowSize + 1):
        onewindow = TSS[i: i + windowSize]
        correspondingcount = count[i: i + windowSize]
        middlecount = correspondingcount[leftIndex]  # This is the FIRST position in window, not the middle!
        foldchange = (middlecount + 1) / (sum(correspondingcount)/len(correspondingcount) + 1)
        storels.append([onewindow[leftIndex], correspondingcount[leftIndex], foldchange])
        
        print(f"  Window {i}: positions {onewindow}, counts {correspondingcount}, "
              f"selected pos {onewindow[leftIndex]}, count {middlecount}, FC {foldchange:.3f}")
    
    # Sort by fold change descending
    foldchangels = [i[2] for i in storels]
    sortindex = sorted(range(len(foldchangels)), key=lambda k: foldchangels[k], reverse=True)
    allsortls = [storels[i] for i in sortindex]
    
    print(f"\nTop 5 candidates by fold-change:")
    for i in range(min(5, len(allsortls))):
        pos, cnt, fc = allsortls[i]
        print(f"  Pos {int(pos)}: count={cnt}, FC={fc:.3f}")
    
    # Apply filters
    keepCTSS = [ele for ele in allsortls if (ele[1] > minCTSSCount) and (ele[2] > minFC)]
    
    print(f"\nApplying filters: count > {minCTSSCount} AND FC > {minFC}")
    print(f"CTSS passing filters: {len(keepCTSS)}")
    
    if keepCTSS:
        print("Passed CTSS:", [(int(item[0]), item[1], item[2]) for item in keepCTSS])
    else:
        print("No CTSS passed filters. Let's analyze why...")
        
        # Analyze the issue
        print(f"\nAnalyzing the problem:")
        print(f"- The algorithm selects the LEFTMOST position in each window (leftIndex=0)")
        print(f"- So for a window containing [975, 976, 977, ...], it evaluates position 975")
        print(f"- Even if 975 has high count, other positions in the window affect fold-change")
        print(f"- For peak at 975, it needs to be in a window where it's the leftmost position")
        print(f"- The fold-change formula: (middlecount + 1) / (avg_in_window + 1)")
        
        # Show what happens with the peak at 975 (index 0 in our filtered TSS)
        peak_idx = 0  # Position 975 is at index 0
        if peak_idx < len(TSS) - windowSize + 1:
            window_start = TSS[peak_idx: peak_idx + windowSize]
            window_counts = count[peak_idx: peak_idx + windowSize]
            peak_count = window_counts[0]  # leftmost position
            avg_count = sum(window_counts) / len(window_counts)
            fc = (peak_count + 1) / (avg_count + 1)
            
            print(f"\nFor peak at 975 (count={peak_count}):")
            print(f"  Window: {window_start}")
            print(f"  Counts: {window_counts}")
            print(f"  Average: {avg_count:.2f}")
            print(f"  Fold-change: {fc:.3f}")
            print(f"  Passes count filter (>10)? {peak_count > minCTSSCount}")
            print(f"  Passes FC filter (>2.0)? {fc > minFC}")

def demonstrate_corrected_approach():
    print("\n" + "="*50)
    print("DEMONSTRATING CORRECTED APPROACH")
    print("="*50)
    
    # Create a better test case where peak is at the start of a window
    positions = list(range(950, 970))  # Smaller region for clarity
    counts = [2] * len(positions)
    counts[0] = 50  # Position 950 has high count
    counts[5] = 45  # Position 955 has high count
    
    cigars = ['50M'] * len(positions)
    cigars[0] = '14S49M'  # Position 950 has unencoded G
    cigars[5] = '15S48M'  # Position 955 has unencoded G
    
    genereads = [(pos, f'cell_{i}', cigar) 
                 for i, (pos, count, cigar) in enumerate(zip(positions, counts, cigars))
                 if ('14S' in cigar) or ('15S' in cigar) or ('16S' in cigar)]
    
    print(f"Test positions with unencoded G: {[r[0] for r in genereads]}")
    print(f"Their counts: {[counts[pos-950] for pos in [r[0] for r in genereads]]}")
    
    TSS, count = np.unique([r[0] for r in genereads], return_counts=True)
    print(f"Sorted TSS: {TSS}, counts: {count}")
    
    windowSize = 6  # Smaller window for this example
    leftIndex = 0
    
    print(f"\nSliding window analysis (windowSize={windowSize}):")
    for i in range(len(TSS) - windowSize + 1):
        onewindow = TSS[i: i + windowSize]
        correspondingcount = count[i: i + windowSize]
        middlecount = correspondingcount[leftIndex]  # First position in window
        avgcount = sum(correspondingcount) / len(correspondingcount)
        foldchange = (middlecount + 1) / (avgcount + 1)
        
        print(f"  Window {i}: pos {onewindow}, counts {correspondingcount}, "
              f"selected {onewindow[leftIndex]} (count={middlecount}), FC={foldchange:.3f}")

if __name__ == "__main__":
    analyze_sliding_window_failure()
    demonstrate_corrected_approach()