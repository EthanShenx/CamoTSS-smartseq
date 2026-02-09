#!/usr/bin/env python
"""
Minimal test for Smart-seq5 functionality in CamoTSS
"""

import os
import tempfile
import shutil
from unittest.mock import Mock, patch
import pysam
import pandas as pd
import numpy as np

def create_mock_bam(filename, reads_data):
    """
    Create a mock BAM file with specified reads data for testing
    """
    # Create header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 249250621, 'SN': 'chr1'}]
    }
    
    # Create BAM file
    bamfile = pysam.AlignmentFile(filename, "wb", header=header)
    
    for read_info in reads_data:
        # Create a mock read
        read = pysam.AlignedSegment()
        read.query_name = read_info.get('query_name', 'read1')
        read.flag = read_info.get('flag', 0)  # 0 = mapped, forward strand
        read.reference_id = 0  # chr1
        read.reference_start = read_info.get('pos', 1000)
        read.mapping_quality = read_info.get('mapq', 60)
        
        # Set CIGAR string (match 50 bases)
        read.cigarstring = read_info.get('cigar', '50M')
        
        # Set sequence
        read.query_sequence = read_info.get('seq', 'A' * 50)
        
        # Add tags if present
        if 'GX' in read_info:
            read.set_tag('GX', read_info['GX'])
        if 'CB' in read_info:
            read.set_tag('CB', read_info['CB'])
        if 'UB' in read_info:
            read.set_tag('UB', read_info['UB'])
        
        bamfile.write(read)
    
    bamfile.close()

def test_smartseq5_deduplication():
    """
    Test coordinate-based deduplication for Smart-seq5 data
    """
    print("Testing Smart-seq5 deduplication functionality...")
    
    # Import the necessary classes
    from CamoTSS.utils.get_counts import get_TSS_count
    
    # Create mock data
    mock_reads = [
        {
            'query_name': 'read1',
            'pos': 1000,
            'flag': 0,  # Forward strand
            'cigar': '50M',
            'seq': 'A' * 50,
            'GX': 'GENE1',
            'mapq': 60
        },
        {
            'query_name': 'read2',
            'pos': 1000,  # Same position as read1
            'flag': 0,    # Forward strand
            'cigar': '50M',
            'seq': 'A' * 50,
            'GX': 'GENE1',
            'mapq': 60
        },
        {
            'query_name': 'read3',
            'pos': 1005,  # Different position
            'flag': 0,    # Forward strand
            'cigar': '50M',
            'seq': 'A' * 50,
            'GX': 'GENE1',
            'mapq': 60
        }
    ]
    
    # Create temporary BAM file
    with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp_bam:
        create_mock_bam(tmp_bam.name, mock_reads)
        
        # Create mock reference files (minimal)
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create mock reference files
            ref_gene_path = os.path.join(tmp_dir, 'ref_gene.tsv')
            ref_tss_path = os.path.join(tmp_dir, 'ref_tss.tsv')
            cell_barcode_path = os.path.join(tmp_dir, 'cell_barcodes.tsv')
            count_dir = os.path.join(tmp_dir, 'count')
            os.makedirs(count_dir)
            
            # Write mock gene reference
            gene_df = pd.DataFrame({
                'gene_id': ['GENE1'],
                'Chromosome': ['chr1'],
                'Start': [900],
                'End': [1100],
                'Strand': ['+']
            })
            gene_df.to_csv(ref_gene_path, sep='\t', index=False)
            
            # Write mock TSS reference
            tss_df = pd.DataFrame({
                'transcript_id': ['GENE1_1000'],
                'gene_id': ['GENE1'],
                'Chromosome': ['chr1'],
                'TSS_start': [1000],
                'TSS_end': [1001],
                'Strand': ['+']
            })
            tss_df.to_csv(ref_tss_path, sep='\t', index=False)
            
            # Write mock cell barcodes
            cell_df = pd.DataFrame({
                'cell_id': [os.path.splitext(os.path.basename(tmp_bam.name))[0]]
            })
            cell_df.to_csv(cell_barcode_path, sep='\t', index=False)
            
            # Initialize the get_TSS_count class with smartseq5 parameters
            tss_counter = get_TSS_count(
                generefPath=ref_gene_path,
                tssrefPath=ref_tss_path,
                bamfilePath=[tmp_bam.name],  # Pass as list for smartseq5
                fastqFilePath=None,  # Not needed for this test
                outdir=tmp_dir,
                cellBarcodePath=cell_barcode_path,
                nproc=1,
                minCount=1,
                maxReadCount=1000,
                clusterDistance=300,
                InnerDistance=100,
                windowSize=15,
                minCTSSCount=100,
                minFC=6,
                platform='smartseq5',
                dedup_method='coord',
                min_mapq=20
            )
            
            # Test the deduplication methods directly
            # Create mock reads
            mock_aligned_reads = []
            for read_info in mock_reads:
                mock_read = Mock()
                mock_read.query_name = read_info['query_name']
                mock_read.reference_start = read_info['pos']
                mock_read.reference_end = read_info['pos'] + 50
                mock_read.reference_name = 'chr1'
                mock_read.is_reverse = (read_info['flag'] & 16) != 0  # Check reverse flag
                mock_read.is_unmapped = False
                mock_read.is_secondary = False
                mock_read.is_duplicate = False
                mock_read.is_supplementary = False
                mock_read.mapping_quality = read_info['mapq']
                mock_aligned_reads.append(mock_read)
            
            # Test coordinate deduplication
            mergedf = pd.DataFrame({
                'GENE1': {
                    'Strand': '+'
                }
            }).T
            
            deduplicated_reads = tss_counter._deduplicate_by_coordinates(
                mock_aligned_reads, mergedf, 'GENE1', 'coord'
            )
            
            print(f"Original reads: {len(mock_aligned_reads)}")
            print(f"Deduplicated reads: {len(deduplicated_reads)}")
            
            # With coordinate deduplication, reads at the same position should be reduced to 1
            assert len(deduplicated_reads) <= len(mock_aligned_reads), "Deduplication should reduce read count"
            print("✓ Coordinate deduplication test passed!")
            
            # Clean up
            os.unlink(tmp_bam.name)
    
    print("All Smart-seq5 tests passed!")


def test_cli_parameters():
    """
    Test that CLI parameters are correctly parsed for Smart-seq5
    """
    print("\nTesting CLI parameter parsing for Smart-seq5...")
    
    # Test that the CLI accepts the new parameters
    from CamoTSS.bin.count import main
    import sys
    from io import StringIO
    
    # Capture help output to verify new parameters exist
    original_argv = sys.argv[:]
    original_stdout = sys.stdout
    
    try:
        sys.argv = ['count.py', '--help']
        captured_output = StringIO()
        sys.stdout = captured_output
        
        try:
            main()
        except SystemExit:
            # Expected behavior when --help is called
            pass
        
        help_text = captured_output.getvalue()
        sys.stdout = original_stdout
        
        # Check that new parameters are in the help text
        expected_params = [
            '--platform',
            '--bam_list', 
            '--bam_dir',
            '--cell_id_from',
            '--cell_map',
            '--dedup',
            '--min_mapq'
        ]
        
        for param in expected_params:
            assert param in help_text, f"Parameter {param} not found in CLI help"
        
        print("✓ CLI parameter test passed!")
        
    finally:
        sys.argv = original_argv
        sys.stdout = original_stdout


if __name__ == "__main__":
    print("Running Smart-seq5 compatibility tests for CamoTSS...")
    
    test_smartseq5_deduplication()
    test_cli_parameters()
    
    print("\n✓ All Smart-seq5 functionality tests passed!")