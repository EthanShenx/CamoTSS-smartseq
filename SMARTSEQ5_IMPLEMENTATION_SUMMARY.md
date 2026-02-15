# Smart-seq5 Support Implementation Summary

## Overview
This implementation adds Smart-seq5 support to the CamoTSS tool while maintaining full backward compatibility with 10x data.

## Key Changes Made

### 1. CLI Parameter Updates (`CamoTSS/bin/count.py`)
- Added `--platform {10x,smartseq5}` parameter with default '10x'
- Added `--bam_list` and `--bam_dir` for multiple BAM file input
- Added `--cell_id_from {filename,tsv}` for cell ID determination
- Added `--cell_map` for custom cell ID mapping
- Added `--dedup {umi,coord,fragment,none}` with platform-dependent defaults
- Added `--min_mapq` for mapping quality filtering
- Added `--tss_read {read1,read2}` to select which mate contains 5' transcript information for TSS calling

### 2. Core Logic Updates (`CamoTSS/utils/get_counts.py`)
- Updated `get_TSS_count` class to accept new parameters
- Rewrote `_get_gene_reads()` method to handle both platforms
- Added `_deduplicate_by_umi()` for 10x-style deduplication
- Added `_deduplicate_by_coordinates()` for Smart-seq5 deduplication
- Implemented proper 5' end position calculation for both strands
- Added quality filtering for reads

### 3. CTSS Module Updates (`CamoTSS/utils/get_ctss.py`)
- Updated `get_CTSS_count` class to accept new parameters

### 4. Documentation Updates (`README.rst`)
- Added comprehensive Smart-seq5 usage documentation
- Included example commands for different use cases
- Documented all new parameters

## Smart-seq5 Specific Features

### A. Coordinate-based Deduplication
- `coord` method: Uses (chrom, five_prime_pos, strand) as dedup key
- `fragment` method: Uses (chrom, fragment_start, fragment_end, strand) for paired-end data
- Proper 5' end position calculation for both strands (+: reference_start, -: reference_end-1)

### B. Cell Processing
- Each BAM file treated as a single cell
- Cell ID derived from BAM filename by default
- Optional custom mapping via TSV file

### C. Quality Control
- Filtering for unmapped, secondary, duplicate, and supplementary reads
- Mapping quality threshold filtering
- Strand-specific processing maintained

## Backward Compatibility
- All existing 10x functionality preserved
- Default behavior unchanged for 10x data
- Same output format maintained (cell × TSS/peak matrices)

## Testing
- Comprehensive tests for deduplication logic
- Integration tests for complete workflow
- Backward compatibility verification
- Edge case handling tests

## Files Modified
- `CamoTSS/bin/count.py` - CLI parameters and main logic
- `CamoTSS/utils/get_counts.py` - Core counting logic
- `CamoTSS/utils/get_ctss.py` - CTSS processing
- `README.rst` - Documentation updates
- `example_smartseq5.sh` - Usage examples
- `test_smartseq5.py` - Basic functionality tests
- `comprehensive_smartseq5_test.py` - Comprehensive tests
- `integration_test.py` - Integration tests

## Usage Examples

### Smart-seq5 with directory of BAM files:
```bash
CamoTSS --platform smartseq5 --gtf annotation.gtf \
        --refFasta reference.fasta \
        --bam_dir /path/to/bam/files/ \
        -o ./output_smartseq5 \
        --mode TC+CTSS
```

### Smart-seq5 with BAM list file:
```bash
CamoTSS --platform smartseq5 --gtf annotation.gtf \
        --refFasta reference.fasta \
        --bam_list bam_list.txt \
        -o ./output_smartseq5 \
        --mode TC+CTSS
```

## Validation
- All tests pass successfully
- Smart-seq5 deduplication reduces duplicate reads appropriately
- 5' position calculation works for both DNA strands
- Backward compatibility with 10x data maintained
- Complete workflow simulation successful
