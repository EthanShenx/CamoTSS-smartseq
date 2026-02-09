# Contributing Smart-seq5 Support to CamoTSS

## Overview
This repository contains a complete implementation of Smart-seq5 support for the CamoTSS tool. The implementation has been thoroughly tested and documented.

## Changes Made
The implementation adds Smart-seq5 support while maintaining full backward compatibility with 10x data. Key features include:

- New CLI parameters for platform selection and BAM file handling
- Coordinate/fragment-based deduplication for Smart-seq5 data
- Proper 5' end position calculation for both DNA strands
- Quality control filters for mapping quality and read flags
- Comprehensive documentation and examples

## How to Contribute These Changes

### Method 1: Using GitHub Pull Request (Recommended)
1. **Fork the repository** on GitHub:
   - Go to https://github.com/EthanShenx/CamoTSS-smartseq
   - Click "Fork" button in the top-right corner

2. **Clone your fork**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/CamoTSS-smartseq.git
   cd CamoTSS-smartseq
   ```

3. **Apply the changes using the patch files**:
   ```bash
   git am 0001-Add-Smart-seq5-support-to-CamoTSS.patch
   git am 0002-Add-comprehensive-Smart-seq5-implementation-with-tes.patch
   ```

4. **Push the changes to your fork**:
   ```bash
   git push origin main
   ```

5. **Create a Pull Request**:
   - Go to your fork on GitHub
   - Click "Pull requests" -> "New pull request"
   - Select your branch and create the PR

### Method 2: Manual Application
If patches don't apply cleanly, manually apply the changes:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/EthanShenx/CamoTSS-smartseq.git
   cd CamoTSS-smartseq
   ```

2. **Apply changes manually** to the following files:
   - `CamoTSS/bin/count.py` - Add CLI parameters and logic
   - `CamoTSS/utils/get_counts.py` - Add Smart-seq5 processing
   - `CamoTSS/utils/get_ctss.py` - Update constructor
   - `README.rst` - Add documentation

3. **Add and commit the test files**:
   - `example_smartseq5.sh` - Usage examples
   - `test_smartseq5.py` - Basic tests
   - `comprehensive_smartseq5_test.py` - Comprehensive tests
   - `integration_test.py` - Integration tests

4. **Push and create PR** as described above

## Verification Steps
Before submitting the PR, verify that:

1. **All tests pass**:
   ```bash
   python test_smartseq5.py
   python comprehensive_smartseq5_test.py
   python integration_test.py
   ```

2. **Backward compatibility is maintained**:
   - 10x data processing works as before
   - Default behavior unchanged

3. **Smart-seq5 functionality works**:
   - Multiple BAM files processed correctly
   - Coordinate-based deduplication works
   - Cell IDs extracted properly from filenames or mapping file

## Files Changed
- **Modified**: `CamoTSS/bin/count.py`, `CamoTSS/utils/get_counts.py`, `CamoTSS/utils/get_ctss.py`, `README.rst`
- **Added**: `example_smartseq5.sh`, `test_smartseq5.py`, `comprehensive_smartseq5_test.py`, `integration_test.py`, `SMARTSEQ5_IMPLEMENTATION_SUMMARY.md`

## Usage Examples
Once merged, users can run Smart-seq5 analysis with:
```bash
# With directory of BAM files
CamoTSS --platform smartseq5 --gtf annotation.gtf --refFasta reference.fasta --bam_dir /path/to/bams/ -o output/ --mode TC+CTSS

# With BAM list file
CamoTSS --platform smartseq5 --gtf annotation.gtf --refFasta reference.fasta --bam_list bam_list.txt -o output/ --mode TC+CTSS
```

## Technical Details
- Coordinate deduplication: `(chrom, five_prime_pos, strand)`
- Fragment deduplication: `(chrom, fragment_start, fragment_end, strand)`
- 5' position calculation: `+ strand: reference_start`, `- strand: reference_end - 1`
- Quality filtering: `MAPQ >= 20`, exclude unmapped/secondary/supplementary reads

## Testing Coverage
- Unit tests for deduplication logic
- Integration tests for complete workflow
- Backward compatibility verification
- Edge case handling tests
- 5' position calculation verification

The implementation is production-ready and thoroughly tested.