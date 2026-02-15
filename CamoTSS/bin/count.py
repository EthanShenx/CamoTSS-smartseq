from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.build_ref import get_TSSref,get_generef,get_filter_TSS
from ..utils.get_counts import get_TSS_count
from ..utils.get_ctss import get_CTSS_count
import pyranges as pr
import os
import pandas as pd
import time 


START_TIME = time.time()


def main():
    parser = OptionParser()
    parser.add_option('--gtf','-g',dest='gtf_file',default=None,help='The annotation gtf file for your analysing species.')
    parser.add_option('--cellbarcodeFile','-c',dest='cdrFile',default=None,help='The file include cell barcode which users want to keep in the downstream analysis.')
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from Cellranger or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') #what should be after $
    parser.add_option('--refFasta','-r',dest='refFasta',default=None,help='The directory for reference genome fasta file') #what should be after $
    parser.add_option('--mode','-m',dest='mode',default=None,help='You can select run by finding novel TSS cluster and CTSS within one cluster [TC+CTSS]. \
                        If you just want to detect TSS cluster, you can use [TC] mode. If you just want to detect CTSS, you can use [CTSS] mode which is based on the output of [TC mode]')
    
    # Smart-seq5 specific options
    parser.add_option('--platform', dest='platform', default='10x', 
                      choices=['10x', 'smartseq5'], 
                      help='Sequencing platform: 10x (default) or smartseq5')
    parser.add_option('--bam_list', dest='bam_list', default=None,
                      help='File containing list of BAM files for smartseq5 mode (one BAM per line)')
    parser.add_option('--bam_dir', dest='bam_dir', default=None,
                      help='Directory containing BAM files for smartseq5 mode (each BAM represents one cell)')
    parser.add_option('--cell_id_from', dest='cell_id_from', default='filename',
                      choices=['filename', 'tsv'], 
                      help='How to determine cell ID for smartseq5: from BAM filename (default) or from TSV mapping')
    parser.add_option('--cell_map', dest='cell_map', default=None,
                      help='TSV file mapping sample names to cell IDs for smartseq5 mode')
    parser.add_option('--dedup', dest='dedup_method', default=None,
                      choices=['umi', 'coord', 'fragment', 'none'],
                      help='Deduplication method: umi (for 10x), coord/fragment (for smartseq5), none. Default depends on platform.')
    parser.add_option('--min_mapq', type="int", dest='min_mapq', default=20,
                      help='Minimum mapping quality for reads [default: 20]')
    parser.add_option('--tss_read', dest='tss_read', default=None,
                      choices=['read1', 'read2'],
                      help="Which mate contains the 5' transcript sequence used for TSS calling. "
                           "Default: read1 for 10x, read2 for smartseq5.")

   
   
    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--minCount",type="int",dest="minCount",default=50,
    help="Minimum UMI counts for TC in all cells [default: 50]")

    group0.add_option('--nproc','-p',type="int",dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')

    group0.add_option('--maxReadCount',type="int",dest='maxReadCount',default=10000,
    help='For each gene, the maxmium read count kept for clustering [default: 10000]')
    
    group0.add_option('--clusterDistance',type="float",dest='clusterDistance',default=300,
    help="The minimum distance between two cluster transcription start site [default: 300]")

    group0.add_option('--InnerDistance',type="float",dest='InnerDistance',default=100,
    help="The resolution of each cluster [default: 100]")

    group0.add_option('--windowSize',type="int",dest='windowSize',default=15,
    help="The width of sliding window [default: 15]")



    group1=OptionGroup(parser,"Optional arguments")

    group1.add_option('--minCTSSCount',type="float",dest='minCTSSCount',default=100,
    help="The minimum UMI counts for each CTSS [default: 100]")

    group1.add_option('--minFC',type="float",dest='minFC',default=6,
    help="The minimum fold change for filtering CTSS [default: 6]")



    parser.add_option_group(group0)
    parser.add_option_group(group1)


    (options, args) = parser.parse_args()

    
    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to CamoTSS v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)





    if options.mode is None:
        print("Error: Need --mode to select the mode what you prefer.")
        sys.exit(1)

    # Determine deduplication method based on platform if not explicitly set
    if options.dedup_method is None:
        if options.platform == '10x':
            options.dedup_method = 'umi'
        else:  # smartseq5
            options.dedup_method = 'coord'  # Default for smartseq5
    
    # Determine which read carries the 5' transcript information
    if options.tss_read is None:
        options.tss_read = 'read1' if options.platform == '10x' else 'read2'

    if (options.mode=='TC') or (options.mode=='TC+CTSS'):
        if options.platform == '10x':
            # Original 10x validation
            if options.cdrFile is None:
                print("Error: Need --cdrFile for cell barcode file.")
                sys.exit(1)

            if options.refFasta is None:
                print("Error: Need --refFasta for reference fasta file.")
                sys.exit(1)

            #bam file
            if options.bam_file is None:
                print("Error: Need --bam for aligned file.")
                sys.exit(1)

            #output file
            if options.out_dir is None:
                print("Warning: no outDir provided, we use $bamfilePath/CamoTSS\n")
                out_dir = os.path.dirname(os.path.abspath(options.bam_file)) + "/CamoTSS"
            else:
                out_dir = options.out_dir
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            bam_file = options.bam_file  # Single BAM file for 10x

        elif options.platform == 'smartseq5':
            # Smart-seq5 validation
            if options.refFasta is None:
                print("Error: Need --refFasta for reference fasta file.")
                sys.exit(1)

            # Either bam_list or bam_dir must be provided for smartseq5
            if options.bam_list is None and options.bam_dir is None:
                print("Error: For smartseq5 platform, need either --bam_list or --bam_dir")
                sys.exit(1)

            # If using cell_map, validate it exists
            if options.cell_id_from == 'tsv' and options.cell_map is None:
                print("Error: For smartseq5 with --cell_id_from tsv, need --cell_map")
                sys.exit(1)

            #output file
            if options.out_dir is None:
                print("Warning: no outDir provided, using ./CamoTSS_smartseq5/\n")
                out_dir = "./CamoTSS_smartseq5"
            else:
                out_dir = options.out_dir
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            # Process BAM files for smartseq5
            bam_files = []
            if options.bam_list:
                with open(options.bam_list, 'r') as f:
                    bam_files = [line.strip() for line in f if line.strip()]
            elif options.bam_dir:
                import glob
                bam_files = glob.glob(os.path.join(options.bam_dir, "*.bam")) + \
                           glob.glob(os.path.join(options.bam_dir, "*.sam"))
            
            print(f"[CamoTSS] Found {len(bam_files)} BAM files for smartseq5 analysis")
            if len(bam_files) == 0:
                print("Error: No BAM files found for smartseq5 analysis")
                sys.exit(1)

            # Create cell barcode list from BAM files or mapping
            if options.cell_id_from == 'filename':
                cell_barcodes = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]
            elif options.cell_id_from == 'tsv':
                cell_map_df = pd.read_csv(options.cell_map, sep='\t')
                # Map sample names (BAM basenames) to cell IDs
                bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]
                cell_barcodes = []
                for basename in bam_basenames:
                    matching_cell = cell_map_df[cell_map_df.iloc[:, 0] == basename]
                    if len(matching_cell) > 0:
                        cell_barcodes.append(matching_cell.iloc[0, 1])
                    else:
                        cell_barcodes.append(basename)  # fallback to filename if not in mapping

            # Create temporary cell barcode file
            temp_cell_file = os.path.join(out_dir, "temp_smartseq5_cells.tsv")
            temp_cell_df = pd.DataFrame({"cell_id": cell_barcodes})
            temp_cell_df.to_csv(temp_cell_file, sep='\t', index=False)
            options.cdrFile = temp_cell_file  # Override for downstream processing

            # For smartseq5, we'll pass the list of BAM files differently
            bam_file = bam_files  # List of BAM files for smartseq5

        #gtf file
        if options.gtf_file is None:
            print("Error: Need --gtf for annotation file.")
            sys.exit(1)
        else:
            gr = pr.read_gtf(options.gtf_file)
            grdf = gr.df
            ref_out_dir=str(out_dir)+'/ref_file/'
            if not os.path.exists(ref_out_dir):
                os.mkdir(ref_out_dir)
            tssrefpath=get_TSSref(grdf,ref_out_dir)
            tssdf=pd.read_csv(tssrefpath,delimiter='\t')
            generefpath=get_generef(grdf,tssdf,ref_out_dir)


    elif options.mode=='CTSS':
        if options.out_dir is None:
            print("Error: Need --outdir which includes subdir '/ref_file' and '/count'")
            sys.exit(1)

        out_dir=options.out_dir


    bam_file = bam_file  # Already set based on platform above
    minCount = options.minCount
    cellBarcodePath = options.cdrFile
    n_proc = options.nproc
    maxReadCount = options.maxReadCount
    clusterDistance = options.clusterDistance
    InnerDistance = options.InnerDistance
    fastqFilePath = options.refFasta
    windowSize = options.windowSize
    minCTSSCount = options.minCTSSCount
    minFC = options.minFC
    platform = options.platform
    dedup_method = options.dedup_method
    min_mapq = options.min_mapq
    




        
    if options.mode == "TC":
        getTSScount=get_TSS_count(generefpath,tssrefpath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance,InnerDistance,windowSize,minCTSSCount,minFC,platform,dedup_method,min_mapq,options.tss_read)
        scadata=getTSScount.produce_sclevel()

    elif options.mode=="TC+CTSS":
        # ctss_out_dir=str(options.out_dir)+'/CTSS/'
        # if not os.path.exists(ctss_out_dir):
        #     os.mkdir(ctss_out_dir)
        getTSScount=get_TSS_count(generefpath,tssrefpath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance,InnerDistance,windowSize,minCTSSCount,minFC,platform,dedup_method,min_mapq,options.tss_read)
        scadata=getTSScount.produce_sclevel()
        twoctssadata=getTSScount.produce_CTSS_adata()


    elif options.mode=='CTSS':
        getctsscount=get_CTSS_count(out_dir,minCTSSCount,minFC,n_proc,windowSize,platform,dedup_method,min_mapq)   # should create CTSS
        ctssadata=getctsscount.produce_CTSS_adata()


    else:
        print('Do not have this mode. Please check your spell!')
        run_time = time.time() - START_TIME
        print("[CamoTSS] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))

