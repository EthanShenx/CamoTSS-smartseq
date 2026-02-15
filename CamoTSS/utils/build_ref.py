import pandas as pd
import pyranges as pr
import numpy as np
from functools import reduce 

def _normalize_chrom_series(chrom: pd.Series) -> pd.Series:
    """Normalize contig names so both UCSC-style (chr1) and Ensembl-style (1) GTFs work.

    We store chromosomes without the leading 'chr' when present. Downstream read fetching
    uses check_pysam_chrom() to resolve either naming convention against the BAM header.
    """
    chrom = chrom.copy()
    chrom = chrom.where(~chrom.isna(), other=np.nan)
    chrom = chrom.astype(str)
    chrom = chrom.str.replace(r'^chr', '', regex=True)
    # astype(str) turns NaN into 'nan' string; convert that back to NaN.
    chrom = chrom.replace('nan', np.nan)
    return chrom

def get_TSSref(grdf,out_dir):
    tssdf=grdf[grdf['Feature']=='transcript']
    tssdf=tssdf[['transcript_id','gene_id','gene_name','Chromosome','Start','End','Strand']]
    tssdf['Chromosome']=_normalize_chrom_series(tssdf['Chromosome'])


    #consider positive and negative strand questions
    tssdf['TSS']=np.where(tssdf['Strand']=='+',tssdf['Start'],tssdf['End'])
    tssdf.drop_duplicates(['Chromosome','TSS'],inplace=True)
    tssdf=tssdf[['transcript_id','gene_id','gene_name','Chromosome','Strand','TSS']]
    #tssdf=tssdf[tssdf.groupby('gene_id').gene_id.transform('count')>1]
    tssdf.dropna(subset=['Chromosome'],axis=0,inplace=True)


    TSSoutput_path=out_dir+'ref_TSS.tsv'
    tssdf.to_csv(TSSoutput_path,index=None,sep='\t')
    return TSSoutput_path

def get_generef(grdf,tssdf,out_dir):
    genedf=grdf[grdf['Feature']=='gene']
    genedf=genedf[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name']]

    genedf['Chromosome']=_normalize_chrom_series(genedf['Chromosome'])
    genedf=genedf[genedf['gene_id'].isin(pd.unique(tssdf['gene_id']))]
    genedf.dropna(subset=['Chromosome'],axis=0,inplace=True)

    geneoutput_path=out_dir+'ref_gene.tsv'
    genedf.to_csv(geneoutput_path,index=None,sep='\t')
    return geneoutput_path


def filter_closer_TSS(tssdf,gene_id):
    specificGenedf=tssdf[tssdf['gene_id']==gene_id]
    specificGenedf.sort_values('TSS',inplace=True)
    specificGenedf['diff']=specificGenedf['TSS'].diff()
    specificGenedf=specificGenedf[(specificGenedf['diff']>5)|(specificGenedf['diff'].isna()==True)]
    if len(specificGenedf)>=2:
        return specificGenedf



def get_filter_TSS(tssdf,out_dir):

    keepTSSls=[]
    for geneid in tssdf['gene_id'].unique():
        specificGenedf=filter_closer_TSS(tssdf,geneid)
        keepTSSls.append(specificGenedf)

    filterTSSdf=reduce(lambda x,y:pd.concat([x,y],axis=0) ,keepTSSls)
    filterTSSdf.dropna(subset=['Chromosome'],inplace=True)


    filterTSSoutput_path=out_dir+'filter_TSS_region.tsv'
    filterTSSdf.to_csv(filterTSSoutput_path,index=None,sep='\t')
    
    return filterTSSoutput_path




    

# if __name__ == '__main__':
    # gtf_path="/storage/yhhuang/users/ruiyan/annotation/human_gene_file/Homo_sapiens.GRCh38.105.chr.gtf"
    # gr = pr.read_gtf(gtf_path)
    # grdf = gr.df
    # tssdf= get_TSSref(grdf)
    # genedf=get_generef(grdf,tssdf)
