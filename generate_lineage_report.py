import bte
import pandas as pd
import numpy as np
import datetime as dt
import argparse

def argparser():
    parser = argparse.ArgumentParser(description="Compute detailed lineage reports for all existing lineages in the tree.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to compute reports from.')
    parser.add_argument("-p", "--proposed", required=True, help='Path to the file containing dumped sublineage proposals.')
    parser.add_argument("-o", "--output", help='Name of the output table.',default=None,required=True)
    parser.add_argument("-m", "--metadata", help="Path to a metadata file matching the protobuf.",required=True)
    parser.add_argument("-f", "--reference", default=None, help="Path to a reference fasta file. Use with -g to annotate amino acid changes in the expanded output.")
    parser.add_argument("-g", "--gtf", default=None, help="Path to a reference gtf file. Use with -f to annotate amino acid changes in the expanded output.")
    args = parser.parse_args()
    return args

def get_date(d):
    try:
        return dt.datetime.strptime(d,"%Y-%m-%d")
    except:
        return np.nan

def is_successive(row):
    if row.EarliestChild > row.EarliestParent and row.LatestChild > row.LatestParent:
        return True
    else:
        return False

def fill_output_table(t,pdf,mdf,fa_file=None,gtf_file=None):
    mdf.set_index('strain',inplace=True)
    def parent_lineage_size(row):
        samples = t.get_leaves_ids(row.parent_nid)
        return len(samples)
    def sublineage_size(row):
        subsamples = t.get_leaves_ids(row.proposed_sublineage_nid)
        return len(subsamples)
    print("Computing sublineage percent")
    pdf['ParentLineageSize'] = pdf.apply(parent_lineage_size,axis=1)
    pdf['SublineageSize'] = pdf.apply(sublineage_size,axis=1)
    pdf['SublineagePercent'] = pdf.SublineageSize/pdf.ParentLineageSize
    def parsimony_parent(row):
        parent_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.parent_nid)])
        return parent_parsimony
    def parsimony_child(row):
        child_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.proposed_sublineage_nid)])
        return child_parsimony
    print("Computing parsimony percent")
    pdf['ParentParsimony'] = pdf.apply(parsimony_parent,axis=1)
    pdf['SublineageParsimony'] = pdf.apply(parsimony_child,axis=1)
    pdf['ParsimonyPercent'] = pdf.SublineageParsimony/pdf.ParentParsimony
    mdf['date'] = mdf.date.apply(get_date)
    def get_start_ends(row):
        parent_samples = set(t.get_leaves_ids(row.parent_nid))
        child_samples = set(t.get_leaves_ids(row.proposed_sublineage_nid))
        parent_only = parent_samples - child_samples
        try:
            parent_dates = mdf.loc[list(parent_only)].date
            child_dates = mdf.loc[list(child_samples)].date
            return min(parent_dates),max(parent_dates),min(child_dates),max(child_dates)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except:
            return np.nan,np.nan,np.nan,np.nan
    print("Computing start and end dates")
    applied_pdf = pdf.apply(lambda row: get_start_ends(row), axis='columns', result_type='expand')
    pdf = pd.concat([pdf, applied_pdf], axis='columns')
    pdf = pdf.rename({0:'EarliestParent',1:'LatestParent',2:'EarliestChild',3:'LatestChild'},axis=1)
    pdf['LogScore'] = np.log10(pdf.proposed_sublineage_score)
    pdf["Successive"] = pdf.apply(is_successive,axis=1)
    print("Doing international")
    def get_regions(nid):
        try:
            return ",".join(list(mdf.loc[t.get_leaves_ids(nid)].country.value_counts().index))
        except:
            return np.nan
    def get_regions_percents(nid):        
        try:
            return ",".join([str(p) for p in mdf.loc[t.get_leaves_ids(nid)].country.value_counts(normalize=True)])
        except:
            return np.nan
    pdf['ChildRegions'] = pdf.proposed_sublineage_nid.apply(get_regions)
    # pdf['ParentRegions'] = pdf.parent_nid.apply(get_regions)
    pdf['ChildRegionPercents'] = pdf.proposed_sublineage_nid.apply(get_regions_percents)
    # pdf['ParentRegionPercents'] = pdf.parent_nid.apply(get_regions_percents)
    def host_jump(row):
        try:
            return mdf.loc[t.get_leaves_ids(row.proposed_sublineage_nid)].host.nunique() > 1
        except:
            return False
    pdf['HostJump'] = pdf.apply(host_jump,axis=1)
    def generate_url(nid):
        mset = t.get_haplotype(nid)
        url = "https://cov-spectrum.org/explore/World/AllSamples/AllTimes/variants?variantQuery=[" + str(len(mset)) + "-of:"
        start = True
        for m in mset:
            if start:
                start = False
            else:
                url += ', '
            url += m
        url += ']'
        return url
    pdf['Links'] = pdf.proposed_sublineage_nid.apply(generate_url)
    def get_separating_mutations(row):
        hapstring = []
        for n in t.rsearch(row.proposed_sublineage_nid,True):
            if n.id == row.parent_nid:
                break
            hapstring.append(",".join(n.mutations))
        return ">".join(hapstring)
    pdf['Mutations'] = pdf.apply(get_separating_mutations,axis=1)
    if gtf_file != None and fa_file != None:
        def get_separating_translation(row):
            translation = t.translate(fasta_file = fa_file, gtf_file = gtf_file)
            hapstring = []
            for n in t.rsearch(row.proposed_sublineage_nid,True):
                if n.id == row.parent_nid:
                    break
                aas = translation.get(n.id,[])
                hapstring.append(",".join([aav.gene+":"+aav.aa for aav in aas]))
            return ">".join(hapstring)
        pdf['AAChanges'] = pdf.apply(get_separating_translation,axis=1)
    return pdf

def main():
    args = argparser()
    mdf = pd.read_csv(args.metadata,sep='\t')
    t = bte.MATree(args.input)
    pdf = pd.read_csv(args.proposed,sep='\t')
    odf = fill_output_table(t,pdf,mdf,args.reference,args.gtf)
    odf.to_csv(args.output,sep='\t',index=False)

if __name__ == "__main__":
    main()