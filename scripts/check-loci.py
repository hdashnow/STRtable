# Validate the STRchive loci JSON file and overwrite it with changes if needed.
# hg38 is assumed to be the default genome build, with hg19 and T2T coordinates being updated as liftovers.

import sys
import argparse
import doctest
import os
import time
import json
import jsbeautifier
from Bio.Seq import Seq
#import pandas as pd
#import re
from pyliftover import LiftOver
import urllib

# This should be overwritten if schema json file provided
list_fields = [
    "reference_motif_reference_orientation",
    "pathogenic_motif_reference_orientation",
    "pathogenic_motif_gene_orientation",
    "benign_motif_reference_orientation",
    "benign_motif_gene_orientation",
    "unknown_motif_reference_orientation",
    "unknown_motif_gene_orientation",
    "inheritance",
    "source",
    "omim",
    "stripy",
    "gnomad",
    "genereviews",
    "mondo",
    "medgen",
    "orphanet",
    "gard",
    "malacard",
    "webstr_hg38",
    "webstr_hg19",
    "tr_atlas",
    "locus_tags",
    "disease_tags",
    "references",
    "additional_literature",
]

citation_fields = ["references", "additional_literature"]

def parse_args():
    parser = argparse.ArgumentParser(description='Validate the STRchive loci JSON file and overwrite it with changes if needed.')
    parser.add_argument('json', default='data/STRchive-loci.json', help='STRchive loci JSON input file path')
    parser.add_argument('--schema', default=None, help='JSON schema file path (optional)')
    parser.add_argument('--out', default=None, help='Defaults to same as input JSON file path (overwriting)')
    parser.add_argument('--pause', default=5, help='Pause time in seconds before overwriting the file. Default: 5')
    parser.add_argument('--lit', default=None, help='json file of literature with at least the fields "id", "additional_literature", "references". This will overwrite these literature fields in the STRchive loci json if they are present')
    return parser.parse_args()

def unique_list(mylist):
    """
    Args:
        lst (list)
    Returns:
        list: unique values in lst
    >>> unique_list([1, 2, 2, 3, 3, 3])
    [1, 2, 3]
    >>> unique_list([])
    []
    >>> unique_list([1, 2, 3])
    [1, 2, 3]
    >>> unique_list([3, 3, 1, 2, 3])
    [3, 1, 2]
    >>> unique_list(['a', 'lazy', 'quick', 'brown', 'fox', 'jumps', 'over', 'the', 'lazy', 'dog'])
    ['a', 'lazy', 'quick', 'brown', 'fox', 'jumps', 'over', 'the', 'dog']
    """
    unique_list = []
    seen = set()

    for item in mylist:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)

    return unique_list

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All circular permutations of x
    >>> circular_permuted('GAG')
    ['GAG', 'AGG', 'GGA']
    >>> circular_permuted('AG')
    ['AG', 'GA']
    >>> circular_permuted('TAGAA')
    ['TAGAA', 'AGAAT', 'GAATA', 'AATAG', 'ATAGA']
    >>> circular_permuted('')
    []
    """
    n = len(x)
    modified_sequences = []
    modified_sequences.extend([x[i:] + x[:i] for i in range(n)])
    return modified_sequences

def normalise_str(in_dna):
    """
    Args:
        in_dna (sequence)
    Returns:
        the normalized output of the string
    Find all possible equivalent STR sequences and return the first alphabetically for each
    >>> normalise_str('ATAG')
    'AGAT'
    >>> normalise_str('NGC')
    'CNG'
    >>> normalise_str('AAG')
    'AAG'
    """
    if in_dna is None or len(in_dna) == 0:
        return ''

    all_possible = []
    for permuted_seq in circular_permuted(in_dna):
        all_possible.append(permuted_seq)

    return min(all_possible)

def get_new_motif(motif, gene_strand):
    """
    Args:
        motif (string)
        gene_strand: either + or -
    Returns:
        the normalized output of the string from ref to gene orientation
    Get the new normalized motif for each row.
    If gene_strand is +, reference orientation = gene orientation
    If gene_strand is -, reverse_complement ref_ori for gene_ori
    >>> get_new_motif('GAG', '+')
    'AGG'
    >>> get_new_motif('GAG', '-')
    'CCT'
    >>> get_new_motif('TCATC', '-')
    'AGATG'
    >>> get_new_motif('TAG', 'plus')
    Traceback (most recent call last):
    ...
    AssertionError: Gene strand plus is not +/-
    """
    if gene_strand == "+":
        normalized_motif = normalise_str(motif)
    elif gene_strand == "-":
        seq = Seq(motif)
        reverse_comp = str(seq.reverse_complement())
        normalized_motif = normalise_str(reverse_comp)
    else:
        raise AssertionError(f'Gene strand {gene_strand} is not +/-')
    return normalized_motif
    
def check_motif_orientation(record):
    """
    Args:
        record (dict): a dictionary containing a single locus from the STRchive json
    Returns:
        record (dict): the record with any motif fields with incorrect orientation updated
    """
    field_pairs = [
        ('pathogenic_motif_reference_orientation', 'pathogenic_motif_gene_orientation'),
        ('benign_motif_reference_orientation', 'benign_motif_gene_orientation'),
        ('unknown_motif_reference_orientation', 'unknown_motif_gene_orientation')
    ]
    for ref_field, gene_field in field_pairs:
        if record[ref_field] is None:
            continue
        old = record[gene_field]
        new = [get_new_motif(x, record['gene_strand']) for x in record[ref_field]]
        if old != new:
            for old_motif, new_motif in zip(old, new):
                if old_motif != new_motif:
                    sys.stderr.write(f'Updating {record['id']} {gene_field} from {old_motif} to {new_motif}\n')
        record[gene_field] = new
    return record

def check_list_fields(record, list_fields = list_fields):
    """
    Args:
        record (dict): a dictionary containing a single locus from the STRchive json
        list_fields (list): a list of names of fields that should be lists
    Returns:
        record (dict): the record with any fields that should be lists converted to lists
    """
    for field in list_fields:
        try:
            old = record[field]
        except KeyError:
            record[field] = []
            continue
        if record[field] == "" or record[field] is None or record[field] == []:
            record[field] = [] #XXX need to decide if Null or empty list is preferred
        elif isinstance(record[field], str):
            if ";" in record[field]:
                record[field] = [x.strip() for x in record[field].split(';')]
            elif "," in record[field]:
                record[field] = [x.strip() for x in record[field].split(',')]
            else:
                record[field] = [record[field].strip()]
        
        if field in citation_fields:
            # Remove leading @
            record[field] = [x.strip('@') for x in record[field]]
            # Ensure citation lists are unique
            record[field] = unique_list(record[field])
            # Remove empty strings
            record[field] = [x for x in record[field] if x.split(':')[-1] != '']

        if old != record[field]:
            sys.stderr.write(f'Updating {record['id']} {field} from {old} to {record[field]}\n')
    return record

def lift_over(loci_data, ref_dir):

    chainfile_t2t_url = "https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz"
    chainfile_t2t = f"{ref_dir}/{os.path.basename(chainfile_t2t_url)}"
    # If directory doesn't exist, create it
    if not os.path.isdir(ref_dir):
        os.makedirs(ref_dir)
    # Check if chain file already exists and download it if not
    if not os.path.exists(chainfile_t2t):
        sys.stderr.write(f"Downloading chain file for T2T from {chainfile_t2t_url}\n")
        urllib.request.urlretrieve(chainfile_t2t_url, chainfile_t2t)

    lo_t2t = LiftOver(chainfile_t2t)

    # Load chain file for hg38 to hg19 conversion
    lo_hg38_to_hg19 = LiftOver('hg38', 'hg19')

    # Process each entry in the JSON
    for entry in loci_data:
        chrom = entry["chrom"]

        # Convert start and stop positions from hg38 to T2T
        converted_start_t2t = lo_t2t.convert_coordinate(chrom, entry["start_hg38"])
        converted_stop_t2t = lo_t2t.convert_coordinate(chrom, entry["stop_hg38"])

        # Convert start and stop positions from hg38 to hg19
        converted_start_hg19 = lo_hg38_to_hg19.convert_coordinate(chrom, entry["start_hg38"])
        converted_stop_hg19 = lo_hg38_to_hg19.convert_coordinate(chrom, entry["stop_hg38"])

        # Store results in the JSON
        if converted_start_t2t and converted_stop_t2t:
            entry["start_t2t"] = converted_start_t2t[0][1]
            entry["stop_t2t"] = converted_stop_t2t[0][1]
        else:
            sys.stderr.write(f"WARNING: No T2T coordinates found for {entry['id']} {chrom}:{entry['start_hg38']}-{entry['stop_hg38']}\n")

        if converted_start_hg19 and converted_stop_hg19:
            entry["start_hg19"] = converted_start_hg19[0][1]
            entry["stop_hg19"] = converted_stop_hg19[0][1]
        else:
            sys.stderr.write(f"WARNING: No hg19 coordinates found for {entry['id']} {chrom}:{entry['start_hg38']}-{entry['stop_hg38']}\n")

    return loci_data

def main(json_fname, json_schema = None, out_json = None, pause = 5, lit = None, ref_dir = "ref-data/"):
    if out_json == json_fname:
        sys.stderr.write(f'WARNING: overwriting {json_fname} in {pause} seconds\n')
        sys.stderr.write('Press Ctrl+C to cancel\n')
        time.sleep(pause)
    else:
        sys.stderr.write(f'Writing {out_json}\n')

    # Check if file exists
    if not os.path.exists(json_fname):
        sys.stderr.write(f'Error: {json_fname} does not exist\n')
        sys.exit(1)
    # Read JSON file
    with open(json_fname, 'r') as json_file:
        data = json.load(json_file)

        if lit:
            with open(lit, 'r') as lit_file:
                lit_data = json.load(lit_file)
                lit_dict = {x['id']: x for x in lit_data}
                for record in data:
                    if record['id'] in lit_dict:
                        record['references'] = lit_dict[record['id']]['references']
                        record['additional_literature'] = lit_dict[record['id']]['additional_literature']

        # Check if the field contains a string that should be a list
        for record in data:
            record = check_list_fields(record)
            record = check_motif_orientation(record)

        # Lift over coordinates
        data = lift_over(data, ref_dir)

        # Sort records by gene name then id
        data = sorted(data, key = lambda x: (x['gene'], x['id']))
        
        # Write JSON file
        with open(out_json, 'w') as out_json_file:
            options = jsbeautifier.default_options()
            options.indent_size = 2
            options.brace_style="expand"
            out_json_file.write(jsbeautifier.beautify(json.dumps(data, ensure_ascii=False), options))

if __name__ == '__main__':
    doctest.testmod()
    args = parse_args()
    if args.out is None:
        args.out = args.json
    main(args.json, args.schema, args.out, args.pause, args.lit)
