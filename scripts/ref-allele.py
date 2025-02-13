import pysam
import re
import sys
import os
import urllib.request

def split_repeat_sequence(motifs, sequence):
    """
    Splits repeat sequence based on the defined motifs

    :param motifs: the list of motifs defined for the locus
    :param sequence: sequence of the locus
    :return: the formatted sequence string split based on the motifs
    """
    
    # initialisers
    position = 0    # search checkpoint index within the sequence
    split_seq = ""

    # regular expressions to search for the motifs
    pattern_str = f'({"|".join([x.replace("N", "[ACGT]") for x in motifs])})'
    pattern = re.compile(pattern_str)
    seq_len = len(sequence)

    while position < seq_len:   # search till the end of the sequence

        mpos = pattern.search(sequence[position:])
        if mpos is None: # if no more motifs are found break
            break
        else:
            
            if sequence[position : mpos.start() + position] != "":
                 # sequence before the identified motif is non-empty add it as a split
                split_seq += f"{sequence[position:mpos.start()+position]} "
            
            # adding the identified motif pattern to the formatted string
            split_seq += f"{sequence[position + mpos.start(): position + mpos.start() + (mpos.end() - mpos.start())]} "
            
            # update the search check point
            position += mpos.start() + mpos.end() - mpos.start()

    # add the remaining stretch of the sequence
    if position < seq_len:
        split_seq += f"{sequence[position:]}"
    
    return split_seq

def read_bed(bed, format):
    """
    Reads the bed file and returns the locus information as a dictionary

    :param bed: the bed file
    :param format: the format of the bed file (strchive or pacbio)
    :return: the list of loci information
    """
    
    with open(bed) as fh:
        for line in fh:
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                ref_motif_index = header.index('reference_motif_reference_orientation')
                path_motif_index = header.index('pathogenic_motif_reference_orientation')
                locusid_index = header.index('id')
                continue
            line = line.strip().split('\t')

            if format == 'strchive':
                motifs = list(dict.fromkeys([x.strip().upper() for x in line[ref_motif_index].split(',')] + [x.strip().upper() for x in line[path_motif_index].split(',')]))
                locusid = line[locusid_index]
            elif format == 'pacbio':
                annotations = line[3].split(';')
                motifs = [x.strip().upper() for x in annotations[1].split('=')[1].split(',')]
                locusid = annotations[0].split('=')[1]

            locus = {
                'chrom': line[0],
                'start': int(line[1]),
                'end': int(line[2]),
                'motifs': motifs,
                'id': locusid
            }
            yield locus

def get_ref(fasta, ref_directory='.'):
    """
    Check if the reference genome can be used remotely or if file exists locally, if not download it
    Returns the reference genome object
    """
    # If directory doesn't exist, create it
    if not os.path.isdir(ref_directory):
        os.makedirs(ref_directory)
    # Try to use the ref genome remotely
    try:
        ref = pysam.Fastafile(fasta)
        return ref
    except:
        sys.stderr.write(f"Warning: Couldn't use open the fasta file remotely: {fasta}. Checking for a local copy\n")
        # Download file if it doesn't exist
    
    # Check if the reference genome file exists locally
    ref_file = os.path.basename(fasta)
    ref_path = ref_directory + ref_file
    if os.path.isfile(ref_path) or os.path.isfile(ref_path.replace('.gz', '')):
        try:
            ref = pysam.Fastafile(ref_path)
            return ref
        except:
            try:
                ref = pysam.Fastafile(ref_path.replace('.gz', ''))
                return ref
            except:
                pass
    # Download the reference genome file
    sys.stderr.write(f"Warning: couldn't find a local copy of reference genome file: {ref_path}. Downloading.\n")
    urllib.request.urlretrieve(fasta, ref_path)
    try:
        ref = pysam.Fastafile(ref_path)
        return ref
    except:
        if ref_path.endswith('.gz'):
            # try to unzip the file
            os.system(f"gzip -d {ref_path}")
            ref = pysam.Fastafile(ref_path.replace('.gz', ''))
            return ref
    # If the file is not found, raise an error
    raise FileNotFoundError(f"Reference genome file not found and couldn't be downloaded: {fasta}")

def main(bed1: str, bed2: str, fasta: str, out: str, storage: str = '.', flank: int = 10):
    """
    Extracts the flanking sequences and the repeat sequence for tandem repeat loci from the reference genome

    :param bed1: Input bed file containing the tandem repeat loci from STRchive
    :param bed2: Input bed file containing the tandem repeat loci from PacBio
    :param fasta: Reference genome fasta file
    :param output: Output file name for the extracted sequences
    :param storage: Reference genome storage directory (default: current directory)
    :param flank: Flanking sequence length to extract (default: 10)
    """
    ref = get_ref(fasta, storage)
    sys.stderr.write(f"Successfully accessed reference genome: {os.path.basename(fasta)}\n")

    # Assuming the bed files are sorted and contain the same loci
    with open(out, 'w') as outfh:
        for strchive, pacbio in zip(read_bed(bed1, 'strchive'), read_bed(bed2, 'pacbio')):
            if strchive['id'] != pacbio['id']:
                raise ValueError(f"IDs do not match: {strchive['id']} != {pacbio['id']}")

            # choose the start of the left flank based on the smallest start coordinate
            lflank_start = strchive['start'] - flank if strchive['start'] < pacbio['start'] else pacbio['start'] - flank
            
            # choose the end of the right flank based on the largest end coordinate
            rflank_end = strchive['end'] + flank if strchive['end'] > pacbio['end'] else pacbio['end'] + flank
            
            # fetching the sequences
            # NOTE: the flank sequences are extracted from the start and end of strchive locus
            lflank            = ref.fetch(strchive['chrom'], lflank_start, strchive['start']).upper()
            rflank            = ref.fetch(strchive['chrom'], strchive['end'], rflank_end).upper()
            strchive_seq = ref.fetch(strchive['chrom'], strchive['start'], strchive['end']).upper()
            pacbio_seq   = ""

            # flank sequence initialisers
            strchive_lflank = ""; strchive_rflank = ""
            pacbio_lflank = ""; pacbio_rflank = ""
                
            if strchive['start'] > pacbio['start']:
                # if strchive start is greater than pacbio start
                diff = strchive['start'] - pacbio['start']
                # pull difference bases into the repeat for pacbio and add gaps in the flank
                pacbio_seq    = lflank[-diff:] + strchive_seq
                pacbio_lflank = lflank[:-diff] + " "*diff
                # add difference gaps in the repeat for pacbio and unchanged flank
                strchive_seq  = " "*diff + strchive_seq
                strchive_lflank = lflank
                
            elif strchive['start'] < pacbio['start']:
                # if pacbio start is greater than strchive start
                diff = pacbio['start'] - strchive['start']
                # add difference gaps to pacbio repeat sequence and different bases to flank
                pacbio_seq = " "*diff + strchive_seq[diff:]
                pacbio_lflank = lflank + strchive_seq[:diff]
                # add difference gaps to strchive flank 
                strchive_lflank = lflank[:-diff] + " "*diff
            else:
                # start coordinates are same
                pacbio_lflank   = lflank
                strchive_lflank = lflank
                pacbio_seq = strchive_seq
                
            if strchive['end'] < pacbio['end']:
                # if strchive end is smaller than pacbio end
                diff = pacbio['end'] - strchive['end']
                # add the difference based from flank to pacbio repeat and gaps in flank
                pacbio_seq = pacbio_seq + rflank[:diff]
                pacbio_rflank = " "*diff + rflank[diff:]
                # add the difference gaps strchive repeat and unchanged flank
                strchive_seq += " "*diff
                strchive_rflank = rflank
            
            elif pacbio['end'] < strchive['end']:
                # if pacbio end is smalled than strchive end
                diff = strchive['end'] - pacbio['end']
                # remove difference bases from pacbio repeat and add gaps and add difference bases to the flank
                pacbio_seq = pacbio_seq[:-diff] + " "*diff
                pacbio_rflank = pacbio_seq[-diff:] + rflank
                # unchanged strchive repeat sequence and add difference gaps to strchive flank
                strchive_rflank = " "*diff + rflank
            else:
                # if the end coordinates are the same
                strchive_rflank = rflank
                pacbio_rflank = rflank


            strchive_seq = split_repeat_sequence(strchive['motifs'], strchive_seq)
            pacbio_seq = split_repeat_sequence(pacbio['motifs'], pacbio_seq)
            
            outfh.write(f"{strchive['id']}\n")
            outfh.write(f"{strchive['chrom']}\t{strchive['start']}\t{strchive['end']}\t{','.join(strchive['motifs'])}\tSTRchive\n")
            outfh.write(f"{pacbio['chrom']}\t{pacbio['start']}\t{pacbio['end']}\t{','.join(pacbio['motifs'])}\tTRGT\n")
            outfh.write(f"{strchive_lflank}\t{strchive_seq}\t{strchive_rflank}\n")
            outfh.write(f"{pacbio_lflank}\t{pacbio_seq}\t{pacbio_rflank}\n")
            outfh.write("\n")

    ref.close()
     
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)