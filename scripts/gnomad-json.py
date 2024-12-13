import sys
import pandas as pd
import json
from statsmodels.stats import proportion
import argparse

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate plotly data (JSON) for pathogenic size ranges")
    parser.add_argument("database_json", help="JSON file of STRchive loci database")
    parser.add_argument("gnomad_tsv", help="TSV file with gnomAD genotype data")
    parser.add_argument("output", help="Output JSON file for plotly library")
    return parser.parse_args()

def cyclical_variations(motifs):
    # builds the set of all possible motifs with a given motif
    # Motifs with N in them are all permuted with combination of all nucleotides
    nucs = ['A', 'C', 'G', 'T']
    variations = []
    for motif in motifs:
        if 'N' in motif:
            for n in nucs:
                variations += cyclical_variations([motif.replace('N', n)])
        variations += [motif[-i:]+motif[:-i] for i in range(len(motif))]
    return variations


def binomial_ci(x, n, confidence=0.95):
    """
    Calculates the confidence interval for a binomial proportion.

    Args:
        x: Number of successes
        n: Number of trials
        confidence: Confidence level (e.g., 0.95 for a 95% CI)

    Returns:
        tuple: Lower and upper bounds of the confidence interval
    """
    
    lower, upper = proportion.proportion_confint(x, n, alpha=1-confidence, method='binom_test')
    return lower, upper

def calc_plotdata(sub_df):
    sorted_poplabels = ["amr", "afr", "ami", "asj", "nfe", "fin", "mid", "sas", "eas", "oth"]
    population_labels = {
        "amr": "Admixed American",
        "afr": "African/African American",
        "asj": "Ashkenazi Jewish",
        "ami": "Amish",
        "eas": "East Asian",
        "fin": "Finnish",
        "nfe": "European (non Finnish)",
        "mid": "Middle Eastern",
        "sas": "South Asian",
        "oth": "Others",
    }
    labels = []
    values = []
    counts = []
    xconf_lowerbound = []
    xconf_upperbound = []
    for population in sorted_poplabels:
        labels.append(f'{population_labels[population]}')
        count = len(sub_df[sub_df['Population']==population])
        counts.append(count)
        if True in sub_df[sub_df['Population']==population]['Pathogenic'].unique():
            pcount = sub_df[sub_df['Population']==population]['Pathogenic'].value_counts()[True]
            values.append((pcount/count)*100)
            conf_lower, conf_upper = binomial_ci(pcount, count)
        else:
            values.append(0)
            conf_lower, conf_upper = binomial_ci(0, count, 0.95)
        xconf_lowerbound.append(conf_lower*100)
        xconf_upperbound.append(conf_upper*100)
    
    return [values, labels, counts, xconf_lowerbound, xconf_upperbound]

def build_JSON(args):
    json_txt = ""
    with open(args.database_json) as fh:
        for line in fh: json_txt += line.strip()
    strchive_info = json.loads(json_txt)
    
    df = pd.read_csv(args.gnomad_tsv, sep='\t')

    json_out = open(args.output, 'w')

    plot_data = {}
    for info in strchive_info:

        # looks for gnomad specific gene name otherwise takes the default gene name
        try: gene = info['gnomad'][0]
        except IndexError: gene = info['gene']
        
        sub_df = df[df['Id'] == gene]
        
        # gene is not in gnomAD dataset
        if len(sub_df) == 0: continue

        chrom = info['chrom']
        pathogenic_motif = info['pathogenic_motif_reference_orientation']
        try: inheritance = info['inheritance'][0]
        except IndexError: inheritance = ''
        
        pathogenic_min = info['pathogenic_min']
        pathogenic_max = info['pathogenic_max']
        
        prevalence = None
        if info['prevalence'] is not None: prevalence = (float(info['prevalence'].split('/')[0].strip())/float(info['prevalence'].split('/')[1].strip()))*100

        if chrom == 'chrX':
            sub_df_XX = sub_df[sub_df['Sex'] == 'XX']
            sub_df_XX = sub_df_XX.assign(MinAllele1  = sub_df_XX["GenotypeConfidenceInterval"].apply(lambda x: min([int(a) for a in x.split('/')[0].split('-')])))
            sub_df_XX = sub_df_XX.assign(MinAllele2  = sub_df_XX["GenotypeConfidenceInterval"].apply(lambda x: min([int(a) for a in x.split('/')[1].split('-')])))
            sub_df_XX = sub_df_XX.assign(Pathogenic1 = [pathogenic_min < x  for i, x in enumerate(sub_df_XX['MinAllele1'])])
            sub_df_XX = sub_df_XX.assign(Pathogenic2 = [pathogenic_min < x  for i, x in enumerate(sub_df_XX['MinAllele2'])])
            sub_df_XX = sub_df_XX.assign(PathogenicMotif = [ x in cyclical_variations(pathogenic_motif) for x in sub_df_XX['Motif']])

            sub_df_XY = sub_df[sub_df['Sex'] == 'XY']
            sub_df_XY = sub_df_XY.assign(MinAllele1  = sub_df_XY["GenotypeConfidenceInterval"].apply(lambda x: min([int(a) for a in x.split('/')[0].split('-')])))
            sub_df_XY = sub_df_XY.assign(MinAllele2  = 'NA')
            sub_df_XY = sub_df_XY.assign(Pathogenic1 = [pathogenic_min < x  for i, x in enumerate(sub_df_XY['MinAllele1'])])
            sub_df_XY = sub_df_XY.assign(Pathogenic2  = 'NA')
            sub_df_XY = sub_df_XY.assign(PathogenicMotif = [x in cyclical_variations(pathogenic_motif) for x in sub_df_XY['Motif']])
            
            if inheritance == 'XD':
                sub_df_XX['Pathogenic'] = sub_df_XX['Pathogenic1'] | sub_df_XX['Pathogenic2']
                sub_df_XY['Pathogenic'] = sub_df_XY['Pathogenic1']
                sub_df = pd.concat([sub_df_XX, sub_df_XY])

            elif inheritance == 'XR':
                sub_df_XX['Pathogenic'] = sub_df_XX['Pathogenic1'] & sub_df_XX['Pathogenic2']
                sub_df_XY['Pathogenic'] = sub_df_XY['Pathogenic1']
                sub_df = pd.concat([sub_df_XX, sub_df_XY])
            
            plot_data[info["id"]] = {'XX': "", 'XY': ""}

            for sex, plot_df in zip(['XX', 'XY'], [sub_df_XX, sub_df_XY]):
                values, labels, counts, xconf_lowerbound, xconf_upperbound = calc_plotdata(plot_df)
                plot_data[info["id"]][sex] = {
                    "id": info["id"],
                    "labels": labels,
                    "values": values,
                    "counts": counts,
                    "confidence_lowerbounds": xconf_lowerbound,
                    "confidence_upperbounds": xconf_upperbound,
                    "title": gene+"_"+sex
                }

            sub_df = pd.concat([sub_df_XX, sub_df_XY])
            values, labels, counts, xconf_lowerbound, xconf_upperbound = calc_plotdata(sub_df)

            plot_data[info["id"]]["id"] = info["id"]
            plot_data[info["id"]]["labels"] = labels
            plot_data[info["id"]]["values"] = values
            plot_data[info["id"]]["counts"] = counts
            plot_data[info["id"]]["confidence_lowerbounds"] = xconf_lowerbound
            plot_data[info["id"]]["confidence_upperbounds"] = xconf_upperbound
            plot_data[info["id"]]["title"] = gene
        
        else:
            sub_df = sub_df.assign(MinAllele1  = sub_df["GenotypeConfidenceInterval"].apply(lambda x: min([int(a) for a in x.split('/')[0].split('-')])))
            sub_df = sub_df.assign(MinAllele2  = sub_df["GenotypeConfidenceInterval"].apply(lambda x: min([int(a) for a in x.split('/')[1].split('-')])))

            if gene == 'VWA1':
                sub_df = sub_df.assign(Pathogenic1 = [info['benign_min'] != x  for i, x in enumerate(sub_df['MinAllele1'])])
                sub_df = sub_df.assign(Pathogenic2 = [info['benign_min'] != x  for i, x in enumerate(sub_df['MinAllele2'])])
            else:
                sub_df = sub_df.assign(Pathogenic1 = [pathogenic_min < x  for i, x in enumerate(sub_df['MinAllele1'])])
                sub_df = sub_df.assign(Pathogenic2 = [pathogenic_min < x  for i, x in enumerate(sub_df['MinAllele2'])])
            sub_df = sub_df.assign(PathogenicMotif = [x in cyclical_variations(pathogenic_motif) for x in sub_df['Motif']])
            
            if   inheritance == 'AD': sub_df = sub_df.assign(Pathogenic = (sub_df['Pathogenic1'] | sub_df['Pathogenic2']) & sub_df['PathogenicMotif'])
            elif inheritance == 'AR': sub_df = sub_df.assign(Pathogenic = (sub_df['Pathogenic1'] & sub_df['Pathogenic2']) & sub_df['PathogenicMotif'])
        
            values, labels, counts, xconf_lowerbound, xconf_upperbound = calc_plotdata(sub_df)
            plot_data[info["id"]] = {
                "id": info["id"],
                "labels": labels,
                "values": values,
                "counts": counts,
                "confidence_lowerbounds": xconf_lowerbound,
                "confidence_upperbounds": xconf_upperbound,
                "title": gene
            }
            
    print(json.dumps(plot_data, indent=4), file=json_out)
    json_out.close()

if __name__ == "__main__":
    args = parse_args()
    build_JSON(args)