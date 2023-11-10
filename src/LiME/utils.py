import argparse


def _parse_args(args, **kwargs):

    parser = argparse.ArgumentParser(
                        prog='Linked_Mutation_Extractor',
                        description=('Classify aligned reads basing on different'
                         			 'nucleotide combinations on a series of' 
                         			 'predefined positions'),
                        epilog='')
    
    parser.add_argument('-i',
                        type=str,
                        help='Import a BAM file.')  
    parser.add_argument('--vcf',
                        type=str,
                        help='Provide positions from a VCF file.')  
    parser.add_argument('--txt',
                        type=str,
                        help='Provide positions from a txt file.')     
    parser.add_argument('-p',
    					default = None,
                        type=str,
                        help='Define the output path for the waffle plot')
    parser.add_argument('-f',
    					default = None,
                        type=str,
                        help='Define the output path for the CSV file')

    return parser.parse_args(args)