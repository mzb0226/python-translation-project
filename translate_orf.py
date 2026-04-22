#! /usr/bin/env python3

import sys
import re
from translate import translate_sequence

def vet_nucleotide_sequence(sequence):
    """
    Return None if `sequence` is a valid RNA or DNA sequence, else raise exception.
    """
    rna_pattern_str = r'^[AUCGaucg]*$'
    dna_pattern_str = r'^[ATCGatcg]*$'

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))


def vet_codon(codon):
    """
    Return None if `codon` is a valid RNA codon, else raise an exception.
    """
    codon_pattern_str = r'^[AUGCaugc]{3}$'
    codon_pattern = re.compile(codon_pattern_str)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))


def find_first_orf(sequence,
        start_codons=['AUG'],
        stop_codons=['UAA', 'UAG', 'UGA']):
    """
    Return the first open-reading frame in the DNA or RNA `sequence`.
    """
    vet_nucleotide_sequence(sequence)

    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    seq = seq.replace('T', 'U')

    start_pattern = '|'.join(starts)
    stop_pattern = '|'.join(stops)
    orf_pattern_str = r'(?:' + start_pattern + r')(?:[AUCG]{3})*?(?:' + stop_pattern + r')'

    orf_pattern = re.compile(orf_pattern_str)
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''


def parse_sequence_from_path(path):
    try:
        file_stream = open(path, 'r')
    except FileNotFoundError as e:
        sys.stderr.write("Sorry, couldn't find path {}".format(path))
        raise e
    except IsADirectoryError as e:
        sys.stderr.write("Sorry, path {} appears to be a directory".format(path))
        raise e
    except:
        sys.stderr.write("Sorry, something went wrong when trying to open {}".format(path))
        raise

    sequence = ''
    for line in file_stream:
        sequence += line.strip()
    return sequence


def main():
    import argparse

    parser = argparse.ArgumentParser()

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    parser.add_argument('sequence',
            metavar='SEQUENCE',
            type=str,
            help=('The sequence to search for an open-reading frame. '
                  'If the path flag (\'-p\'/\'--path\') is specified, '
                  'then this should be a path to a file containing the '
                  'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action='store_true',
            help=('The sequence argument should be treated as a path to a '
                  'file containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type=str,
            action='append',
            default=None,
            help=('A start codon. This option can be used multiple times '
                  'if there are multiple start codons. '
                  'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type=str,
            action='append',
            default=None,
            help=('A stop codon. This option can be used multiple times '
                  'if there are multiple stop codons. '
                  'Default: {0}.'.format(" ".join(default_stop_codons))))

    args = parser.parse_args()

    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    orf = find_first_orf(sequence=sequence,
            start_codons=args.start_codon,
            stop_codons=args.stop_codon)

    if not orf:
        sys.exit("No ORF found in sequence")

    amino_acid_sequence = translate_sequence(orf, genetic_code)
    sys.stdout.write('{}\n'.format(amino_acid_sequence))


if __name__ == '__main__':
    main()
