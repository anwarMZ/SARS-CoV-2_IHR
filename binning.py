#!/usr/bin/env python

import argparse as ap
import pandas as pd
import subprocess
import os
import tempfile
import statistics
import logging


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Runs MASH and computes pairwise distance and distance against the reference')

    parser.add_argument('-ll', '--loglevel', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Set the logging level')
    parser.add_argument('-t', '--time', type=int, default=14, help='Number of days within a bin (default: 14)')
    parser.add_argument('-w', '--window', type=int, default=1, help='Number of days between bins (default: 1)')
    parser.add_argument('-c', '--clean', action='store_true', help='Remove all temporary MASH files after execution')
    parser.add_argument('-k', '--kmer', type=int, default=17, help='kmer size for MASH (default: 17)')
    parser.add_argument('-p', '--threads', type=int, default=1, help='Number of threads to use for MASH (default: 1)')
    parser.add_argument('--prefix', type=str, default='out', help='Output file prefix (default: out)')
    parser.add_argument('meta', metavar='METADATA', help='GISAID metadata')
    parser.add_argument('fasta', metavar='FASTA', help='FASTA of sequences')
    parser.add_argument('ref', metavar='REF', help='Reference genome in FASTA format')  
    parser.add_argument('outdir', metavar='DIR', help='Output directory')

    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) -> %(levelname)s: %(message)s',
                        datefmt='%I:%M:%S %p')

    logger = logging.getLogger()

    # Check if MASH is installed
    mash_rc = subprocess.call(['mash'], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if mash_rc == 127:
        raise FileNotFoundError("Unable to locate \'MASH\'")

    # Check if output directory exists
    tmpdir = os.path.join(args.outdir, 'tmp_mash')
    tmpdir_obj = None
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    if args.clean:
        tmpdir_obj = tempfile.TemporaryDirectory()
        tmpdir = tmpdir_obj.name
    elif not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    # Read and preprocess metadata file
    logger.info('Reading and preprocessing metadata file')
    metadata = pd.read_csv(args.meta, sep='\t')
    metadata['date'] = pd.to_datetime(metadata['date'])
    metadata = metadata.sort_values('date')
    metadata = metadata.reset_index()

    start_date = metadata['date'][0]
    end_date = start_date + pd.DateOffset(days=args.time)

    out_dates = []
    out_ref_dist = []
    out_ref_min = []
    out_ref_max = []
    out_pair_dist = []
    out_pair_min = []
    out_pair_max = []
    bin_size = []

    prefix = os.path.join(tmpdir, args.prefix)

    # Create sketch 
    logger.info('Run MASH to obtain distance against reference')
    subprocess.run(['mash sketch -k {0} -p {1} -o {2}.msh -i {3}'.format(args.kmer, args.threads, prefix, args.fasta)], shell=True, stdout=subprocess.DEVNULL)

    # Run MASH distance calculation
    subprocess.run(['mash dist -t {0}.msh {1} >> {2}.dist'.format(prefix, args.ref, prefix)], shell=True, stdout=subprocess.DEVNULL)

    ref_dist_mat = pd.read_csv('{0}.dist'.format(prefix), sep='\t', index_col=0, header=0)
    ref_dist_mat.columns = ['strain', 'distance']

    logger.info('Going through each time bin to run pairwise MASH and aggregate results')
    while end_date <= max(metadata['date']):
        sub_meta = metadata.query('date >= @start_date and date <= @end_date')
        ids = sub_meta['strain']

        prefix = '{0}/{1}_{2}'.format(tmpdir, start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'))

        # Create fasta file of isolates within time bin
        for strain in ids:
            subprocess.run(['grep -A1 {0} {1} >> {2}.fasta'.format(strain, args.fasta, prefix)], shell=True)

        # Create sketch 
        subprocess.run(['mash sketch -k {0} -p {1} -o {2}.msh -i {3}.fasta'.format(args.kmer, args.threads, prefix, prefix)], shell=True, stdout=subprocess.DEVNULL)

        # Run MASH distance calculation
        subprocess.run(['mash dist -t {0}.msh {1}.msh >> {2}.dist'.format(prefix, prefix, prefix)], shell=True, stdout=subprocess.DEVNULL)

        # Load distance matrix
        pair_dist_mat = pd.read_csv('{0}.dist'.format(prefix), sep='\t', index_col=0, header=0)
        pair_distances = []

        # Parse distance matrix (might be a better way to do this)
        for i in ids:
            for j in ids:
                if i == j:
                    break
                pair_distances.append(pair_dist_mat.loc[i, j])

        logger.info("Current bin has {0} isolates, {1} pairwise comparisons".format(len(ids), len(pair_distances)))

        sub_ref_dist = ref_dist_mat[dist_mat['strain'].isin(ids)]
        out_ref_dist.append(sub_ref_dist['distance'].mean())
        out_ref_max.append(sub_ref_dist['distance'].max())
        out_ref_min.append(sub_ref_dist['distance'].min())

        bin_count = len(ids)
        out_dates.append('{0} - {1}'.format(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
        out_pair_dist.append(statistics.mean(pair_distances))
        out_pair_min.append(min(pair_distances))
        out_pair_max.append(max(pair_distances))
        bin_size.append(bin_count)

        start_date += pd.DateOffset(days=args.window)
        end_date += pd.DateOffset(days=args.window)

    out_df = pd.DataFrame({'date_range': out_dates, 
                           'average_ref_distance': out_ref_dist, 
                           'min_ref_distance': out_ref_min, 
                           'max_ref_distance': out_ref_max,
                           'average_pair_distance': out_pair_dist,
                           'min_pair_distance': out_pair_min,
                           'max_pair_distance': out_pair_max, 
                           'sample_size': bin_size})
    out_df.to_csv(os.path.join(args.outdir, '{0}.tsv'.format(args.prefix)), sep='\t', index=False)

    if args.clean and tmpdir_obj is not None:
        tmpdir_obj.cleanup()
