#!/usr/bin/env python

import argparse as ap
import pandas as pd
import subprocess
import os
import tempfile
import statistics

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Converts a distance matrix into time bins')

    subparser = parser.add_subparsers(dest='subcommand', required=True, help='Run pairwise MASH or against reference')

    ref_parser = subparser.add_parser('reference', help='Run MASH against a reference')
    ref_parser.add_argument('ref', metavar='REF', help='Reference genome in FASTA format')

    pair_parser = subparser.add_parser('pairwise', help='Run pairwise MASH')

    for sub in [ref_parser, pair_parser]:   
        sub.add_argument('-t', '--time', type=int, default=14, help='Number of days within a bin (default: 14)')
        sub.add_argument('-w', '--window', type=int, default=1, help='Number of days between bins (default: 1)')
        sub.add_argument('-c', '--clean', action='store_true', help='Remove all temporary MASH files after execution')
        sub.add_argument('-k', '--kmer', type=int, default=17, help='kmer size for MASH (default: 17)')
        sub.add_argument('-p', '--threads', type=int, default=1, help='Number of threads to use for MASH (default: 1)')
        sub.add_argument('--prefix', type=str, default='summary', help='Output file prefix (default: summary)')
        sub.add_argument('meta', metavar='METADATA', help='GISAID metadata')
        sub.add_argument('fasta', metavar='FASTA', help='FASTA of sequences')
        sub.add_argument('outdir', metavar='DIR', help='Output directory')

    args = parser.parse_args()

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
    metadata = pd.read_csv(args.meta, sep='\t')
    metadata['date'] = pd.to_datetime(metadata['date'])
    metadata = metadata.sort_values('date')
    metadata = metadata.reset_index()

    start_date = metadata['date'][0]
    end_date = start_date + pd.DateOffset(days=args.time)

    out_dates = []
    out_dist = []
    out_min = []
    out_max = []
    bin_size = []

    if args.subcommand == 'reference':
        prefix = os.path.join(tmpdir, args.prefix)

        # Create sketch 
        subprocess.run(['mash sketch -k {0} -p {1} -o {2}.msh -i {3}'.format(args.kmer, args.threads, prefix, args.fasta)], shell=True)

        # Run MASH distance calculation
        subprocess.run(['mash dist -t {0}.msh {1} >> {2}.dist'.format(prefix, args.ref, prefix)], shell=True)

        dist_mat = pd.read_csv('{0}.dist'.format(prefix), sep='\t', index_col=0, header=0)
        dist_mat.columns = ['strain', 'distance']

        while end_date <= max(metadata['date']):
            sub_meta = metadata.query('date >= @start_date and date <= @end_date')
            ids = sub_meta['strain']
            sub_dist = dist_mat[dist_mat['strain'].isin(ids)]
            avg_dist = sub_dist['distance'].mean()
            min_dist = sub_dist['distance'].min()
            max_dist = sub_dist['distance'].max()
            bin_count = sub_dist['distance'].count()
            out_dates.append('{0} - {1}'.format(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
            out_dist.append(avg_dist)
            out_min.append(min_dist)
            out_max.append(max_dist)
            bin_size.append(bin_count)

        start_date += pd.DateOffset(days=args.window)
        end_date += pd.DateOffset(days=args.window)

    else:
        while end_date <= max(metadata['date']):
            sub_meta = metadata.query('date >= @start_date and date <= @end_date')
            ids = sub_meta['strain']

            prefix = '{0}/{1}_{2}'.format(tmpdir, start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'))

            # Create fasta file of isolates within time bin
            for strain in ids:
                subprocess.run(['grep -A1 {0} {1} >> {2}.fasta'.format(strain, args.fasta, prefix)], shell=True)

            # Create sketch 
            subprocess.run(['mash sketch -k {0} -p {1} -o {2}.msh -i {3}.fasta'.format(args.kmer, args.threads, prefix, prefix)], shell=True)

            # Run MASH distance calculation
            subprocess.run(['mash dist -t {0}.msh {1}.msh >> {2}.dist'.format(prefix, prefix, prefix)], shell=True)

            # Load distance matrix
            dist_mat = pd.read_csv('{0}.dist'.format(prefix), sep='\t', index_col=0, header=0)
            distances = []

            # Parse distance matrix (might be a better way to do this)
            for i in ids:
                for j in ids:
                    if i == j:
                        break
                    distances.append(dist_mat.loc[i, j])

            avg_dist = statistics.mean(distances)
            min_dist = min(distances)
            max_dist = max(distances)
            bin_count = len(distances)
            out_dates.append('{0} - {1}'.format(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
            out_dist.append(avg_dist)
            out_min.append(min_dist)
            out_max.append(max_dist)
            bin_size.append(bin_count)

            start_date += pd.DateOffset(days=args.window)
            end_date += pd.DateOffset(days=args.window)

    out_df = pd.DataFrame({'date_range': out_dates, 'average_distance': out_dist, 'min_distance': out_min, 'max_distance': out_max, 'sample_size': bin_size})
    out_df.to_csv(args.out, sep='\t', index=False)

    if args.clean and tmpdir_obj is not None:
        tmpdir_obj.cleanup()