#!/usr/bin/env python

import argparse as ap
import pandas as pd

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Converts a distance matrix into time bins')

    parser.add_argument('dist', metavar='MATRIX', help='MASH distance matrix (against a reference)')
    parser.add_argument('meta', metavar='METADATA', help='GISAID metadata')
    parser.add_argument('out', metavar='OUTPUT', help='Output file')
    parser.add_argument('-t', '--time', type=int, default=14, help='Number of days within a bin (default: 14)')
    parser.add_argument('-w', '--window', type=int, default=1, help='Number of days between bins (default: 1)')

    args = parser.parse_args()

    dist_mat = pd.read_csv(args.dist, sep='\t')
    dist_mat.columns = ['strain', 'distance']

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

    while end_date <= max(metadata['date']):
        sub_meta = metadata.query('date >= @start_date and date <= @end_date')
        ids = sub_meta['strain']
        sub_dist = dist_mat[dist_mat['strain'].isin(ids)]
        avg_dist = sub_dist['distance'].mean()
        min_dist = sub_dist['distance'].min()
        max_dist = sub_dist['distance'].max()
        bin_size = sub_dist['distance'].count()
        out_dates.append('{0} - {1}'.format(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
        out_dist.append(avg_dist)
        out_min.append(min_dist)
        out_max.append(max_dist)

        start_date += pd.DateOffset(days=args.window)
        end_date += pd.DateOffset(days=args.window)

    out_df = pd.DataFrame({'date_range': out_dates, 'average_distance': out_dist, 'min_distance': out_min, 'max_distance': out_max, 'sample_size': bin_size})
    out_df.to_csv(args.out, sep='\t', index=False)
