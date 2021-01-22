#!/usr/bin/env python

import argparse as ap
import pandas as pd

if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Converts a distance matrix into time bins')

    parser.add_argument('dist', metavar='MATRIX', help='MASH distance matrix (against a reference)')
    parser.add_argument('meta', metavar='METADATA', help='GISAID metadata')
    parser.add_argument('out', metavar='OUTPUT', help='Output file')
    parser.add_argument('-t', '--time', type=int, default=14, help='Number of days within a bin (default: 14)')

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

    while end_date <= max(metadata['date']):
        sub_meta = metadata.query('date >= @start_date and date <= @end_date')
        ids = sub_meta['strain']
        sub_dist = dist_mat[dist_mat['strain'].isin(ids)]
        avg_dist = sub_dist['distance'].mean()
        min_dist = sub_dist['distance'].min()
        max_dist = sub_dist['distance'].max()
        out_dates.append('{0} - {1}'.format(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))
        out_dist.append(avg_dist)
        out_min.append(min_dist)
        out_max.append(max_dist)

        start_date += pd.DateOffset(days=1)
        end_date += pd.DateOffset(days=1)

    out_df = pd.DataFrame({'date_range': out_dates, 'average_distance': out_dist, 'min_distance': out_min, 'max_distance': out_max})
    out_df.to_csv(args.out, sep='\t', index=False)
