#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import csv
import logging
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filters records from GISAID Metadata '
                    'and Sequences')
    parser.add_argument('-ll', '--loglevel', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                 'CRITICAL'],
                        help='Set the logging level')
    parser.add_argument('--table', type=str, default=None,
                        help='GISAID Metadata file (.tsv) format')
    parser.add_argument('--fasta', type=str, default=None,
                        help='GISAID Sequence file in .fasta or .fa '
                             'format')
    parser.add_argument('--outdir', type=str, default='./',
                        help='Output directpry')
    parser.add_argument('--location', type=str, default=None,
                        help='location of interest from GISAID e.g. '
                             'region:North America, country:Canada or '
                             'division:England etc')
    parser.add_argument('--startdate', type=str, default='2020-01-10',
                        help='Start date yyyy-mm-dd')
    parser.add_argument('--enddate', type=str, default=None,
                        help='End date yyyy-mm-dd')
    return parser.parse_args()


def filter_fasta(id_list, loc_id, fastain, out_path):
    """ Filter fasta file """
    fin = open(fastain, 'r')
    fout = open(out_path + loc_id + "_" + str(args.startdate) +
                "_" + str(args.enddate) +
                "_GISAID.fasta", 'w')

    for record in SeqIO.parse(fin, 'fasta'):
        for item in id_list:
            if item.strip() == record.id:
                fout.write(">" + record.id + "\n")
                fout.write(str(record.seq) + "\n")

    fin.close()
    fout.close()


def filter_metadata(df, startdate, enddate, loctype, loc):
    """ Filter Metadata file
        host = human
        length >=29kb
        startdate <= date >= enddate
        location = region/country/division
     """
    df1 = df[(df['host'].str.lower() == 'human') & (
            df[loctype].str.lower() == loc.lower()) & (
                     df['length'] >= 29000)]
    sdate = pd.to_datetime(startdate).date()
    edate = pd.to_datetime(enddate).date()
    df2 = df1[df1["date"].isin(pd.date_range(sdate, edate))]

    return df2


if __name__ == '__main__':
    args = parse_args()
    logging.basicConfig(level=args.loglevel,
                        format='%(asctime)s (%(relativeCreated)d ms) '
                               '-> %(levelname)s: %(message)s',
                        datefmt='%I:%M:%S %p')

    logger = logging.getLogger()

    if args.location is None:
        logger.error(
            'No Region, Country or division provided to subset')
        sys.exit(0)

    if args.enddate is None:
        logger.error('No end date provided')
        sys.exit(0)

    Metadata = pd.read_csv(args.table, sep="\t", low_memory=False,
                           parse_dates=['date'])
    Metadata['date'] = pd.to_datetime(Metadata.date, format='%Y-%m-%d',
                                      errors='coerce')

    if args.location.split(":")[0] == 'region':
        loc = 'region'
    elif args.location.split(":")[0] == 'country':
        loc = 'country'
    else:
        loc = 'division'

    df = filter_metadata(df=Metadata, startdate=args.startdate,
                         enddate=args.enddate, loctype=loc,
                         loc=str(args.location.split(":")[1]))

    df.to_csv(
        args.outdir + str(args.location.split(":")[1]).replace(" ",
                                                               "") + "_" + str(
            args.startdate) + "_" + str(args.enddate) + "_GISAID.tsv",
        sep="\t", quoting=csv.QUOTE_NONE, index=False, header=True)

    filter_fasta(id_list=df['strain'].tolist(),
                 loc_id=str(args.location.split(":")[1]).replace(" ",
                                                                 ""),
                 fastain=args.fasta, out_path=args.outdir)
