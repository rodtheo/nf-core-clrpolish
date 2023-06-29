#!/usr/bin/env python

import json
import argparse
import logging
import sys
import csv
from pathlib import Path
import re
from collections import defaultdict

logger = logging.getLogger()

base_dict_multiqc_line = {
    "id": "custom_data_lineplot",
    "section_name": "Polishing Results",
    "description": "Graphed are the consensus quality (QV*) obtained after running merfin -hist (yaxis).",
    "plot_type": "linegraph",
    "pconfig": {
        "id": "custom_data_linegraph",
        "title": "Polishing Rounds Evaluation",
        "ylab": "Consensus quality (QV*)",
        "xDecimals": False,
        # "yPlotLines": [
        #     {
        #         "color": "#FF0000",
        #         "value": 10
        #     }
        # ],
    },
    # "data": {
    #     "sample_1": {
    #         "1": 12,
    #         "2": 14,
    #         "3": 10,
    #         "4": 7,
    #         "5": 16
    #     },
    #     "sample_2": {
    #         "1": 9,
    #         "2": 11,
    #         "3": 15,
    #         "4": 18,
    #         "5": 21
    #     }
    # }
}

result_table = """# plot_type: 'table'
# section_name: 'Table - Summary Results'
# description: 'Main metrics that evaluates genome polish process'
# pconfig:
#     namespace: 'Cust Data'
# headers:
#     Polish_Round:
#         title: 'Polish_Round'
#         description: 'Number of iteration applied. Round 0 means the unpolished assembly.'
#     QV*:
#         title: 'QV*'
#         description: 'QV* of the assembly as obtained by running merfin -hist'
#         format: '{:,.2f}'
#     K*:
#         title: 'K*'
#         description: 'K* completeness obtained by running merfin -completeness'
#         format: '{:,.4f}'
#     FASTA:
#         title: 'FASTA'
#         description: 'Fasta name of the evaluated assembly'
"""



def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "sample_info",
        metavar="SAMPLE_INFO",
        type=str,
        help="List of lists input format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    my_lists = json.loads(args.sample_info.encode('utf8').decode('utf8'), parse_float=True)
    my_dict_out = defaultdict(lambda: {})
    for rec in my_lists:
        genome_id_re = re.search(r"((round_\d+_)*)(.*)$", rec[2])
        genome_id = genome_id_re.groups()[2]
        my_dict_out[genome_id][str(int(rec[0]) - 1)] = float(rec[1])
    base_dict_multiqc_line['data'] = my_dict_out
    with args.file_out.open('w') as outf: 
        json.dump(base_dict_multiqc_line, outf)
    with open('table_results_mqc.txt', 'wt') as out_file_txt:
        out_file_txt.write(result_table)
        tsv_writer = csv.writer(out_file_txt, delimiter='\t')
        header = ["Polish_Round", "QV*", "K*", "FASTA"]
        tsv_writer.writerow(header)
        for rec in my_lists:
            tsv_writer.writerow([str(int(rec[0]) - 1), rec[1], rec[3], rec[2]])
        # tsv_writer.writerows(my_lists)
	


if __name__ == "__main__":
    sys.exit(main())