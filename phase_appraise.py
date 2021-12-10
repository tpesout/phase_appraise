#!/usr/bin/env python3
from __future__ import print_function
import argparse
import glob
import sys
import matplotlib.pyplot as plt
import pysam
import collections
import numpy as np
import os

OUTPUT_LOG_FILE_NAME = None

HP_TAG = "HP"
UNCLASSIFIED = 'u'
IN_CIS = 'c'
IN_TRANS = 't'
UNKNOWN = 'k'

def parse_args(args = None):
    parser = argparse.ArgumentParser("Compares phasing for haplotagged reads.  Calculates local accuracy and switch errors")
    parser.add_argument('--input_bam', '-i', dest='input_bam', default=None, required=True, type=str,
                        help='Haplotagged BAM file ')
    parser.add_argument('--truth_hap1', '-1', dest='truth_hap1', default=None, required=True, type=str,
                       help='Truth Hap1 readset')
    parser.add_argument('--truth_hap2', '-2', dest='truth_hap2', default=None, required=True, type=str,
                       help='Truth Hap2 readset')
    parser.add_argument('--region', '-r', dest='region', default=None, required=False, type=str,
                       help='Region to run on')
    parser.add_argument('--phaseset_bed', '-p', dest='phaseset_bed', default=None, required=False, type=str,
                       help='File describing phase sets')
    parser.add_argument('--high_conf_bed', '-C', dest='high_conf_bed', default=None, required=False, type=str,
                       help='If set, will only count regions within this BED file')
    parser.add_argument('--title', '-t', dest='title', default=None, required=False, type=str,
                       help='Figure title')
    parser.add_argument('--spacing', '-s', dest='spacing', default=1000, required=False, type=int,
                       help='What depth should be used on y axes')
    parser.add_argument('--output_name', '-o', dest='output_name', default=None, required=False, type=str,
                       help='Output name (will save if set)')
    parser.add_argument('--output_name_bam', '-O', dest='output_name_bam', default=False, required=False, action='store_true',
                       help='Output name should be based off of BAM file (-o overrides this)')
    parser.add_argument('--max_depth', '-d', dest='max_depth', default=72, required=False, type=int,
                       help='What depth should be used on y axes')
    parser.add_argument('--only_natural_switch', '-N', dest='only_natural_switch', default=False, required=False, action='store_true',
                       help='Only plot the natural switch')
    parser.add_argument('--global_data', '-G', dest='global_data', default=False, required=False, action='store_true',
                        help='Output tsv file containing contig, start, end, cis, trans, unk, unc, correct ratio and switch for every bucket ')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)
    global OUTPUT_LOG_FILE_NAME
    if OUTPUT_LOG_FILE_NAME is not None:
        with open(OUTPUT_LOG_FILE_NAME, 'a') as out:
            print(msg, file=out)

def write_global_data(msg, fig_name_base):
    '''
    Writes a list of strings (msg) as tab-separated values global data tsv file
    '''
    filename= fig_name_base+"_global_data.tsv"
    with open(filename, 'a') as out:
        print(*msg, file=out, sep='\t', end='\n')

def smooth_values(values,size=5,should_sum=False):
    new_values = []
    for i in range(0,len(values)):
        s=max(0,i-size)
        e=min(len(values),i+size)
        new_values.append(sum(values[s:e]) if should_sum else np.mean(values[s:e]))
    return new_values


def plot_only_natural_switch(classification_data, args, phasesets=None, fig_name=None):

    start_idx = min(classification_data.keys())
    end_idx = max(classification_data.keys())
    x = []
    right = []
    rong = []
    for i in range(start_idx, end_idx + 1):
        x.append(i / (1000000.0 / args.spacing))
        ri = classification_data[i][IN_CIS]
        ro = classification_data[i][IN_TRANS]

        right.append(min(args.max_depth, ri))
        rong.append(-1 * min(args.max_depth, ro))

    # get plots
    # fig, (ax1) = plt.subplots(nrows=1,ncols=1)
    fig, (ax1) = plt.subplots(nrows=1,ncols=1, figsize=(8, 1.75))
    lw = .5

    ax1.set_ylabel('Concordant/Discordant\nDepth')
    ax1.set_xlabel('Chr1 Positions (Mb)')
    ax1.fill_between(x, right, color='tab:red', linewidth=lw)
    ax1.fill_between(x, rong, color='tab:blue', linewidth=lw)
    ax1.set_ylim(-1 * args.max_depth, args.max_depth)

    if phasesets is not None:
        top = True
        for ps in phasesets:
            modifier = 1000000.0
            start = ps[0] / modifier
            end = ps[1] / modifier
            # ax1.plot(range(start, end), [2 if top else -2 for _ in range(start, end)],
            #          color='black', alpha=.65, linewidth=2)
            # start = ps[0] // SPACING
            # end = ps[1] // SPACING
            ps_range = list(map(lambda x: x/modifier, range(int(start*modifier), int(end*modifier), 10000)))
            ax1.plot(ps_range, [2 if top else -2 for _ in ps_range], color='black', alpha=1, linewidth=.75)
            top = not top

    fig.tight_layout()
    # fig.set_size_inches(12, 3)
    if fig_name is not None:
        plt.savefig(fig_name + ".png", format='png', dpi=50)
        plt.savefig(fig_name + ".pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()



def plot_full(classification_data, args, fig_name_base, fig_name=None, phasesets=None, title=None, highconf_positions=None, region=None):
    log("Generating plots for region {}".format("all" if region is None else region))
    start_idx = min(classification_data.keys())
    end_idx = max(classification_data.keys())
    x = []
    in_cis = []
    in_trans = []
    total_classified = []
    total_reads = []
    unknown = []
    unclassified = []
    correct_ratio = []
    switches = []
    prev_val = None
    for i in range(start_idx, end_idx + 1):
        x.append(i)
        cis = classification_data[i][IN_CIS]
        trans = classification_data[i][IN_TRANS]
        unc = classification_data[i][UNCLASSIFIED]
        unk = classification_data[i][UNKNOWN]
        if highconf_positions is not None and i not in highconf_positions:
            cis = 0
            trans = 0
            unc = 0
            unk = 0
        total = cis+trans+unc+unk

        in_cis.append(cis)
        in_trans.append(-1 * trans)
        unknown.append(unc)
        unclassified.append(unk)
        total_classified.append(cis + trans + unc)
        total_reads.append(total)
        curr_correct_ratio = None if cis + trans == 0 else 100 * abs(max(cis, trans) / (cis + trans))
        correct_ratio.append(curr_correct_ratio)
        curr_correct_ratio_abs = 0.0 if cis + trans + unc == 0 else abs(max(cis, trans) / (cis + trans + unc))

        val = IN_CIS if cis > trans else IN_TRANS if trans > cis else None
        switch = False
        if val is not None:
            if prev_val is not None and val != prev_val:
                switches.append(i)
                switch = True
            prev_val = val

        # Write record to global data file if requested
        if args.global_data:
            contig = "unknown" if region is None else region
            start = i * args.spacing
            end = start + args.spacing
            record = [contig, start, end, cis, trans, unk, unc, 0 if curr_correct_ratio is None else curr_correct_ratio,curr_correct_ratio_abs ,switch]
            write_global_data(record,fig_name_base)

    # get averages
    avg_unknown = np.mean(list(map(lambda y: y[1], filter(lambda x: x[0] != 0, zip(total_reads, unknown)))))
    avg_classified = np.mean(list(map(lambda y: y[1], filter(lambda x: x[0] != 0, zip(total_reads, total_classified)))))
    avg_total_reads = np.mean(list(filter(lambda x: x != 0, total_reads)))
    avg_correct = np.mean(list(filter(lambda x: x is not None, correct_ratio)))
    std_correct = np.std(list(filter(lambda x: x is not None, correct_ratio)))

    # smooth values
    smoothed_unknown = smooth_values(unknown, size=20)
    smoothed_total_classified = smooth_values(total_classified, size=20)
    smoothed_total_reads = smooth_values(total_reads, size=20)

    # switches
    switch_count = len(switches)
    total_evaluated_buckets = len(list(filter(lambda x: x is not None, correct_ratio)))
    switch_error_rate = 1.0 * switch_count / total_evaluated_buckets

    # get plots
    fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=3,ncols=1,sharex='all',gridspec_kw={'height_ratios': [2, 1, 1]})
    lw = .5

    ax1.set_ylabel('Phasing Partitions')
    ax1.fill_between(x, in_cis, color='tab:red', linewidth=lw)
    ax1.fill_between(x, in_trans, color='tab:blue', linewidth=lw)
    ax1.set_ylim(-1 * args.max_depth, args.max_depth)
    ax1.scatter(switches, [0 for _ in switches], marker="*", color='black')
    log("Switches:\n\tTotal Switches: {}\n\tTotal Buckets: {}\n\tSwitch Rate: {}".format(switch_count, total_evaluated_buckets, switch_error_rate))
    ax1.annotate("Switch Rate: {:6.4f}".format(switch_error_rate), (x[0], args.max_depth * .9), fontfamily='monospace', fontsize=12,weight="bold")

    ax2.set_ylabel('Correct Ratio')
    c = ax2.scatter(x, correct_ratio, c=correct_ratio, s=lw, cmap=plt.cm.get_cmap('RdYlGn'), alpha=.5)
    # c = ax2.scatter(x, [102.5 if c is not None else 0 for c in correct_ratio], c=correct_ratio, s=lw, cmap=plt.cm.get_cmap('RdYlGn'), alpha=.5)
    ax2.plot(x, [avg_correct for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Correct ratio:\n\tAvg: {}\n\tStd: {}".format(avg_correct, std_correct))
    ax2.annotate("Average Accuracy: {:5.2f}".format(avg_correct), (x[0], avg_correct+1), fontfamily='monospace', fontsize=12,weight="bold")
    ax2.set_ylim(45, 105)

    if phasesets is not None:
        concordancy_in_phasesets = list()
        top = True
        for ps in phasesets:
            start = ps[0] // args.spacing
            end = ps[1] // args.spacing
            if not any([highconf_positions is None or x in highconf_positions for x in range(start, end) ]):
                continue
            ax2.plot(range(start, end), [48 if top else 46 for _ in range(start, end)],
                     color='black', alpha=.65, linewidth=2)
            ax1.plot(range(start, end), [2 if top else -2 for _ in range(start, end)],
                     color='black', alpha=.65, linewidth=2)
            if start != end:
                top = not top
                concordancy_in_phasesets.extend(list(filter(lambda x: x is not None, correct_ratio[start:end])))
        avg_ps_correct = np.mean(concordancy_in_phasesets)
        log("Correct ratio in phasesets:\n\tAvg: {}\n\tStd: {}".format(avg_ps_correct,
                                                                       np.std(concordancy_in_phasesets)))

    ax3.set_ylabel('Classified Depth')
    ax3.plot(x, smoothed_unknown, color='red', alpha=.25, linewidth=lw)
    ax3.plot(x, smoothed_total_classified, color='blue', alpha=.25, linewidth=lw)
    ax3.plot(x, smoothed_total_reads, color='black', alpha=.25, linewidth=lw)

    ax3.plot(x, [avg_unknown for _ in x], color='red', alpha=.5, linewidth=lw)
    log("Total Unknown:\n\tAvg: {}\n\tStd: {}".format(avg_unknown, np.std(unknown)))
    ax3.annotate("{:3d}x Unknown".format(int(avg_unknown)), (x[0], avg_unknown+1), fontfamily='monospace', color='darkred', fontsize=12, weight="bold")

    ax3.plot(x, [avg_classified for _ in x], color='blue', alpha=.5, linewidth=lw)
    log("Total Classified:\n\tAvg: {}\n\tStd: {}".format(avg_classified, np.std(total_classified)))
    ax3.annotate("{:3d}x Classified".format(int(avg_classified)), (x[0], avg_classified+1), fontfamily='monospace', color='darkblue', fontsize=12, weight="bold")

    ax3.plot(x, [avg_total_reads for _ in x], color='black', alpha=.5, linewidth=lw)
    log("Total Reads:\n\tAvg: {}\n\tStd: {}".format(avg_total_reads, np.std(total_reads)))
    ax3.annotate("{:3d}x Total Reads".format(int(avg_total_reads)), (x[0], max(avg_total_reads, avg_classified+8)+1), fontfamily='monospace', color='black', fontsize=12,weight="bold")

    ax3.set_ylim(-.05 * args.max_depth, args.max_depth)

    if title is not None:
        ax1.set_title(title)
    fig.tight_layout()
    fig.set_size_inches(18, 12)
    if fig_name is not None:
        plt.savefig(fig_name+".png", format='png', dpi=50)
        plt.savefig(fig_name+".pdf", format='pdf', dpi=300)
    plt.show()
    plt.close()


def read_phaseset_bed(bed_file):
    phasesets = list()
    with open(bed_file) as bed:
        for line in bed:
            parts = line.split("\t")
            assert(len(parts) >= 3)
            phasesets.append((int(parts[1]), int(parts[2])))
    return phasesets


def get_highconf_positions(bed_file, args):
    highconf_positions = set()
    with open(bed_file) as bed:
        for line in bed:
            parts = line.split("\t")
            assert(len(parts) >= 3)
            for i in range(int(int(parts[1]) / args.spacing), int(int(parts[2]) / args.spacing) + 1):
                highconf_positions.add(i)
    return highconf_positions


def get_position_classifications(bam_location, truth_h1_ids, truth_h2_ids, args, region=None, verbose=True):
    # get read phasing pairs
    log("Reading reads for region {}".format("all" if region is None else region))
    samfile = None
    read_count = 0
    missing_hp_count = 0
    position_classifications = collections.defaultdict(
        lambda : collections.defaultdict(lambda : 0)
    )
    analyzed_lengths = []
    try:
        samfile = pysam.AlignmentFile(bam_location, 'rb' if bam_location.endswith("bam") else 'r')
        for read in samfile.fetch(region=region):
            read_count += 1
            id = read.query_name
            spos = read.reference_start
            epos = read.reference_end

            if not read.has_tag(HP_TAG):
                missing_hp_count += 1
                classifier = UNCLASSIFIED
            else:
                hp = read.get_tag(HP_TAG)
                if hp == 0:
                    classifier = UNCLASSIFIED
                elif id not in truth_h1_ids and id not in truth_h2_ids:
                    classifier = UNKNOWN
                elif hp == 1 and id in truth_h1_ids:
                    classifier = IN_CIS
                elif hp == 2 and id in truth_h2_ids:
                    classifier = IN_CIS
                else:
                    classifier = IN_TRANS

            if classifier != UNCLASSIFIED:
                analyzed_lengths.append(epos - spos)

            while spos <= epos:
                pos = int(spos / args.spacing)
                position_classifications[pos][classifier] += 1
                spos += args.spacing


    finally:
        if samfile is not None: samfile.close()

    if verbose:
        log("Classified Read Lengths{}:".format("" if region is None else " for {}".format(region)))
        log("\tmean:   {}".format(np.mean(analyzed_lengths)))
        log("\tmedian: {}".format(np.median(analyzed_lengths)))
        analyzed_lengths.sort()
        len_total = sum(analyzed_lengths)
        len_curr = 0
        for l in analyzed_lengths:
            len_curr += l
            if len_curr > len_total/2:
                log("\tN50:    {}".format(l))
                break

    return position_classifications


def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # # for "professional" plots
    # if args.output_name is not None or args.output_name_bam is not None:
    #     plt.style.use('ggplot')
    #     plt.rcParams.update({'font.size': 8})
    #     plt.rcParams['pdf.fonttype'] = 42
    #     plt.switch_backend('agg')

    # get figure name
    fig_name_base = args.output_name
    if fig_name_base is None and args.output_name_bam:
        fig_name_base = os.path.basename(args.input_bam) + (
            ".natural_switch" if args.only_natural_switch else ("" if args.high_conf_bed is None else ".highconf"))

    # update logging
    if fig_name_base is not None:
        global OUTPUT_LOG_FILE_NAME
        OUTPUT_LOG_FILE_NAME = fig_name_base + ("" if args.region is None else "." + args.region.replace(":", "-")) + ".log"
        with open(OUTPUT_LOG_FILE_NAME, 'w') as out:
            pass

    # Write header of global data file if option is provided
    if args.global_data:
        write_global_data(["##input=".format(args.input_bam)], fig_name_base)
        write_global_data(["##hap1=".format(args.truth_hap1)], fig_name_base)
        write_global_data(["##hap2=".format(args.truth_hap2)], fig_name_base)
        write_global_data(["#contig", "start", "end", "cis", "trans", "unknown", "unclassified", "correct_ratio", "correct_ratio_abs","switch"], fig_name_base)

    # get truth reads
    truth_h1 = set()
    truth_h2 = set()
    with open(args.truth_hap1, 'r') as fin:
        for line in fin:
            truth_h1.add(line.split(',')[0].strip())
    with open(args.truth_hap2, 'r') as fin:
        for line in fin:
            truth_h2.add(line.split(',')[0].strip())
    log("Found {} truth H1 reads and {} truth H2 reads".format(len(truth_h1), len(truth_h2)))

    # classify positions for reads
    position_classifications = dict()

    # case of user selecting only one contig
    if args.region is not None:
        position_classifications[args.region.replace(":", "-")] = get_position_classifications(args.input_bam, truth_h1, truth_h2, args, args.region)
    # case of using all contigs in bam
    else:
        # find all contigs in bam
        samfile = pysam.AlignmentFile(args.input_bam, 'rb' if args.input_bam.endswith("bam") else 'r')
        contigs=samfile.references
        samfile.close()
        # get position classifications for every contig - triply nested dict
        for region in contigs:
            position_classifications[region] = get_position_classifications(args.input_bam, truth_h1, truth_h2, args)

    #TODO needs to be contig/region aware
    highconf_positions = None
    if args.high_conf_bed is not None:
        highconf_positions = get_highconf_positions(args.high_conf_bed, args)

    #TODO needs to be contig/region aware
    phasesets = None
    if args.phaseset_bed is not None:
        phasesets = read_phaseset_bed(args.phaseset_bed)

    for region in sorted(list(position_classifications.keys())):
        # get region data
        region_position_classifications = position_classifications[region]
        # TODO get after making contig/region aware
        region_phasesets = None
        region_highconf = None

        # get title
        if args.title is None:
            title = os.path.basename(args.input_bam)
            if region is not None:
                title = "{}.{}".format(title, region)
        else:
            title = args.title

        # get fig name
        if region is None:
            fig_name = fig_name_base
        else:
            fig_name = "{}.{}".format(fig_name_base, region)

        # for only natural switch
        if (args.only_natural_switch):
            plot_only_natural_switch(position_classifications, args, phasesets, fig_name)
            continue

        # plot it
        plot_full(region_position_classifications, args, fig_name_base, fig_name=fig_name, phasesets=region_phasesets,
                title=title, highconf_positions=region_highconf, region=region)

if __name__ == "__main__":
    main()











