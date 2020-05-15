from collections import defaultdict

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def validate_name(df, gene_name):

    gene_df = df.loc[df['feature'] == 'gene']
    d_dict = defaultdict(list)
    for gene, id in zip(gene_df.gene_name, gene_df.gene_id):
        d_dict[gene].append(id)

    if(len(d_dict[gene_name]) != 1):
        print('Error - {} contains more than one corresponding id ({})'.format(
            gene_name, ','.join(d_dict[gene_name])))
        exit()

    else:
        return d_dict[gene_name][0]


def get_snos(gene_id):
    file = '/home/danx/Documents/projects/snoRNA_network_V2/ref/modified/sno_host.tsv'
    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.host_id == gene_id]
    return list(df.sno_id)


def load_ref(gene_name, df):

    gene_id = validate_name(df, gene_name)

    snos_list = get_snos(gene_id)

    sno_file = '/home/danx/Documents/projects/snoRNA_network_V2/ref/modified/human_ensembl_87_wo_dup_v3.BB_v3_PARSED_snoRNAs.bed'
    snos = pd.read_csv(sno_file, sep='\t')

    snos = snos[['gene_name', 'start', 'end']]
    # print(snos)

    df = df.loc[df.gene_id == gene_id]

    # print(df.columns)
    df = df[['seqname', 'feature', 'start', 'end', 'strand', 'exon_id',
             'exon_number', 'gene_biotype', 'gene_name', 'gene_id',
             'transcript_id', 'transcript_name', 'transcript_biotype', 'transcript_support_level']]
    df = df.loc[(df['feature'] == 'gene') | (df['transcript_support_level'] < 6)]
    df = df.loc[df.feature.isin(['gene', 'exon', 'transcript'])]

    # df.to_csv("test", sep='\t', index=False)

    return df, snos


def draw_sno_in_host(ax, start_, t_start, t_end,
                     sno_in_host, yloc, BOX_HEIGHT, high):
    snos = []
    for i in sno_in_host.index:
        w = sno_in_host.at[i, 'end'] - sno_in_host.at[i, 'start']
        s = sno_in_host.at[i, 'start'] - start_
        name = sno_in_host.at[i, 'gene_name'].replace('SNOR', '')
        # print('{}: {} {}'.format(sno_in_host.at[i, 'gene_name'], s, w+s))

        if s >= t_start and w + s <= t_end:
            rect = Rectangle((s, yloc-BOX_HEIGHT/4), w, BOX_HEIGHT/2)
            snos.append(rect)
            pc = PatchCollection(snos, facecolor='grey', alpha=1,
                                 edgecolor=None, zorder=100)
            ax.add_collection(pc)
            ax.text(s,
                    yloc-(BOX_HEIGHT * (0.2714 + high * 0.0342857)),  # (high / 30)
                    name,
                    fontname='Comic sans ms',
                    fontsize=6)


def prepare_fig(ref_df_, sno_in_host, gathered_trans):

    FONTSIZE = 10
    BOX_HEIGHT = 0.8
    STEPS = 1.5

    ref_df = ref_df_.copy(deep=True)
    ref_df.reset_index(inplace=True)
    gene = ref_df.loc[ref_df.feature == 'gene']
    start_ = gene.at[0, 'start']
    start = 0
    end = gene.at[0, 'end'] - start_
    gene_name = gene.at[0, 'gene_name']

    # Create figure and axes
    fig, axes = plt.subplots(2)

    offset = 0.05 * end

    titles = ['NMD transcripts', 'tsl1 transcript (protein_coding)']

    for i, b in enumerate([True, False]):

        transcripts = sorted(list({x for x in ref_df.transcript_id if str(x) != 'nan'}))

        if b:
            # ADDED for NMD viz ------------------------------------------------------------------------------
            transcripts = sorted([x for x in transcripts if x in gathered_trans])[::-1]
            # ADDED for NMD viz ------------------------------------------------------------------------------
        else:
            # print(ref_df.columns)
            tsl1_df = ref_df.loc[(ref_df['transcript_support_level'] < 2)
                                 & (ref_df['feature'] == 'transcript')
                                 & (ref_df['transcript_biotype'] == 'protein_coding')]
            transcripts = sorted(list(tsl1_df.transcript_id))[::-1]

        axes[i].set_xlim(start - offset/2, end + offset/2)

        yloc = 1
        exons = []
        high = len(transcripts) + 1
        for t in transcripts:
            tmp = ref_df.loc[ref_df.transcript_id == t]
            tmp.reset_index(inplace=True, drop=True)

            t_id = tmp.at[0, 'transcript_id']
            biot = tmp.at[0, 'transcript_biotype']
            tsl = tmp.at[0, 'transcript_support_level']
            t_start = tmp.at[0, 'start'] - start_
            t_end = tmp.at[0, 'end'] - start_

            exon_list = [(tmp.at[i, 'start'] - start_,
                          tmp.at[i, 'end'] - tmp.at[i, 'start']) for i in tmp.index[1:]]

            axes[i].hlines(yloc, t_start, t_end, color='#fb9a99')
            axes[i].text(0, yloc + 0.4 * STEPS,
                         '{} ({}) tsl{}'.format(t_id, biot, int(tsl)),
                         fontsize=FONTSIZE)

            for e in exon_list:
                s = e[0]
                w = e[1]

                rect = Rectangle((s, yloc-BOX_HEIGHT/2), w, BOX_HEIGHT)
                exons.append(rect)

            # DRAW THE SNORNAs OF THE HOST !!
            draw_sno_in_host(axes[i], start_, t_start, t_end,
                             sno_in_host, yloc, BOX_HEIGHT, high)

            # Create patch collection with specified colour/alpha
            pc = PatchCollection(exons, facecolor='#fb9a99', alpha=.95,
                                 edgecolor=None, zorder=200)

            # Add collection to axes
            axes[i].add_collection(pc)

            yloc += STEPS
            exons.clear()

        axes[i].set_ylim(0, yloc)
        axes[i].get_yaxis().set_visible(False)
        axes[i].set_title('{} - {}'.format(gene_name, titles[i]))

    axes[0].get_xaxis().set_visible(False)
    plt.tight_layout(pad=0.5)
    # plt.savefig('/data/labmeetings/6mars_michelle/VIPR1_transcripts.svg',
    #             format='svg', transparent=True)
    plt.show()
