import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

gene_tpm = 'tpm.tsv'
transcript_tpm = 'transcript_tpm.tsv'
ref = 'gene_name.txt'


def load_df(file, corresp_dict, GOI):
    df = pd.read_csv(file, sep='\t',
                     dtype={'ctrl': float,
                               'XRN1': float,
                            'SMG6-XRN1': float,
                            'UPF1-XRN1': float})
    # for c in df.columns[1:]:
    #     df[c] = df[c].map(int)

    df['gene'] = df['transcript'].map(corresp_dict)
    df = df.loc[df.gene == GOI]
    df.set_index('transcript', inplace=True)
    df.drop('gene', axis=1, inplace=True)

    df.sort_index(inplace=True)
    print(df)

    return df


def trans_of_genes(file):
    corresp_dict = {}
    with open(file, 'r') as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            corresp_dict[fields[0]] = fields[1]
    return corresp_dict


def load_process_colombo(corresp_dict, GOI):

    THRES_MIN = 5
    MIN_RATIO = 1

    file = '/home/danx/Documents/projects/NMD_RNA_seq/results/colombo/kallisto/NAME_CORRECTED_transcript_tpm.tsv'
    df = pd.read_csv(file, sep='\t')
    df.set_index('transcript', inplace=True)
    df = df.transpose()

    new_index = []
    for idx in df.index:
        if 'ctrl' in idx:
            new_index.append('ctrl')
        elif 'KD' in idx:
            new_index.append(idx[:-1])
        elif 'double' in idx:
            new_index.append(idx[:-2])
        else:
            new_index.append(idx[:-1])

    df.index = new_index
    df.reset_index(inplace=True)
    df = df.groupby('index').mean().transpose()

    # full
    # df = df[['ctrl', 'UPF1_KD', 'SMG6_KD', 'SMG7_KD', 'double_KD', 'UPF1_res',
    #          'SMG6_res', 'SMG7_res', 'double_resSMG6', 'double_resSMG7']]

    # partial
    df = df[['ctrl', 'UPF1_KD', 'SMG6_KD', 'SMG7_KD', 'UPF1_res', 'SMG6_res', 'SMG7_res']]

    df['gene'] = df.index.map(corresp_dict)
    df = df.loc[df.gene == GOI]
    df.drop(columns=['gene'], inplace=True)
    df.sort_index(inplace=True)
    print(df)

    df['avg'] = df.mean(axis=1)
    df = df.loc[df['avg'] > THRES_MIN]
    df.drop(columns=['avg'], inplace=True)

    df['max'] = df[['UPF1_KD', 'SMG6_KD', 'SMG7_KD']].max(axis=1)
    df['ratio'] = df['max'] / df['ctrl']
    df = df.loc[df['ratio'] > MIN_RATIO]
    df.drop(columns=['max', 'ratio'], inplace=True)

    print(df)

    return df


def graph_trans(transcript_df, yaxix_label):

    labels = transcript_df.index
    ctrl = transcript_df.ctrl
    XRN1 = transcript_df.XRN1
    SMG6_XRN1 = transcript_df['SMG6-XRN1']
    UPF1_XRN1 = transcript_df['UPF1-XRN1']

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    axes = []
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], wspace=0.025, hspace=0.6, bottom=0.15)
    fig = plt.figure()
    axes.append(plt.subplot(gs[0]))
    axes.append(plt.subplot(gs[1]))

    rects1 = axes[0].bar(x - 2*width/2 + width/4, ctrl, width/2, label='ctrl')
    rects2 = axes[0].bar(x - width/2 + width/4, XRN1, width/2, label='XRN1')
    rects3 = axes[0].bar(x + width/4, SMG6_XRN1, width/2, label='SMG6-XRN1')
    rects4 = axes[0].bar(x + width/2 + width/4, UPF1_XRN1, width/2, label='UPF1-XRN1')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    axes[0].set_ylabel(yaxix_label, fontsize=12)
    axes[0].set_title('Transcripts relative changes', fontsize=15)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    if len(transcript_df) != 0:
        axes[0].legend()

    plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # def autolabel(rects):
    #     """Attach a text label above each bar in *rects*, displaying its height."""
    #     for rect in rects:
    #         height = rect.get_height()
    #         axes[0].annotate('{}'.format(height),
    #                     xy=(rect.get_x() + rect.get_width() / 2, height),
    #                     xytext=(0, 3),  # 3 points vertical offset
    #                     textcoords="offset points",
    #                     ha='center', va='bottom')
    # autolabel(rects1)
    # autolabel(rects2)
    # autolabel(rects3)
    # autolabel(rects4)

    # fig.subplots_adjust(bottom=0.5)

    return axes, fig


def process(transcript_df):

    # graph_trans(transcript_df, "TPM")

    filter_df = transcript_df.copy(deep=True)
    filter_df['mean'] = filter_df.mean(axis=1)
    filter_df = filter_df.loc[filter_df['mean'] > 2]

    filter_df['max'] = filter_df[['SMG6-XRN1', 'UPF1-XRN1']].values.max(axis=1)
    filter_df['min'] = filter_df[['XRN1', 'ctrl']].values.min(axis=1)
    filter_df['ratio'] = filter_df['max'] / filter_df['min']

    filter_df = filter_df.loc[(filter_df['ratio'] > 2)]

    filter_df.drop(columns=['mean', 'max', 'min', 'ratio'], inplace=True)

    # print(filter_df)
    #
    # relative_df = filter_df.copy(deep=True)
    # for i in relative_df.index:
    #     relative_df.at[i, 'XRN1'] = relative_df.at[i, 'XRN1'] / relative_df.at[i, 'ctrl']
    #     relative_df.at[i, 'SMG6-XRN1'] = relative_df.at[i, 'SMG6-XRN1'] / relative_df.at[i, 'ctrl']
    #     relative_df.at[i, 'UPF1-XRN1'] = relative_df.at[i, 'UPF1-XRN1'] / relative_df.at[i, 'ctrl']
    #     relative_df.at[i, 'ctrl'] = relative_df.at[i, 'ctrl'] / relative_df.at[i, 'ctrl']
    # print(relative_df)
    #
    # graph_trans(relative_df, "Relative TPM\nnormalized to ctrl")

    return filter_df


def graph_col(df, ax, fig):

    x = np.arange(len(df.index))  # the label locations
    width = 0.1  # the width of the bars
    start = (len(df.columns)/2)*width

    for i, col in enumerate(df.columns):
        ax.bar(x - start + width*i + 0.5*width, list(df[col]), width, label=col)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('TPM', fontsize=12)
    ax.set_title('Transcripts relative changes', fontsize=15)
    ax.set_xticks(x)
    ax.set_xticklabels(df.index)
    if len(df) != 0:
        ax.legend()

    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor")


def graph(df_lykke, df_colombo):

    axes, fig = graph_trans(df_lykke, "TPM")

    graph_col(df_colombo, axes[1], fig)

    # plt.tight_layout(pad=5)
    plt.show()


def main():

    corresp_dict = trans_of_genes(ref)

    GOI = 'RPL17'

    transcript_df = load_df(transcript_tpm, corresp_dict, GOI)

    proc_trans_df = process(transcript_df)
    print(proc_trans_df)
    print('------------------------------------------------------------------')

    colombo_trans_df = load_process_colombo(corresp_dict, GOI)

    graph(proc_trans_df, colombo_trans_df)


if __name__ == '__main__':
    main()
