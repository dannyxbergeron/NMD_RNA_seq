import os
import pandas as pd
import numpy as np
import matplotlib as plt

from viz_transcripts import load_ref, prepare_fig

files = [x for x in os.listdir() if x.endswith('.csv')]
ref = 'gene_name.txt'


def trans_of_genes(file):
    corresp_dict = {}
    with open(file, 'r') as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            corresp_dict[fields[0]] = fields[1]
    return corresp_dict


def get_tpms():
    file = '/home/danx/Documents/projects/NMD_RNA_seq/results/colombo/' + \
        'kallisto/NAME_CORRECTED_transcript_tpm.tsv'
    df = pd.read_csv(file, sep='\t')
    df.set_index('transcript', inplace=True)

    new_index = []
    for col in df.columns:
        if 'ctrl' in col:
            new_index.append('ctrl')
        else:
            new_index.append(col[:-1].replace('_', '').replace('KD', '_KD').replace('res', '_res'))

    df.columns = new_index
    df = df.transpose()
    df = df.reset_index()
    df = df.groupby("index").mean().transpose()

    for c in df.columns:
        df[c] = df[c].map(int)

    return df


def load_files(corresp_dict, tpms, gene):

    # filters !!!
    p_value = 0.05
    fc_treshold = 1

    colnames = ['transcript', 'baseMean', 'log2FoldChange',
                'lfcSE', 'stat', 'pvalue', 'padj']

    dfs = {}
    for f in files:
        name = f.replace('_.csv', '').replace('.csv', '').replace('ctrl2_', 'ctrl')
        cond1, cond2 = name.split('-')[0], name.split('-')[1]

        df = pd.read_csv(f)
        df.columns = colnames
        df['gene'] = df.transcript.map(corresp_dict)
        df.dropna(axis=0, subset=['padj'], inplace=True)

        df = df.loc[df.gene == gene]
        df = df.loc[df.padj < p_value]
        df = df.loc[(df.log2FoldChange <= -fc_treshold) | (df.log2FoldChange >= fc_treshold)]
        df = df.drop(columns=['lfcSE', 'stat', 'pvalue', 'gene', 'baseMean'])
        df.sort_values('padj', axis=0, inplace=True)

        df[cond1] = df.transcript.map(tpms[cond1].to_dict())
        df[cond2] = df.transcript.map(tpms[cond2].to_dict())
        # print('{}, number of rows: {}'.format(name, len(df)))
        # print(df)
        # print('\n----------------------------------------------------------------------------------\n')

        dfs[name] = df

    return dfs


def analysis(dfs, gene, ref_df_):

    # Loading the annotation data ! -------
    ref_df, sno_in_host = load_ref(gene, ref_df_)
    # -------------------------------------

    # print('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')

    affected_transcripts = []
    ctrl_list = [x for x in dfs.keys() if x.startswith('ctrl')]
    for name in ctrl_list:

        df = dfs[name].loc[dfs[name].log2FoldChange < 0]
        transcript_set = set(df.transcript)

        second_name = name.split('-')[1]
        res = [x for x in dfs.keys() if second_name in x and x != name]
        for other in res:

            other_df = dfs[other].loc[dfs[other].log2FoldChange > 0]
            transcript_list = transcript_set.intersection(set(other_df.transcript))
            affected_transcripts += transcript_list

            if len(transcript_list) != 0:

                for factor in ['double', 'SMG6', 'UPF1', 'SMG7']:
                    if factor in other:
                        print('Condition:', factor)
                        if factor != 'double':
                            break

                # print('KD -------------------------------------------')
                new_df = df.copy(deep=True)
                new_df = new_df.loc[new_df.transcript.isin(transcript_list)]
                new_df.sort_values(by='transcript', inplace=True)
                # print(new_df)
                # print('Rescue =======================================')
                new_other_df = other_df.copy(deep=True)
                new_other_df = new_other_df.loc[new_other_df.transcript.isin(transcript_list)]
                new_other_df.sort_values(by='transcript', inplace=True)
                # print(new_other_df)
                # print('End+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

                tpm_name = new_other_df.columns[-1]
                new_df[tpm_name] = new_df.transcript.map(
                    dict(zip(new_other_df.transcript, new_other_df[tpm_name])))
                new_df['log2FC_res'] = new_df.transcript.map(
                    dict(zip(new_other_df.transcript, new_other_df['log2FoldChange'])))
                new_df['padj_res'] = new_df.transcript.map(
                    dict(zip(new_other_df.transcript, new_other_df['padj'])))
                print(new_df)
                print('----------------------------------------------------------')

    t = set(affected_transcripts)
    if len(t) != 0:
        prepare_fig(ref_df, sno_in_host, t)
    print('done')


def main():

    file = '/home/danx/Documents/projects/snoRNA_network_V2/ref/modified/' + \
        'human_ensembl_87_wo_dup_v3.BB_v3_PARSED.tsv'
    ref_df = pd.read_csv(file, sep='\t')

    # neg_correlated
    # genes = [
    #     'ATP5B',
    #     'MIR664A',
    #     'HSPA9',
    #     'RABGGTB',
    #     'RP11-1094M14.11',
    #     'RC3H2',
    #     'NCL',
    #     'PRRC2B',
    #     'PRRC2B',
    #     'PUM1',
    #     'SNHG6',
    #     'CHD8',
    #     'EIF5',
    #     'PPAN',
    #     'EEF1B2',
    #     'EIF4A2',
    #     'TSR1',
    #     'EIF4A2',
    #     'SF3B3',
    #     'RABGGTB',
    #     'RABGGTB',
    #     'TSR1',
    #     'TEX14',
    #     'SF3B3',
    #     'CHD8',
    #     'NCL',
    #     'CCT6P1',
    #     'RPL4',
    #     'NCL',
    #     'HSPA8',
    #     'WDR43',
    #     'CCAR1',
    #     'OSCP1',
    #     'PUM1',
    #     'SNHG5',
    # ]

    # pos RPL genes
    genes = [
        'RPS2',
        'RPL10',
        'RPL27A',
        'RPL5',
        'RPL23A',
        'RPL17',
        'RPL32',
        'RPL30',
        'RPL13A',
        'RPL17',
        'RPL17',
        'RPL3',
        'RPL5',
        'RPS3',
        'RPLP2',
        'RPSA',
        'RPL23A',
        'RPL18A',
        'RPL12',
        'RPL32P3',
        'RPS8',
        'RPL21',
        'RPS20',
        'RPS2',
        'RPS12',
        'RPL7A',
        'RPL27A',
        'RPL3',
        'RPL4',
        'RPL23',
        'RPS8',
        'RPS3A'
    ]

    ref_df = ref_df.loc[ref_df.gene_name.isin(genes)]

    genes = sorted(list(set(genes)))
    print(len(genes))

    for gene in genes:

        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("----------------------{}--------------------".format(gene))
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")

        corresp_dict = trans_of_genes(ref)

        tpms = get_tpms()

        dfs = load_files(corresp_dict, tpms, gene)

        analysis(dfs, gene, ref_df)


if __name__ == '__main__':
    main()
