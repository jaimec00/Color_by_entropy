# ----------------------------------------------------------------------------------------------------------------------
# Colorme.py, by Jaime Cardenas
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
from Bio import pairwise2
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def colorme(obj, thresh):
    # ------------------------------------------------------------------------------------------------------------------
    # Get consensus sequence
    consensus_file = pd.read_csv('con_seq_entropy.csv')
    consensus_seq_old = consensus_file.iloc[:, 0]
    consensus_seq = ''
    for residue in consensus_seq_old:
        consensus_seq += residue
    # ------------------------------------------------------------------------------------------------------------------
    # Get all chains
    fasta = cmd.get_fastastr('model ' + obj)
    fasta = fasta.split()
    all_chains = []
    for line in fasta:
        if '>' in line:
            j = -1
            temp_chains = ''
            while line[j] != '_':
                temp_chains += line[j]
                j -= 1
            all_chains += [temp_chains]
    # ------------------------------------------------------------------------------------------------------------------
    # Get info per chain
    for chainletter in all_chains:
        object_sequence_dict = {'object_sequence': []}
        cmd.iterate_state(-1, "model %s and name ca and chain %s" % (obj, chainletter), \
                          "object_sequence.append([oneletter, resi, chain])", space=object_sequence_dict)
        # --------------------------------------------------------------------------------------------------------------
        # Get object sequence
        obj_seq = ''
        for i in object_sequence_dict['object_sequence']:
            obj_seq += i[0]
            
        # --------------------------------------------------------------------------------------------------------------
        # Parameters for alignment
        gap_open_consensus = -10000
        gap_open_obj = -6
        gap_extend_obj = -0.5
        # --------------------------------------------------------------------------------------------------------------
        # Align locally to test sequence
        local_alignments = pairwise2.align.localxd(sequenceA=consensus_seq, sequenceB=obj_seq,
                                                   openA=gap_open_consensus,
                                                   extendA=gap_extend_obj,
                                                   openB=gap_open_obj,
                                                   extendB=gap_extend_obj)
        # --------------------------------------------------------------------------------------------------------------
        if local_alignments[0][2] <= 40:
            cmd.color('red', 'model %s and chain %s' % (obj, chainletter))
        else:
            messy_alignments = pairwise2.align.globalxd(sequenceA=consensus_seq, sequenceB=obj_seq,
                                                        openA=gap_open_consensus,
                                                        extendA=gap_extend_obj,
                                                        openB=gap_open_obj,
                                                        extendB=gap_extend_obj)

            alignment = pairwise2.format_alignment(*messy_alignments[0])
            # ----------------------------------------------------------------------------------------------------------
            # Get aligned object sequence
            obj_seq_aligned_long = alignment[(len(consensus_seq) * 2 + 2):]
            obj_seq_aligned = ''
            character = 0
            while obj_seq_aligned_long[character:character + 2] != 'Sc':
                obj_seq_aligned += obj_seq_aligned_long[character]
                character += 1
            # ----------------------------------------------------------------------------------------------------------
            # Find the positions of the aligned object sequence residues
            idx = []
            k = 0
            i = 0
            while k < len(obj_seq) and i < len(obj_seq_aligned):
                if obj_seq_aligned[i] == obj_seq[k]:
                    idx += [i]
                    k += 1
                i += 1
            # ----------------------------------------------------------------------------------------------------------
            # Get the entropy values of each residue in sample sequence by using real_sample_idx
            entropy = []
            for num in idx:
                entropy += [consensus_file.iloc[num, 1]]
            # ----------------------------------------------------------------------------------------------------------
            # Get the chain, one-letter code, and position for each residue
            chain, name, position = [], [], []
            for i in object_sequence_dict['object_sequence']:
                name += i[0]
                position += i[1:-1]
                chain += i[-1]
            # ----------------------------------------------------------------------------------------------------------
            # Color residues based on entropy values
            for val in range(len(entropy)):
                thresh = float(thresh)
                if entropy[val] > thresh:
                    col = 'red'
                else:
                    perc = entropy[val] / thresh

                    startc = 550
                    endc = 615
                    formperc = str(int((endc - startc) * (1 - perc) + startc))

                    col = 'c' + formperc
                # ------------------------------------------------------------------------------------------------------
                cmd.color(col, 'model %s and chain %s and resi %s' % (obj, chain[val], position[val]))


# ----------------------------------------------------------------------------------------------------------------------
cmd.extend('colorme', colorme)
