#!/usr/bin/python
#Author: Jasen M. Jackson, Loyola '19
import os.path
import itertools
import sys
from collections import defaultdict
from io_utils import type_fa, type_fasta, fasta_check

class RedundancyMap():
    """Redundancy maps organize redundant integration site data
       and store sequences representing non-redundant integrations.
       It assumes that the FASTA files it is built on contain
       reads that the flanking site is on the 3' end of the read,
       and the repetitive element has already been trimmed away."""

    def __init__(self, fasta_path, end_size):
        self.fasta = fasta_path
        self.end_size = end_size

        ## check for fasta format
        err = "Error (redundancy_map): "
        if type_fa(fasta_path):
            self.out = fasta_path[:-3] + ".duplicates_removed.fa"
        elif type_fasta(fasta_path):
            self.out = fasta_path[:-6] + ".duplicates_removed.fasta"
        else: print(err+"fasta format required")

        ## validate end_size parameter and create redundancy map
        if not isinstance(end_size, int): print(err+"invalid end_size type")
        elif not end_size > 0: print(err+"end_size must be greater than 0")
        elif fasta_check(fasta_path):
            print("\tCreating redundancy map from "+self.fasta)
            self.rm = self.rm_make(fasta_path, end_size)

        self.unique_count = len(self.rm)
        self.duplicate_count = self.total_count - self.unique_count
        ## convert redundancy map to fasta file
        print("\t\tSaving redundancy map to "+self.out)
        if self.rm is not None:
            outfile = open(self.out, 'w')
            features = self.get_features(self.rm, 'header', 'sequence')
            for i in range(len(self.rm)):
                if features[1][i] and len(features[1][i]) >= 20:
                    outfile.write(features[0][i]+'\n')
                    outfile.write(features[1][i]+'\n')
        else: print(err+"could not create redundancy map")

    def rm_add(self, seq, seq_len, header, map, end_size):
        ## adds non-redundant sequence info to redundancy map data structure
        new_key = seq[-1*end_size:]
        new_value = {'header': header, 'sequence': seq, 'length': seq_len, 'depth': 1}
        if not map[new_key]: map[new_key] = new_value
        elif seq_len > map[new_key]['length']:
            new_value['depth'] = map[new_key]['depth'] + 1
            map[new_key] = new_value
        else: map[new_key]['depth'] += 1

    def rm_make(self, input, end_size):
        ## make redundancy map
        nested_dict = lambda: defaultdict(nested_dict)
        rm = nested_dict() # rm = redundancy map
        fi = open(input, 'r')
        seq_len, self.total_count = int(), int()
        while True:
            header = fi.readline()[:-1]
            seq = fi.readline()[:-1]
            if not seq: break
            self.total_count += 1
            seq_len = len(seq)
            if (seq_len > end_size): self.rm_add(seq,seq_len,header,rm,end_size)
        return rm

    def print_map(self):
        print("\nRedundancy map ("+str(self.unique_count)+" of "+str(self.unique_count)+" shown) \n")
        print_count = int()
        for p_key, p_feats in self.rm.items():
            print(p_key)
            for feat in p_feats: print("\t"+feat+':'+str(p_feats[feat]))
            print(""); print_count += 1

    def print_map_head(self, n):
        print("\nRedundancy map ("+str(n)+" of "+str(self.unique_count)+" shown) \n")
        print_count = int()
        for p_key, p_feats in self.rm.items():
            if print_count < n:
                print(p_key)
                for feat in p_feats: print("\t"+feat+':'+str(p_feats[feat]))
                print(""); print_count += 1
            else: break

    def get_features(self, map, *features):
        ## return list of features lists // TODO: error handling
        feat_lists = list()
        for feat in features:
            feat_list = [value[feat] for key, value in map.items()]
            feat_lists.append(feat_list)
        return feat_lists


'''
TODO:
0. modify rm_feats to return list of features
1. Convert rm into class
2. Count duplicates & uniques, store as rm attribute
3. Use dups on trimmed_file
4. Align trimmed_dups_removed file to reference
5. Clean up i/o with argparse
'''
