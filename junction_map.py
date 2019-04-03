#!/usr/bin/python
#Author: Jasen M. Jackson, Loyola '19
import os.path
import itertools
import sys
from collections import defaultdict
from feature_search import *
from io import *

class junction_map():
    """Redundancy maps organize redundant integration site data
       and store sequences representing non-redundant integrations.
       It assumes that the FASTA files it is built on contain
       reads that the flanking site is on the 3' end of the read,
       and the repetitive element has already been trimmed away."""

    def __init__(self, r1_path, r2_path, ltr, end_size, verbose):
        self.fastq_r1 = r1_path
        self.fastq_r2 = r2_path
        self.ltr = ltr
        self.end_size = end_size
        self.verbose = verbose
        self.unique_count = int()
        self.duplicate_count = int()
        self.total_count = int()
        self.out1 = self.fastq_r1[:-6] + ".junctions.FASTQ"
        self.out2 = self.fastq_r2[:-6] + ".junctions.FASTQ"
        if not self.end_size > 0:
            print(err+"end_size must be greater than 0")
            return
        else: self.build_from_paired()

        ## compute counts
        self.unique_count = len(self.map)
        self.duplicate_count = self.total_count - self.unique_count


    def build_from_merged(self, file):
        if self.verbose == True: print("Building genome junction map from reads")
        nested_dict = lambda: defaultdict(nested_dict)
        self.map = nested_dict()
        reads = open(file, 'r')
        while True:
            header = reads.readline()[:-1]
            seq = reads.readline()[:-1]
            if not seq: break
            self.total_count += 1

            ## locate genome-LTR junction & add to map
            if (self.verbose == True):
                print("Identified "+ str(self.total_count) + " junctions")
            ltr_pos, ltr_dist = kmer_search(seq, self.ltr, 1, 0)
            self.add_junction(seq,header,ltr_pos)

    def build_from_paired(self):
        if self.verbose == True: print("Building genome junction map from reads")
        nested_dict = lambda: defaultdict(nested_dict)
        self.map = nested_dict()
        r1, r2 = open(self.fastq_r1, 'r'), open(self.fastq_r2, 'r')
        while True:
            junction = {}; self.total_count += 1

            ## update user
            if (self.verbose == True and self.total_count % 1000 == 0):
                print("Identified "+ str(self.total_count) + " junctions. ")

            ## grab read data
            junction["r1_header"] = r1.readline()[:-1]
            junction["r2_header"] = r2.readline()[:-1]
            junction["r1_sequence"] = r1.readline()[:-1]
            junction["r2_sequence"] = r2.readline()[:-1]
            comment = r1.readline(); comment = r2.readline()
            junction["r1_quality"] = r1.readline()[:-1]
            junction["r2_quality"] = r2.readline()[:-1]
            if not comment: break

            ## identify genome-LTR junction and add to junction map
            r1_pos, r1_dist = kmer_search(junction['r1_sequence'], self.ltr, 1, 0)
            r2_pos, r2_dist = kmer_search(junction['r2_sequence'], self.ltr, 1, 0)
            if r1_dist <= r2_dist:
                junction = self.add_junction(self.paired_junction(junction, r1_pos, "r1"))
            if r2_dist < r1_dist:
                junction = self.add_junction(self.paired_junction(junction, r2_pos, "r2"))

    def paired_junction(self, junction, pos, read):
        ## return an appropriate dictionary for a paired junction
        ltr_read = read + "_sequence"
        junction['flank'] = junction[ltr_read][:pos]
        junction['flank_length'] = len(junction['flank'])
        junction['ltr_position'] = pos
        junction['ltr'] = junction[ltr_read][pos:]
        junction['ltr_read'] = read
        junction['depth'] = 1
        return junction

    def add_junction(self, junction):
        key = junction['flank'][-1*self.end_size:]
        value = junction
        if junction['flank_length'] < self.end_size:
            return
        if not self.map[key]:
            self.map[key] = junction
        elif junction['flank_length'] > self.map[key]['flank_length']:
            junction['depth'] = self.map[key]['depth'] + 1
            self.map[key] = junction
        else: self.map[key]['depth'] += 1

    def remove_junction(self, kill_sequence):
        for key in self.map:
            if (key == kill_sequence): print(key)

    def save(self,*out):
        if out: out_path = out[0]
        else: out_path = self.out
        print("\tSaving to "+out_path)
        if self.map is not None:
            outfile = open(out_path, 'w')
            features = self.get_features(self.map, 'header', 'sequence')
            for i in range(len(self.map)):
                outfile.write(features[0][i]+'\n')
                outfile.write(features[1][i]+'\n')
        else: print(err+"could not save redundancy map")

    def print_map(self):
        print("\nRedundancy map ("+str(self.unique_count)+" of "+str(self.unique_count)+" shown) \n")
        print_count = int()
        for p_key, p_feats in self.map.items():
            print(p_key+" * "+self.map[p_key]['ltr'])
            for feat in p_feats: print("\t"+feat+':'+str(p_feats[feat]))
            print(""); print_count += 1

    def print_map_head(self, n):
        n_shown = int()
        if n > self.unique_count: n_shown = self.unique_count # can't show more junctions than we have
        else: n_shown = n
        print("\nRedundancy map ("+str(n_shown)+" of "+str(self.unique_count)+" shown) \n")
        print_count = int()
        for p_key, p_feats in self.map.items():
            if print_count < n:
                print(p_key)
                for feat in p_feats: print("\t"+feat+': '+str(p_feats[feat]))
                print(""); print_count += 1
            else: break

    def get_features(self, map, *features):
        ## return list of features lists // TODO: error handling
        feat_lists = list()
        for feat in features:
            feat_list = [value[feat] for key, value in map.items()]
            feat_lists.append(feat_list)
        return feat_lists


def main():

    ## build hash map of genome junctions (takes a while)
    lib = junction_map("r1.fastq","r2.fastq","GMR30",20, True)
    #lib = junction_map("data/Glycine_GMR30/HL-2_S2_L001_R1_001.fastq","data/Glycine_GMR30/HL-2_S2_L001_R2_001.fastq","TGTTAGCCCATA",20)

    lib.print_map()
    print("Uniques: "+ str(lib.unique_count))
    print("Duplicates: "+ str(lib.duplicate_count))
    print("Total reads: "+str(lib.total_count))

    ## TODO remove internal sequences (error sensitive)
    lib.remove_junction("Polypurine-Tract")

    ## TODO remove adapters
    lib.remove_adapter("Adapter")

    ## save map as FASTQ w/ LTR
    lib.save(mode="junction","HL2.GMR30.junctions.r1.FASTQ", "HL2.GMR30.junctions.r2.FASTQ")

    ## save map as FASTQ w/o LTR
    lib.save(mode="flanking","HL2.GMR30.flanking.r1.FASTQ","HL2.GMR30.flanking.r2.FASTQ")

if __name__ == "__main__":
    main()
