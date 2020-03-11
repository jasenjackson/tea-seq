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

    def __init__(self, r1_path, r2_path, ltr, end_size, verbose=True):
        self.fastq_r1 = r1_path
        self.fastq_r2 = r2_path
        self.ltr = ltr
        self.end_size = end_size
        self.verbose = verbose
        self.total_count = int()
        self.sorted_keys = defaultdict(list)# key = feature, value = list
        self.out1 = self.fastq_r1[:-6] + ".junctions.FASTQ"
        self.out2 = self.fastq_r2[:-6] + ".junctions.FASTQ"
        if not self.end_size > 0:
            print(err+"end_size must be greater than 0")
            return
        else: self.build_from_paired()

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
                print("Processed " + str(self.total_count) + " reads. " + str(len(self.map)) + " junctions identified. ")

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
        libmap_copy = {}
        for junction in self.map:
            if self.map[junction]['ltr_read'] == 'r1':
                if (kill_sequence not in self.map[junction]['r1_sequence']):
                    libmap_copy[junction] = self.map[junction]
                else: pass
            if self.map[junction]['ltr_read'] == 'r2':
                if (kill_sequence not in self.map[junction]['r2_sequence']):
                    libmap_copy[junction] = self.map[junction]
                else: pass
        self.map = libmap_copy

    def get_sorted_keys(self, feature):
        #if self.sorted_keys[feature]: return self.sorted_keys[feature]
        key_list, feat_list = list(), list()
        for key in self.map:
            key_list.append(key)
            feat_list.append(self.map[key][feature])
        key_list = [x for _,x in sorted(zip(feat_list,key_list))]
        key_list.reverse()
        self.sorted_keys[feature] = key_list
        return self.sorted_keys[feature]

    def saveMap(self,*out):
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

    def writeFASTQ(self):

        fastq_output1 = "" # output for R1 file
        fastq_output2 = "" # output for R2 file
        for junction in self.map:
            # generate output for R1
            fastq_output1 = fastq_output1 + self.map[junction]['r1_header'] + "\n"
            fastq_output1 = fastq_output1 + self.map[junction]['r1_sequence'] + "\n"
            fastq_output1 = fastq_output1 + "+" + "\n"
            fastq_output1 = fastq_output1 + self.map[junction]['r1_quality'] + "\n"
            # generate output for R2
            fastq_output2 = fastq_output2 + self.map[junction]['r2_header'] + "\n"
            fastq_output2 = fastq_output2 + self.map[junction]['r2_sequence'] + "\n"
            fastq_output2 = fastq_output2 + "+" + "\n"
            fastq_output2 = fastq_output2 + self.map[junction]['r2_quality'] + "\n"

        # write to outfile1
        outfile1 = open(self.out1, 'w')
        outfile1.write(fastq_output1)
        outfile1.close()

        # write to outfile2
        outfile2 = open(self.out2, 'w')
        outfile2.write(fastq_output2)
        outfile2.close()

    def writeFlanks(self): # not tested yet
        fasta_output = "" # output for R1 file
        for junction in self.map:
            fasta_output = fasta_output + ">" + self.map[junction]['r1_header'] + "\n"
            fasta_output = fasta_output + self.map[junction]['r1_s'] + "\n"
        out_path = self.out1[:-6] + ".fasta"
        outfile = open(out_path, 'w')
        outfile.write(fasta_output)
        outfile1.close()

    def map_view(self, max_lines=20, feature='depth', view="head", max_char=100):
        '''print (top/bottom) junctions ranked by feature'''

        ## personalize print settings based on parameters
        lines_displayed = int()
        if max_lines > len(self.map): lines_displayed = len(self.map)
        else: lines_displayed = max_lines
        keys_list = self.get_sorted_keys(feature)
        if view=="head":
            display = "top"
            keys_list = keys_list[:max_lines]
        if view=="tail":
            display = "bottom"
            keys_list = keys_list[-1*max_lines:]

        ## display print parameters
        print("~~~ map "+view+" ~~~\n"+
              "Junction map ("+ display+" "+str(lines_displayed)+" of "+
              str(len(self.map))+" shown)\n"+
              "Sorted by "+feature+" in descending order.")

        ## print junction map
        for key in keys_list:
            print("\n"+key)
            for feat in self.map[key]:
                if isinstance(self.map[key][feat],str):
                    feat_value=self.map[key][feat][:max_char]
                else: feat_value=self.map[key][feat]
                print("\t"+feat+": "+str(feat_value))
        print("~~~~~~~~~~~~")

    def get_features(self, map, *features):
        ## return list of features lists // TODO: error handling
        feat_lists = list()
        for feat in features:
            feat_list = [value[feat] for key, value in map.items()]
            feat_lists.append(feat_list)
        return feat_lists

def main():

    ## build hash map of genome-LTR junctions
    lib = junction_map("data/Glycine_GMR30/HL-2_S2_L001_R1_001.fastq","data/Glycine_GMR30/HL-2_S2_L001_R2_001.fastq","TGTTAGCCCATA",30,verbose=True)

    ## look at head of library
    lib.map_view(max_lines=20,feature='depth', view="head",max_char=100); print("")
    lib.map_view(max_lines=20,feature='depth', view="tail",max_char=100); print("")

    ## remove poly-purine tract (TODO: identify better kill sequence)
    lib.remove_junction("AATCGAAGCAGACATTTTTTGGGGGCAATT")

    ## write FASTQ (one entry per junction)
    lib.writeFASTQ()

    ## write FASTA file containing sequence of flank (LTR removed)
    lib.writeFlanks()

    ###### PROGRESS ######
    ## TODO remove adapters ***
    #lib.remove_adapter("Adapter")

    ## TODO quality trim (skip for now)
    #lib.quality_trim()

    ## align junctions with bowtie ***
    ## add nearby genes with gtf file
    ## compare to known integrations
    ## export results & feed into Rplot ***

    ## align flanks with bowtie
    ## classify each junction as 'composite','unknown_repeat'or'non-repetitive'
    ## export results & feed into rplot
    ## add nearby genes with gtf file
    ## compare to known integrations

    ## build functions to study each gene type / repeat type
    ## build functions to compare junction maps

if __name__ == "__main__":
    main()
