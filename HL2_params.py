WORDSIZE = 8 #kmers
#MAXDIST = 6 # max hamming distance for string comparison
#BREAKDIST = 3 # minimum hamming distance to break slidingSequence() loop

DATA_DIR = 'data/Glycine_GMR30'
TEMP_DIR = 'temp'
BIN_DIR = 'bin'

R1_FILES = [
            'HL-2_S2_L001_R1_001.fastq'
           ]

R2_FILES = [
            'HL-2_S2_L001_R2_001.fastq'
           ]
LIBS = [
        'HL-2'
       ]

#Feature name,sequence, coverage threshold, required/optional parameter
FEATURES = [
    ["SPLINK1","CGTGGCTGAATGAGACTGGTGTCGACACTAGTGGT",0.6, 'adapter'],
    #["SPLINK1_RC","ACCACTAGTGTCGACACCAGTCTCATTCAGCCACG",0.6, 'adapter'],    
    ["GMR30_LTR_END", "TGTTAGCCCATA", 1.0, 'element'],
    ["GMR30_LTR_END_RC", "TATGGGCTAACA", 1.0, 'element'],
    ["Polypurine Track", "GGGCAATT", 1.0, 'remove'],
    ["Polypurine Track RC", "AATTGCCC", 1.0, 'remove'],
    ["Primer 1404 FW","ACATGCTCCTGAGATTCACTAGT",0.9,'primer'],
    ["Primer 1404 RC", "ACTAGTGAATCTCAGGAGCATGT", 0.9, 'primer'],
    ["Primer 1405 FW", "GCCCTCTTTCTTGTCTTAGCCTT", 0.9, 'primer'],
    ["Primer 1405 RC", "AAGGCTAAGACAAGAAAGAGGGC", 0.9, 'primer'],
    ["Primer 1406 FW", "GTCCTCTTTCTTGTCTTAGCTTT", 0.9, 'primer'],
    ["Primer 1406 RC", "AAAGCTAAGACAAGAAAGAGGAC", 0.9, 'primer'],
    ["Primer 1407 FW", "GCCCTCTTTCTTGTCCTAGCTT", 0.9, 'primer'],
    ["Primer 1407 RC", "AAGCTAGGACAAGAAAGAGGGC", 0.9, 'primer']
]
