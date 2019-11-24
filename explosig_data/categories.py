from itertools import product
import math
from Bio.Seq import reverse_complement

from .constants import *


'''
Single-base substitution categories
'''
# Names
def SBS_6_category_name(row):
    return SBS_category_name_helper('', row[COLNAME.REF.value], row[COLNAME.VAR.value], '', 0)[1:4]

def SBS_12_category_name(row):
    if row[COLNAME.TSTRAND.value] == "+":
        return '%s>%s' % (row[COLNAME.REF.value], row[COLNAME.VAR.value])
    elif row[COLNAME.TSTRAND.value] == "-":
        return '%s>%s' % (BASE_PAIR[row[COLNAME.REF.value]], BASE_PAIR[row[COLNAME.VAR.value]])
    else:
        return None
        #raise ValueError("SBS_12 categorization is not possible without transcription strand information")

def SBS_96_category_name(row):
    return SBS_category_name_helper(row[COLNAME.FPRIME.value], row[COLNAME.REF.value], row[COLNAME.VAR.value], row[COLNAME.TPRIME.value], 1)

def SBS_192_category_name(row):
    five_prime = row[COLNAME.FPRIME.value][-1:]
    three_prime = row[COLNAME.TPRIME.value][:1]
    # if row[COLNAME.TSTRAND.value] is TSTRAND_VALS.PLUS: doesn't work - should we worry about casting?
    if row[COLNAME.TSTRAND.value] == "+":
        return '%s[%s>%s]%s' % (five_prime, row[COLNAME.REF.value], row[COLNAME.VAR.value], three_prime)
    elif row[COLNAME.TSTRAND.value] == "-":
        return '%s[%s>%s]%s' % (reverse_complement(three_prime), BASE_PAIR[row[COLNAME.REF.value]], BASE_PAIR[row[COLNAME.VAR.value]], reverse_complement(five_prime))
    else:
        return None

def SBS_1536_category_name(row):
    return SBS_category_name_helper(row[COLNAME.FPRIME.value], row[COLNAME.REF.value], row[COLNAME.VAR.value], row[COLNAME.TPRIME.value], 2)

def SBS_category_name_helper(five_prime, ref, variant, three_prime, flanking_size):
    assert len(five_prime) >= flanking_size
    assert len(three_prime) >= flanking_size
    five_prime = five_prime[-flanking_size:]
    three_prime = three_prime[:flanking_size]
    # Check that this is an SBS
    if (ref in BASES) and (variant in BASES):
        # Reverse complement if needed
        if ref in PURINES:
            ref = BASE_PAIR[ref]
            variant = reverse_complement(variant)
            five_prime_orig = five_prime
            five_prime	= reverse_complement(three_prime)
            three_prime	= reverse_complement(five_prime_orig)
        cat_name = '%s[%s>%s]%s' % (five_prime, ref, variant, three_prime)
    else:
        raise ValueError("Received a mutation that is not a single base substitution.")
    return cat_name

# Lists
def SBS_6_category_list():
    return ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

def SBS_12_category_list():
    return ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'G>A', 'G>C', 'G>T', 'A>C', 'A>G', 'A>T']

def SBS_96_category_list():
    return SBS_category_list_helper(flanking_size=1)

def SBS_192_category_list():
    return [ '%s[%s]%s' % (five, substitution, three)
    			for substitution in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'G>A', 'G>C', 'G>T', 'A>C', 'A>G', 'A>T']
                for five in [ ''.join(list(t)) for t in product(BASES, repeat=1) ]
                for three in [ ''.join(list(t)) for t in product(BASES, repeat=1) ] ]

def SBS_1536_category_list():
    return SBS_category_list_helper(flanking_size=2)

def SBS_category_list_helper(flanking_size):
    return [ SBS_category_name_helper(five, ref, variant, three, flanking_size)
                for ref, variant in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
                for five in [ ''.join(list(t)) for t in product(BASES, repeat=flanking_size) ]
                for three in [ ''.join(list(t)) for t in product(BASES, repeat=flanking_size) ] ]


'''
Doublet-base substitution categories
'''
# Names
def DBS_10_category_name(row):
    ref = row[COLNAME.REF.value]
    if len(ref) == 2: # doublet-base substitution
        if ref in ['AA', 'AG', 'CA', 'GA', 'GG', 'GT']:
            ref = reverse_complement(ref)
        # Implicit else: no reverse complement needed.
        cat_name = ('%s>NN' % (ref))
    else:
        raise ValueError("Received a mutation that is not a doublet base substitution.")
    return cat_name

def DBS_78_category_name(row):
    ref = row[COLNAME.REF.value]
    variant = row[COLNAME.VAR.value]
    if len(ref) == 2 and len(variant) == 2: # doublet-base substitution
        if ref in ['AA', 'AG', 'CA', 'GA', 'GG', 'GT']:
            ref = reverse_complement(ref)
            variant = reverse_complement(variant)
        elif (ref == 'AT' and variant in ['TG', 'GG', 'TC']) or \
            (ref == 'TA' and variant in ['AG', 'CC', 'AC']) or \
            (ref == 'CG' and variant in ['AC', 'AA', 'GA']) or \
            (ref == 'GC' and variant in ['CT', 'TT', 'TG']):
            # ref == reverse complement of ref
            variant = reverse_complement(variant)
        # Implicit else: no reverse complements needed.
        # See https://www.synapse.org/#!Synapse:syn11801895 for table of DBS mutations and their reverse complements
        cat_name = ('%s>%s' % (ref, variant))
    else:
        raise ValueError("Received a mutation that is not a doublet base substitution.")
    return cat_name

# Lists
def DBS_10_category_list():
    return DBS_category_list_helper(ref_only=True)

def DBS_78_category_list():
    return DBS_category_list_helper(ref_only=False)

def DBS_category_list_helper(ref_only=False):
    refs = set([''.join(tup) for tup in product(BASES, repeat=2)]) - set(['AA', 'AG', 'CA', 'GA', 'GG', 'GT'])
    def dbs_cats_for_ref(ref):
        vars = set([''.join(tup) for tup in product(set(BASES) - set([ref[0]]), set(BASES) - set([ref[1]]))])
        if ref[1] == BASE_PAIR[ref[0]]: # AT, CG, TA, GC (6-item var)
            if ref == 'AT': vars -= set(['TG', 'GG', 'TC'])
            if ref == 'TA': vars -= set(['AG', 'CC', 'AC'])
            if ref == 'CG': vars -= set(['AC', 'AA', 'GA'])
            if ref == 'GC': vars -= set(['CT', 'TT', 'TG'])
        return zip([ref] * len(vars), vars)
    cats = []
    if ref_only:
        # Only list the 10 AC>NN categories
        cats = [ DBS_10_category_name({REF:ref}) for ref in refs ]
    else:
        # List all 78
        for ref in refs:
            cats += [ DBS_78_category_name({REF:inner_ref, VAR:inner_var}) for inner_ref, inner_var in dbs_cats_for_ref(ref) ]
    return cats


'''
Insertion/deletion categories
'''
# Names

# certain insertions and deletions are characterized
def INDEL_Haradhvala2018_8_category_name(row):
    ref = row[COLNAME.REF.value]
    variant = row[COLNAME.VAR.value]
    # print(ref)
    # print(variant)
    if (ref == "-") or (len(variant) > len(ref)): # then insertion
        # TODO - determine the correct way to handle below issue (and the same issue with deletions)
        # most insertions have "-" as the ref and XXXX as the variant
        # but some insertions have XXXX as the ref and YYYXXXYYYYY as the variant
        # I'm currently treating the second type as one insertion of the length of both insertion sides
        # I could treat it as multiple insertions
        if ref == "-":
            insertion_size = len(variant)
        else:
            insertion_size = len(variant) - len(ref)
        assert insertion_size > 0
        if insertion_size == 1:
            return "INS1"
        elif insertion_size == 2:
            return "INS2"
        elif insertion_size == 3:
            return "INS3"
        else:
            return "INS4"
    elif (variant == "-") or (len(ref) > len(variant)):
        if variant == "-":
            deletion_size = len(ref)
        else:
            deletion_size = len(ref) - len(variant)
        assert deletion_size > 0
        if deletion_size == 1:
            return "DEL1"
        elif deletion_size == 2:
            return "DEL2"
        elif deletion_size == 3:
            return "DEL3"
        else:
            return "DEL4"
    else:
        raise ValueError("mutation type is {} > {}, it must be either INS or DEL".format(ref, variant))

def INDEL_Alexandrov2018_16_category_name(row):
    return INDEL_Alexandrov2018_category_name_helper(row[COLNAME.FPRIME.value], row[COLNAME.REF.value], row[COLNAME.VAR.value], row[COLNAME.TPRIME.value])[0]

def INDEL_Alexandrov2018_83_category_name(row):
    return INDEL_Alexandrov2018_category_name_helper(row[COLNAME.FPRIME.value], row[COLNAME.REF.value], row[COLNAME.VAR.value], row[COLNAME.TPRIME.value])[1]

# Constants for Alexandrov2018 indel categories
ALEXANDROV_INDEL_RANGE_MAX = 5
ALEXANDROV_REPEAT_DEL_RANGE = range(0, 6)
ALEXANDROV_REPEAT_INS_RANGE = range(0, 6)
def INDEL_Alexandrov2018_category_name_helper(five_prime, ref, variant, three_prime):
    # Reverse complement if needed
    if ref in PURINES:
        ref = BASE_PAIR[ref]
        variant = reverse_complement(variant)
        five_prime_orig = five_prime
        five_prime	= reverse_complement(three_prime)
        three_prime	= reverse_complement(five_prime_orig)

    # TODO put this logic into its own function
    # TODO discuss - I'm not sure this logic is sound - see my code above
    if ref == '-' and len(variant) >= 1: # insertion
        # 5' and 3' flanking must be >= 5*len(variant) long to determine repeats
        if (len(five_prime) >= max(ALEXANDROV_REPEAT_INS_RANGE)*len(variant)) and (len(three_prime) >= max(ALEXANDROV_REPEAT_INS_RANGE)*len(variant)):
            indel_length = len(variant)
            n_repeat_units = min(ALEXANDROV_REPEAT_INS_RANGE)
            if len(variant) == 1: # 1bp insertion
                # do delayed reverse complement stuff
                if variant in PURINES:
                    variant = BASE_PAIR[variant]
                    five_prime_orig = five_prime
                    five_prime	= reverse_complement(three_prime)
                    three_prime	= reverse_complement(five_prime_orig)
                # determine homopolymer length (before the insertion)
                for flanking_bp in three_prime:
                    if flanking_bp == variant:
                        n_repeat_units += 1
                    else:
                        break
                for flanking_bp in five_prime[::-1]:
                    if flanking_bp == variant:
                        n_repeat_units += 1
                    else:
                        break
                n_repeat_units = min(max(ALEXANDROV_REPEAT_INS_RANGE), n_repeat_units)

                cat_name = 'INS_' + variant + '_1'
                subcat_name = cat_name + '_' + str(n_repeat_units)
                if n_repeat_units == max(ALEXANDROV_REPEAT_INS_RANGE):
                    subcat_name += '+'
            else: # >1bp insertion
                indel_length = min(ALEXANDROV_INDEL_RANGE_MAX, indel_length)
                # determine number of repeat units (before the insertion)
                n_repeat_units = min(ALEXANDROV_REPEAT_INS_RANGE)
                n_repeat_units_tprime = 0
                n_repeat_units_fprime = 0
                for three_i in ALEXANDROV_REPEAT_INS_RANGE:
                    if three_prime.startswith(variant * three_i):
                        n_repeat_units_tprime = three_i
                    else:
                        break
                for five_i in ALEXANDROV_REPEAT_INS_RANGE:
                    if five_prime.endswith(variant * five_i):
                        n_repeat_units_fprime = five_i
                    else:
                        break
                n_repeat_units += n_repeat_units_tprime + n_repeat_units_fprime
                n_repeat_units = min(max(ALEXANDROV_REPEAT_INS_RANGE), n_repeat_units)

                cat_name = 'INS_repeats_' + str(indel_length)
                if indel_length == ALEXANDROV_INDEL_RANGE_MAX:
                    cat_name += '+'
                subcat_name = cat_name + '_' + str(n_repeat_units)
                if n_repeat_units == max(ALEXANDROV_REPEAT_INS_RANGE):
                    subcat_name += '+'
        else:
            raise ValueError('Flanking base pair lengths too short for indel classification')

    # TODO - put this logic into its own function
    elif len(ref) >= 1 and variant == '-': # deletion
        # 5' and 3' flanking must be >= 6*len(ref) long to determine repeats
        if (len(five_prime) >= max(ALEXANDROV_REPEAT_DEL_RANGE)*len(ref)) and (len(three_prime) >= max(ALEXANDROV_REPEAT_DEL_RANGE)*len(ref)):
            indel_length = len(ref)
            n_repeat_units = min(ALEXANDROV_REPEAT_DEL_RANGE)
            if len(ref) == 1: # 1bp deletion
                # determine homopolymer length (before the deletion)
                for flanking_bp in three_prime:
                    if flanking_bp == ref:
                        n_repeat_units += 1
                    else:
                        break
                for flanking_bp in five_prime[::-1]:
                    if flanking_bp == ref:
                        n_repeat_units += 1
                    else:
                        break
                n_repeat_units = min(max(ALEXANDROV_REPEAT_DEL_RANGE), n_repeat_units)

                cat_name = 'DEL_' + ref + '_1'
                subcat_name = cat_name + '_' + str(n_repeat_units)
                if n_repeat_units == max(ALEXANDROV_REPEAT_DEL_RANGE):
                    subcat_name += '+'
            else: # >1bp deletion
                indel_length = min(ALEXANDROV_INDEL_RANGE_MAX, indel_length)
                # determine number of repeat units (before the deletion)
                n_repeat_units = min(ALEXANDROV_REPEAT_DEL_RANGE)
                n_repeat_units_tprime = 0
                n_repeat_units_fprime = 0
                for three_i in ALEXANDROV_REPEAT_DEL_RANGE:
                    if three_prime.startswith(ref * three_i):
                        n_repeat_units_tprime = three_i
                    else:
                        break
                for five_i in ALEXANDROV_REPEAT_DEL_RANGE:
                    if five_prime.endswith(ref * five_i):
                        n_repeat_units_fprime = five_i
                    else:
                        break
                n_repeat_units += n_repeat_units_tprime + n_repeat_units_fprime
                n_repeat_units = min(max(ALEXANDROV_REPEAT_DEL_RANGE), n_repeat_units)

                cat_name = 'DEL_repeats_' + str(indel_length)
                if indel_length == ALEXANDROV_INDEL_RANGE_MAX:
                    cat_name += '+'
                subcat_name = cat_name + '_' + str(n_repeat_units)
                if n_repeat_units == max(ALEXANDROV_REPEAT_DEL_RANGE):
                    subcat_name += '+'

                if n_repeat_units == min(ALEXANDROV_REPEAT_DEL_RANGE):
                    # Check for deletion with microhomology
                    n_overlap_tprime = 0
                    n_overlap_fprime = 0
                    # Check 3' overlap
                    for tp_bp, mut_bp in zip(list(three_prime[:len(ref)]), list(ref)):
                        if tp_bp == mut_bp:
                            n_overlap_tprime += 1
                        else:
                            break
                    # Check 5' overlap
                    for fp_bp, mut_bp in zip(list(five_prime[-len(ref):][::-1]), list(ref[::-1])):
                        if fp_bp == mut_bp:
                            n_overlap_fprime += 1
                        else:
                            break
                    n_overlap = max(n_overlap_tprime, n_overlap_fprime)
                    if n_overlap >= n_repeat_units:
                        # Falls into microhomology category, set name
                        cat_name = 'DEL_MH_' + str(indel_length)
                        if indel_length == 5:
                            cat_name += '+'
                        subcat_name = cat_name + '_' + str(n_overlap)
                        if n_overlap == 5:
                            subcat_name += '+'
        else:
            raise ValueError('Flanking base pair lengths too short for indel classification')
    return (cat_name, subcat_name)

# Lists
def INDEL_Alexandrov2018_16_category_list():
    cats = []
    # Insertions (1bp)
    cats += ['INS_C_1', 'INS_T_1']
    # Insertions (at repeats)
    cats += ['INS_repeats_2', 'INS_repeats_3', 'INS_repeats_4', 'INS_repeats_5+']
    # Deletions (1bp)
    cats += ['DEL_C_1', 'DEL_T_1']
    # Deletions (at repeats)
    cats += ['DEL_repeats_2', 'DEL_repeats_3', 'DEL_repeats_4', 'DEL_repeats_5+']
    # Deletions with microhomology
    cats += ['DEL_MH_2', 'DEL_MH_3', 'DEL_MH_4', 'DEL_MH_5+']
    return cats

def INDEL_Alexandrov2018_83_category_list():
    cats = []
    # Insertions (1bp)
    cats += ['INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+']
    cats += ['INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+']
    # Insertions (at repeats)
    cats += ['INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4', 'INS_repeats_2_5+']
    cats += ['INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3', 'INS_repeats_3_4', 'INS_repeats_3_5+']
    cats += ['INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2', 'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+']
    cats += ['INS_repeats_5+_0', 'INS_repeats_5+_1', 'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+']
    # Deletions (1bp)
    cats += ['DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+']
    cats += ['DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+']
    # Deletions (at repeats)
    cats += ['DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4', 'DEL_repeats_2_5+']
    cats += ['DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3', 'DEL_repeats_3_4', 'DEL_repeats_3_5+']
    cats += ['DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2', 'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+']
    cats += ['DEL_repeats_5+_0', 'DEL_repeats_5+_1', 'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+']
    # Deletions with microhomology
    cats += ['DEL_MH_2_1']
    cats += ['DEL_MH_3_1', 'DEL_MH_3_2']
    cats += ['DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3']
    cats += ['DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+']
    return cats


def INDEL_Haradhvala2018_8_category_list():
    return ["INS1", "INS2", "INS3", "INS4", "DEL1", "DEL2", "DEL3", "DEL4"]
