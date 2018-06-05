#cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False

import numpy as _np
import cython

# Notes/assumptions for nucleotide counters:
# 1) Contigs ONLY consists of values 64, 67, 71, 84, NOTHING ELSE
# 2) They are at least 4 bp long
# 3) A pseudocount of 1 is added to all unseen Kmers for Markov normalization

# Lookup in this array gives the index of the canonical tetranucleotide.
# E.g CCTA is the 92nd alphabetic 4mer, whose reverse complement, TAGG, is the 202nd.
# So the 92th and 202th value in this array is the same.
# Hence we can map 256 4mers to 136 normal OR reverse-complemented ones
cdef unsigned char[:] complementer = _np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 11, 31, 32,
        33, 34, 35, 36, 37, 38, 39, 40, 41, 23, 42, 43, 44, 7, 45, 46,
        47, 48, 49, 50, 51, 34, 52, 53, 54, 19, 55, 56, 57, 3, 58, 59,
        60, 57, 61, 62, 63, 44, 64, 65, 66, 30, 67, 68, 69, 14, 70, 71,
        72, 54, 73, 74, 75, 41, 76, 77, 78, 26, 79, 80, 66, 10, 81, 82,
        83, 51, 84, 85, 86, 37, 87, 88, 75, 22, 89, 90, 63, 6, 91, 92,
        93, 47, 94, 95, 83, 33, 96, 97, 72, 18, 98, 99, 60, 2, 100,
        101, 99, 56, 102, 103, 90, 43, 104, 105, 80, 29, 106, 107, 68,
        13, 108, 109, 97, 53, 110, 111, 88, 40, 112, 113, 77, 25, 114,
        105, 65, 9, 115, 116, 95, 50, 117, 118, 85, 36, 119, 111, 74,
        21, 120, 103, 62, 5, 121, 122, 92, 46, 123, 116, 82, 32, 124,
        109, 71, 17, 125, 101, 59, 1, 126, 125, 98, 55, 127, 120, 89,
        42, 128, 114, 79, 28, 129, 106, 67, 12, 130, 124, 96, 52, 131,
        119, 87, 39, 132, 112, 76, 24, 128, 104, 64, 8, 133, 123, 94,
        49, 134, 117, 84, 35, 131, 110, 73, 20, 127, 102, 61, 4, 135,
        121, 91, 45, 133, 115, 81, 31, 130, 108, 70, 16, 126, 100, 58, 0], dtype=_np.uint8)

# This is the same, but for both 4mers, 3mers and 2mers
cdef unsigned char[:] complementer2 = _np.array([  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
        13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
        26,  27,  28,  29,  30,  11,  31,  32,  33,  34,  35,  36,  37,
        38,  39,  40,  41,  23,  42,  43,  44,   7,  45,  46,  47,  48,
        49,  50,  51,  34,  52,  53,  54,  19,  55,  56,  57,   3,  58,
        59,  60,  57,  61,  62,  63,  44,  64,  65,  66,  30,  67,  68,
        69,  14,  70,  71,  72,  54,  73,  74,  75,  41,  76,  77,  78,
        26,  79,  80,  66,  10,  81,  82,  83,  51,  84,  85,  86,  37,
        87,  88,  75,  22,  89,  90,  63,   6,  91,  92,  93,  47,  94,
        95,  83,  33,  96,  97,  72,  18,  98,  99,  60,   2, 100, 101,
        99,  56, 102, 103,  90,  43, 104, 105,  80,  29, 106, 107,  68,
        13, 108, 109,  97,  53, 110, 111,  88,  40, 112, 113,  77,  25,
       114, 105,  65,   9, 115, 116,  95,  50, 117, 118,  85,  36, 119,
       111,  74,  21, 120, 103,  62,   5, 121, 122,  92,  46, 123, 116,
        82,  32, 124, 109,  71,  17, 125, 101,  59,   1, 126, 125,  98,
        55, 127, 120,  89,  42, 128, 114,  79,  28, 129, 106,  67,  12,
       130, 124,  96,  52, 131, 119,  87,  39, 132, 112,  76,  24, 128,
       104,  64,   8, 133, 123,  94,  49, 134, 117,  84,  35, 131, 110,
        73,  20, 127, 102,  61,   4, 135, 121,  91,  45, 133, 115,  81,
        31, 130, 108,  70,  16, 126, 100,  58,   0, 136, 137, 138, 139,
       140, 141, 142, 143, 144, 145, 146, 143, 147, 148, 149, 139, 150,
       151, 152, 149, 153, 154, 155, 146, 156, 157, 155, 142, 158, 159,
       152, 138, 160, 161, 159, 148, 162, 163, 157, 145, 164, 163, 154,
       141, 165, 161, 151, 137, 166, 165, 158, 147, 167, 164, 156, 144,
       167, 162, 153, 140, 166, 160, 150, 136, 168, 169, 170, 171, 172,
       173, 174, 170, 175, 176, 173, 169, 177, 175, 172, 168], dtype=_np.uint8)

cdef int reverse_complement_fourmer(int num):
    """Reverse-complements a fourmer.
    The input and output fourmers are represented by their alphabetic order.
    E.g. 92 (CCTA is 92nd 4mer) -> 202 (TAGG, the 202nd)"""
    
    cdef int result = 0
    cdef int inverse = num ^ 0b11111111
    
    # Here, from rightmost base to leftmose base, invert it (complement it),
    # Then move it to the right position. Each base is added independently to result
    result += (inverse & 0b11) << 6
    result += ((inverse >> 2) & 0b11) << 4
    result += ((inverse >> 4) & 0b11) << 2
    result += (inverse >> 6) & 0b11
    
    return result

cdef void tetnuccounts(unsigned char[:] bytesarray, int[:] counts):
    """Count tetranucleotides of contig and put them in counts vector.
    
    The bytearray is expected to be np.uint8 of bytevalues of the contig.
    The counts is expected to be an array of 256 32-bit integers with value 0"""
    
    cdef int kmer = 0
    cdef int character
    cdef int charvalue
    cdef int i
    cdef int countdown = 3
    cdef int contiglength = len(bytesarray)
            
    for i in range(contiglength):
        character = bytesarray[i]

        if character == 65: # A = 0b00 = 0
            charvalue = 0
        elif character == 67: # C = 0b01 = 1
            charvalue = 1
        elif character == 71: # G = 0b10 = 2
            charvalue = 2
        elif character == 84: # T = 0b11 = 3
            charvalue = 3
        else:
            kmer = 0
            countdown = 3
            continue
        
        kmer += charvalue
        
        if countdown == 0:
            counts[kmer] += 1
            kmer &= 0b111111 # Remove leftmost base 
        
        else:
            countdown -= 1
            
        kmer <<= 2 # Shift to prepare for next base

cdef void _tnf(unsigned char[:] bytearray, int[:] counts, float[:] result):
    """Puts TNF freqs in the results vector.
    
    The result is expected to be a 136 32-bit float vector
    The counts is expected to be an array of 256 32-bit integers
    """
    
    cdef int countsum = 0
    cdef int i
    cdef unsigned char[:] converter = complementer
    
    for i in range(256):
        counts[i] = 0
    
    tetnuccounts(bytearray, counts)
    
    for i in range(256):
        countsum += counts[i]
        
    cdef float floatsum = <float>countsum
        
    for i in range(256):
        result[converter[i]] += counts[i] / floatsum

# Assining these arrays for each sequence takes about 6% longer time than
# having assigned them once in userspace. Worth it.
cpdef tnf(unsigned char[:] sequence):
    """Returns 136-length float32 vector of TNF of DNA sequences.
    The DNA must be uppercase, represented as np.uint8 vector"""
    
    cdef int[:] counts = _np.empty(256, dtype=_np.int32)
    result = _np.zeros(136, dtype=_np.float32)
    cdef float[:] resultmem = result
    
    _tnf(sequence, counts, resultmem)
    return result
        
cdef void kmer432counts(unsigned char[:] bytesarray, int[:] counts):
    """Count tetranucleotides of contig and put them in counts vector.
    
    The bytearray is expected to be np.uint8 of bytevalues of the contig.
    The counts is expected to be an array of 336 32-bit integers with value 0.
    The first 256 is 4-mers, the next 64 is 3-mers and the last 16 2-mers."""
    
    cdef int kmer = 0
    cdef int character
    cdef int charvalue
    cdef int i
    cdef int countdown = 3
    cdef int contiglength = len(bytesarray)
        
    for i in range(contiglength):
        character = bytesarray[i]

        if character == 65:
            charvalue = 0
        elif character == 67:
            charvalue = 1
        elif character == 71:
            charvalue = 2
        elif character == 84:
            charvalue = 3
        else:
            kmer = 0
            countdown = 3
            continue
            
        kmer += charvalue
        
        if countdown < 3:
            counts[(kmer & 0b1111) + 256 + 64] += 1
        if countdown < 2:
            counts[(kmer & 0b111111) + 256] += 1
        if countdown < 1:
            counts[kmer] += 1
            kmer &= 0b111111 # Remove leftmost base
            countdown += 1
            
        countdown -= 1
        
        kmer <<= 2 # Shift to prepare for next base
        
cdef void _freq_432mers(unsigned char[:] bytearray, int[:] counts, float[:] result):
    cdef unsigned char[:] converter = complementer2
    cdef float fourfactor = 0
    cdef float threefactor = 0
    cdef float twofactor = 0
    cdef int i
    
    kmer432counts(bytearray, counts)
    
    for i in range(336):
        if i < 256:
            fourfactor += counts[i]
        elif i < 256 + 64:
            threefactor += counts[i]
        else:
            twofactor += counts[i]
            
    if twofactor != 0:
        twofactor = 1 / twofactor
    if threefactor != 0:
        threefactor = 1 / threefactor
    if fourfactor != 0:
        fourfactor = 1 / fourfactor

    for i in range(336):
        if i < 256:
            result[converter[i]] += counts[i] * fourfactor
        elif i < 256 + 64:
            result[converter[i]] += counts[i] * threefactor
        else:
            result[converter[i]] += counts[i] * twofactor
            
cpdef freq_432_mers(unsigned char[:] sequence):
    """Returns 136-length float32 vector of 4, 3, and 2 mers of DNA sequences.
    The DNA must be uppercase, represented as np.uint8 vector"""
    
    cdef int[:] counts = _np.zeros(336, dtype=_np.int32)
    result = _np.zeros(178, dtype=_np.float32)
    cdef float[:] resultmem = result
    
    _freq_432mers(sequence, counts, resultmem)
    return result
        
cdef void _markov_normalized_tnf(unsigned char[:] bytearray, int[:] counts, float[:] frequencies, float[:] result):
    """Puts TNF freqs in the results vector.
    
    The result is expected to be a 136 32-bit float vector
    The counts is expected to be an array of 336 32-bit integers
    The frequencues is expected to be an array of 336 32-bit floats
    For counts and frequencies, the first 256 is 4-mers, the next 64 is 3-mers and the last 16 2-mers.
    """
    
    cdef int foursum = 0
    cdef int threesum = 0
    cdef int twosum = 0
    cdef int i
    cdef int reverse
    cdef float value
    cdef unsigned char[:] converter = complementer
    
    for i in range(336):
        counts[i] = 0
        frequencies[i] = 0.0
    
    kmer432counts(bytearray, counts)

    for i in range(336):
        if counts[i] == 0: # prevent division by zero
            counts[i] = 1 
        
        if i >= 256+64:
            twosum += counts[i]
        elif i >= 256:
            threesum += counts[i]
        else:
            foursum += counts[i]
    
    # Invert here so when we normalize we do multiply (cheap) not divide (expensive)
    cdef float fourfactor = 1 / <float>foursum
    cdef float threefactor = 1 / <float>threesum
    cdef float twofactor = 1 / <float>twosum

    # Normalize by counts here  (the factors are inverted above)
    for i in range(336):
        if i >= 256+64:
            frequencies[i] = counts[i] * twofactor
        elif i >= 256:
            frequencies[i] = counts[i] * threefactor
        else:
            frequencies[i] = counts[i] * fourfactor
            
    for i in range(256):
        reverse = reverse_complement_fourmer(i)            
        
        # This is the actual Markov normalization:
        # 1: f(n1n2n3n4)
        value = frequencies[i] + frequencies[reverse]
        # 2: Times f(n2n3)
        value *= frequencies[((i & 0b111100) >> 2) + 320] + frequencies[((reverse & 0b111100) >> 2) + 320]
        # 3: Divided by f(n1n2n3)
        value /= frequencies[((i & 0b11111100) >> 2) + 256] + frequencies[((reverse & 0b11111100) >> 2) + 256]
        # 4: Divided by f(n2n3n4)
        value /= frequencies[(i & 0b111111) + 256] + frequencies[(reverse & 0b111111) + 256]
        
        result[converter[i]] += value

# Only about 14% slower than function tnf(sequence)
cpdef markov_normalized_tnf(unsigned char[:] sequence):
    """Returns 136-length float32 vector of TNF of DNA sequences.
    The DNA must be uppercase, represented as np.uint8 vector"""
    
    cdef int[:] counts = _np.empty(336, dtype=_np.int32)
    cdef float[:] frequencies = _np.empty(336, dtype=_np.float32)
    result = _np.zeros(136, dtype=_np.float32)
    cdef float[:] resultmem = result
    
    _markov_normalized_tnf(sequence, counts, frequencies, resultmem)
    return result