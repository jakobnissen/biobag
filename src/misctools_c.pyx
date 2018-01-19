import numpy as _np
import cython

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

@cython.boundscheck(False)
@cython.cdivision(True)
cpdef void tetnucfq_unsafe(unsigned char[:] bytesarray, long[:] counts, float[:] result):
    """Given an np.uin8 bytearray of only the ACGT sequence, a 136 int array,
    and a 136 np.float32 array to write results to, writes TNF to the latter."""
    
    cdef int fourmer = 0
    cdef int character
    cdef int charvalue
    cdef int i
    cdef int countdown = 3
    cdef int sumofcounts = 0
    cdef int contiglength = len(bytesarray)
    cdef unsigned char[:] localcomplementer = complementer
    
    # First clear out counts
    for i in range(136):
        counts[i] = 0

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
            countdown = 3
            fourmer = 0
            continue

        fourmer += charvalue
        
        if countdown == 0:
            counts[localcomplementer[fourmer]] += 1
            fourmer &= 0b111111
        
        else:
            countdown -= 1
            
        fourmer <<= 2
        
    for i in range(136):
        sumofcounts += counts[i]
    
    for i in range(136):
        result[i] = (<float>counts[i]) / sumofcounts

cpdef tetnucfq_safe(unsigned char[:] bytearray):
    counts = _np.zeros(136, dtype=_np.int)
    result = _np.empty(136, dtype=_np.float32)
    tetnucfq_unsafe(bytearray, counts, result)
    return result
