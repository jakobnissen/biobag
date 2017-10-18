"""Reads a SAM file from stdin to stdout.
If both members of a mate pair maps to same gene, discards the first mate.
"""

import sys

for line in sys.stdin:
    # Write if line begins with @ (i.e. it is a header), always write
    if line.startswith('@'):
        sys.stdout.write(line)
        
    else:
        fields = line.split('\t')
    
        # If paired (1) and first mate (64), but neither unmapped (4) nor
        # mate unmapped (8), and mates mapped to same gene (=), do not write.
        if int(fields[1]) & 77 == 65 and fields[6] == '=':
            pass
        
        else:
            sys.stdout.write(line)
