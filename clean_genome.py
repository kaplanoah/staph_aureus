import re

genome = ''

numbered_line_with_bases = re.compile(r'^( +)\d+ [a-z ]+$')
invalid_bases = re.compile(r'[^acgtrykmswbdhvn]')

with open('staph_aureus_genome_raw.txt',   'r') as raw_genome_file, \
     open('staph_aureus_genome_clean.txt', 'w') as clean_genome_file:

    for line in raw_genome_file:
        if re.match(numbered_line_with_bases, line):

            bases = re.sub(invalid_bases, '', line)

            if len(bases) != 60:
                print 'WARNING: unexpected number of valid bases'
                print line.strip()
