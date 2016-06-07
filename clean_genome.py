import re

numbered_line_with_bases     = re.compile(r'^( +)\d+ [a-zA-Z ]+$')
expected_non_base_characters = re.compile(r'[ 0-9\n]')
not_valid_base               = re.compile(r'[^acgtrykmswbdhvn]')

with open('staph_aureus_genome_raw.txt',   'r') as raw_genome_file, \
     open('staph_aureus_genome_clean.txt', 'w') as clean_genome_file:

    for line in raw_genome_file:

        if not re.match(numbered_line_with_bases, line):
            continue

        bases = re.sub(expected_non_base_characters, '', line)

        if re.search(not_valid_base, bases):
            raise Exception('Invalid base: %s' % line.strip())

        clean_genome_file.write(bases)
