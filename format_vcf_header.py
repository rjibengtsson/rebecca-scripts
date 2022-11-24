import sys

# Read vcf file
in_handle = list(open(sys.argv[1], 'r').read().splitlines())

# Read codex header list
codex_lst = list(open(sys.argv[2], 'r').read().splitlines())

# convert contig names and make new list
new_lst = []


# format contig names
for line in in_handle:
    # print(line)
    if line.startswith('##contig=<ID=03123017924237'):
        line_split = line.split('-')
        line_end = line.split(',')[-1]
        new_contig = str(line_split[0] + '-' + line_split[1] + '-' + line_split[2] + ',' + line_end + '\n')
        new_lst.append(new_contig)
    else:
        new_lst.append(line)


tmp_lst = []

# Check if contig name is in codex header list
for contig in new_lst:
    if contig.startswith('##contig=<ID='):
        contig_split = contig.split(',')[0]
        item = contig_split.replace('##contig=<ID=', '')
        if item in codex_lst:
            tmp_lst.append(contig)
    else:
        tmp_lst.append(contig)

tmp_lst = [x[:-1] for x in tmp_lst]

# write tmp_lst into a new file
with open(sys.argv[3], 'w') as fp:
    for i in tmp_lst:
        fp.write("%s\n" % i)