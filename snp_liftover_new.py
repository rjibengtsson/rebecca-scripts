import argparse
import os
import re
import shutil

import gfapy


'''
Author: Benedict
So long, and thanks for all the fish!

Desc: Identify base differences between a given input strain reference and the mapped poistion on a template
Usage: python -m HelperScripts.snp_liftover -infolder [ROOT_TEMPLATE_OUTPUT_FOLDER] -minstrains [MIN_STRAINS] -minlength [MIN_LENGTH] -outfolder [FOLDER_FOR_OUTPUT]
Output: One VCF per input strain into the SNP_LIFTOVER subdirectory of the root template directory
NB: Does not consider INDELS at this time
'''

parser = argparse.ArgumentParser(description='Filter blocks')
parser.add_argument('-infolder','--infolder', help='Root folder for template pipeline run', required=True)
parser.add_argument('-minlength','--minlength', help='Min length cutoff', required=True)
parser.add_argument('-minstrains','--minstrains', help='Min no.strains cutoff', required=True)
parser.add_argument('-outfolder','--outfolder', help='Output folder vcfs', required=True)
args = vars(parser.parse_args())

'''
Given a multiple sequence alignment, extract the individual sequence names and their lengths.
'''
def get_block_info(block_filename:str, strain=None):
    info = []
    block_filename = '/scratch/rjibengtsson/templates/Paeruginosa_g1-V0.3.14/Template/' + block_filename.split('Template')[-1]
    with open(block_filename,'r') as blockfile:
        #print(blockfile)
        #exit(0)
        for line in blockfile:
            if line.startswith(">") and (strain is None or strain in line):
                seq_range = line.split('>')[1].rstrip()
                if seq_range[-1] == '+' or seq_range[-1] == '-':
                    name, seq_range = seq_range[:-1].split(':')
                else:
                    name, seq_range = seq_range.split(':')
                seq_range = seq_range.split('/')[0].split('-')
                lower = abs(int(seq_range[0]))
                upper = abs(int(seq_range[1]))
                length = upper - lower + 1
                info.append({'strain':name, 'length':length, 'start':lower, 'end':upper, 'name': block_filename})
    
    return info if strain is None else info[0]


def _get_fasta_length(fasta_file_path:str):
    sequence_str = ''
    with open(fasta_file_path, 'r') as fastafile:
        _ = fastafile.readline()
        for line in fastafile.readlines():
            if line.startswith('>'):
                continue
            sequence_str += ''.join(line.rstrip())
    return len(sequence_str)

def _read_fasta(fasta_file_path:str):
    sequence_str = ''
    with open(fasta_file_path, 'r') as fastafile:
        header = fastafile.readline().rstrip()
        for line in fastafile.readlines():
            sequence_str += ''.join(line.rstrip())
    return sequence_str, header

def _read_seq_from_block(fname:str, strain:str):
    sequence_str = ''
    start = False
    with open(fname, 'r') as fastafile:
        for line in fastafile.readlines():
            if strain in line:
                start = True
                continue
            elif start and '>' in line:
                return sequence_str 
            elif start:
                sequence_str += ''.join(line.rstrip())
    return sequence_str

def _write_fasta(sequence:str, header:str, filename:str, append=False, wrap=False):
    open_for = 'a' if append else 'w'
        
    with open(filename, open_for) as outfile:
        if header:
            outfile.write(header + "\n") 
            
        if wrap:
            for i in range(0, len(sequence), 70):
                if i + 70 > len(sequence) - 1:
                    seq_to_write = sequence[i:len(sequence)]
                else:
                    seq_to_write = sequence[i:i+70]
                    
                outfile.write(seq_to_write + "\n")
        else:
            outfile.write(sequence)

def _get_blocks_ordered_by_strain(gfa_graph):
    strains = {}
    for path in gfa_graph.paths:
        strain_name = str(path).split('\t')[1]
        list_of_links = []

        for link in path.links:
            list_of_links.append(str(link.line).split('\t')[1])
            list_of_links.append(str(link.line).split('\t')[3])

        strains[strain_name] = [get_block_info(block, strain=strain_name) for block in list(dict.fromkeys(list_of_links))]
    return strains

def _get_M_regions(mapping:str, strain:str):
    regions = re.split('(?<=M)|(?<=D)|(?<=N)', mapping)[:-1]
    offset_strain = 0
    offset_codex = 0
    M_regions_strain = []
    M_regions_codex = []
    for region in regions:
        if region[-1] == 'M':
            M_regions_codex.append((offset_codex, offset_codex + int(region[:-1])))
            M_regions_strain.append((offset_strain, offset_strain + int(region[:-1])))
            offset_strain += int(region[:-1])
            offset_codex += int(region[:-1])
        if region[-1] == 'N':
            offset_codex += int(region[:-1])
        if region[-1] == 'D':
            offset_strain += int(region[:-1])

    return M_regions_strain, M_regions_codex

def _get_mapping(block_name:str, strain:str):

    sequence_str = ''
    with open(block_name, 'r') as mapfile:
        start = False
        for line in mapfile.readlines():
            line = line.rstrip()
            if line.startswith('>') and strain in line:
                start = True
                continue
            elif line.startswith('>') and strain not in line: 
                start = False
                continue
            if start:
                sequence_str += ''.join(line.rstrip())

    return sequence_str

def _get_sequence(file):
    data = ''
    with open(file, 'r') as fastafile:
        for line in fastafile:
            if line.startswith('>'):
                continue
            data += line.replace('\n', '')
    return data

def _call_snps(block_name:str, strain_seq:str, codex_seq:str, codex_region_start:int):
    snps = []
    idx = 0

    for strain_base, codex_base in zip(strain_seq, codex_seq):
        if strain_base != codex_base:   
            snps.append((block_name, codex_base, strain_base, codex_region_start + idx + 1))
        idx += 1
    return snps

def _build_VCF(strain:str, snps:list, outdir:str):
    '''
    Construct a VCF from variants identified in liftover
    '''
    outfile = os.path.join(outdir, f"{strain}.vcf")
    with open(outfile, 'w') as vcf_file:
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for snp in snps:
            id = '.'
            vcf_file.write(f'{snp[0]}\t{snp[3]}\t{id}\t{snp[1]}\t{snp[2]}\n')

def _build_template(blocks_for_template:list, outdir:str):
    '''
    Construct a multi-contig fasta from blocks passing n_strain and n_length filters
    '''
    template_fasta_path = os.path.join(outdir, 'template_snp_liftover.fasta')
    for block in blocks_for_template:
        parts = block['name'].split('/')[-1].split('-')
        header = "{}-{}-{}".format(parts[0], parts[1], parts[2])
        sequence, _ = _read_fasta(block['name'] + '.mafft.cons.trimmed')
        _write_fasta(sequence, f'>{header}', template_fasta_path, append = True, wrap = True)


def get_blocks(root_dir):
    blocks = []
    strains = []
    
    gfa_graph = gfapy.Gfa.from_file(os.path.join(root_dir, 'template.gfa'))
    
    for path in gfa_graph.paths:
        strain_name = str(path).split('\t')[1]
        #print(strain_name)
        if strain_name not in strains:
            strains.append(strain_name)
        list_of_links = []

        for link in path.links:
            list_of_links.append(str(link.line).split('\t')[1])
            list_of_links.append(str(link.line).split('\t')[3])

        [blocks.append(get_block_info(block)) for block in list(dict.fromkeys(list_of_links)) if get_block_info(block) not in blocks]
    return strains, blocks

def main():

    ROOT_DIR = args.get('infolder', '')
    MIN_NUM_STRAINS = int(args.get('minstrains', 3))
    MIN_LENGTH = int(args.get('minlength', 500))
    OUTPUT_FOLDER = args.get('outfolder', 'my_outfolder')
    
    #if os.path.exists(OUTPUT_FOLDER):
        #shutil.rmtree(OUTPUT_FOLDER)
    os.mkdir(OUTPUT_FOLDER)
    
    # gfa_graph = gfapy.Gfa.from_file(os.path.join(ROOT_DIR, 'template.gfa'))
    # blocks_ordered_by_strain = _get_blocks_ordered_by_strain(gfa_graph)
    # blocks_for_template = []
    print("Load blocks")
    strains, blocks = get_blocks(ROOT_DIR)
    for strain in strains:
        snps = []
        print("strain: ", strain)
        for block in blocks:
            #print("Block: ", block)
            #for i in block:
            #    print(i)
            strain_info_for_block = [i for i in block if i.get('strain') == strain]
            if len(strain_info_for_block) == 0:
                continue
            else:
                strain_info_for_block = strain_info_for_block[0]
            
            block_name = strain_info_for_block.get('name')
            mapping = _get_mapping(f"{block_name}.mafft.cons.trimmed.map", strain)
            m_regions_strain, m_regions_codex = _get_M_regions(mapping, strain)
            m_regions = list(zip(m_regions_strain, m_regions_codex))


            for m_region in m_regions:
                codex_region = (m_region[1][0], m_region[1][1])
                codex_seq = _get_sequence(f'{block_name}.mafft.cons.trimmed')[codex_region[0] : codex_region[1]]

                strain_seq = _read_seq_from_block(block_name, strain)
                strain_seq = strain_seq[m_region[0][0] : m_region[0][1]]
                
                snps.extend(_call_snps(block_name, strain_seq, codex_seq, m_region[1][0]))
            
        _build_VCF(strain, snps, OUTPUT_FOLDER)

if __name__ == '__main__':
    main()
