from Bio import AlignIO, SeqIO, Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from cStringIO import StringIO
from run_ext_prog import runExtProg
import Bio
import os
import re

from Bio.Align import AlignInfo
import shutil

MATCH_TREE_NODE = re.compile("(\d+\.)(\d+)(\.\d+)")

dataDir = "/home/sw167/Postdoc/Project_Srp/simulateData/"
softwareDir = "/home/sw167/Postdoc/Project_Srp/software/"
metaSimDir = softwareDir + "metaSim/"
shrimpDir = softwareDir + "SHRiMP/bin/"
bsscDir = softwareDir + "BSSC/"


def runBSSC(workingFile, srp_hap_file, srp_tree_file):

    bssc = runExtProg(bsscDir + "BSSC_original", pdir=bsscDir, length=3)
    bssc.set_param_at("-f", 1)
    bssc.set_param_at(workingFile + ".par", 2)
    bssc.set_param_at(1, 3)
    bssc.run()

    bssc_result = workingFile + ".paup"
    while not os.path.exists(bssc_result):
        bssc.run()


    input_handle = open(bssc_result, "rU")
    output_handle = open(workingFile + ".fasta", "w")
    cons_handle = open(workingFile + ".cons", "w")

    sequences = AlignIO.read(input_handle, "nexus")

    summary_align = AlignInfo.SummaryInfo(sequences)
    consensus = summary_align.dumb_consensus(0.1, "N")
    cons_handle.write(">Cons\n")
    cons_handle.write(str(consensus))
    cons_handle.close()

    new_seqs = []
    for i, seq in enumerate(sequences):
        seq.id = "hap_" + str(i)
        new_seqs.append(seq)
    SeqIO.write(new_seqs, output_handle, "fasta")

    output_handle.close()
    input_handle.close()

    shutil.copy(workingFile + ".fasta", srp_hap_file)

    input_handle = open(bssc_result, "rU")
    for line in input_handle:
        if line.find("tree true_tree_1") > 0:
            line = line.strip()
            start = line.index("U]") + 3
            treeString = line[start:]
            tree = Phylo.read(StringIO(treeString), 'newick')
    input_handle.close()

    for clade in tree.find_clades():
        if clade.name:
            match = re.match(MATCH_TREE_NODE, clade.name)
            if match:
                index = match.group(2)
                clade.name = "hap_" + str(int(index) - 1)

    Phylo.write(tree, srp_tree_file, "newick")
#    print os.path.exists(bssc_result)
#    print os.remove(bssc_result)
#    print os.path.exists(bssc_result)

def runMetaSim(workingFile, workingDir):
    metaSim_infile = workingFile + ".fasta"

    if os.path.exists(workingFile + "-454.fna"):
        os.remove(workingFile + "-454.fna")


    metaSim = runExtProg(metaSimDir + "./MetaSim", pdir=metaSimDir, length=9)
    metaSim.set_param_at("cmd", 1)
    metaSim.set_param_at("--454", 2)
    metaSim.set_param_at("-f250", 3)
    metaSim.set_param_at("-t25", 4)
    metaSim.set_param_at("-r700", 5)
    metaSim.set_param_at("-c", 6)
    metaSim.set_param_at("-d", 7)
    metaSim.set_param_at(workingDir, 8)
    metaSim.set_param_at(metaSim_infile, 9)
    metaSim.run()

# ./MetaSim cmd --454 -f 250 -t 25 -r300 -c -d /home/sw167/Postdoc/Project_Srp/simulateData/H6a /home/sw167/Postdoc/Project_Srp/simulateData/H6a/H6a.fasta


def runShrimp(workingFile, srp_file):
    shrimp_infile = workingFile + "-454.fna"
    consensus_file = workingFile + ".cons"
    sam_outfile = workingFile + ".sam"
    srp_temp_file = workingFile + "_temp.fasta"


    shrimp = runExtProg(shrimpDir + "gmapper-ls", length=3)
    shrimp.set_param_at(shrimp_infile, 1)
    shrimp.set_param_at(consensus_file, 2)
    shrimp.set_param_at("-h 30%", 3)
#    shrimp.set_param_at("-P", 3)
    shrimp.run(0)

    outfile_handle = open(sam_outfile, "w")
    outfile_handle.write(shrimp.output)
#    print shrimp.output


    sam2fasta = runExtProg(softwareDir + "./sam2fasta.py", length=3)
    sam2fasta.set_param_at(consensus_file, 1)
    sam2fasta.set_param_at(sam_outfile, 2)
    sam2fasta.set_param_at(srp_temp_file, 3)
    sam2fasta.run()

    temp_handle = open(srp_temp_file, "rU")
    srp_handle = open(srp_file, "w")

    for i, line in enumerate(temp_handle):
        if i > 1:
            srp_handle.write(line)

    temp_handle.close()
    srp_handle.close()


def main():

    prefix = "H7"
    workingDir = dataDir + prefix + "/"
    workingFile = workingDir + prefix
#    base_bssc_infile = workingFile + ".par"

    for index in range(50, 100):
        str_index = str(index)
#        workingDir = dataDir + prefix + "/" + prefix + "_" + str_index + "/"
#        workingFile = workingDir + prefix + "_" + str_index

#        if not os.path.exists(workingDir):
#            os.mkdir(workingDir)
#        shutil.copy(base_bssc_infile, workingFile + ".par")

        resultDir = workingDir + prefix + "_" + str_index + "/"
        if not os.path.exists(resultDir):
            os.mkdir(resultDir)
        srp_file = resultDir + prefix + "Srp.fasta"
        srp_tree_file = resultDir + prefix + "Srp.tree"
        srp_hap_file = resultDir + prefix + "Srp_fullHaplotype.fasta"

        print index, resultDir


        runBSSC(workingFile, srp_hap_file, srp_tree_file)
        runMetaSim(workingFile, workingDir)
        runShrimp(workingFile, srp_file)

if __name__ == '__main__':
    main()





