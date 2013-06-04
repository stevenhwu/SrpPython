from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import AlignInfo
from Bio.Phylo.BaseTree import Tree, Clade
from cStringIO import StringIO
from run_ext_prog import runExtProg
import Bio
import os
import re
import shutil


MATCH_TREE_NODE = re.compile("(\d+\.)(\d+)(\.\d+)")

dataDir = "/home/sw167/Postdoc/Project_Srp/simulateData/"
softwareDir = "/home/sw167/Postdoc/Project_Srp/software/"
metaSimDir = softwareDir + "metaSim/"
seqgenDir = softwareDir + "SeqGen/"
shrimpDir = softwareDir + "SHRiMP/bin/"
bsscDir = softwareDir + "BSSC/"


def runBSSC(workingFile, srp_tree_file):

    bssc = runExtProg(bsscDir + "BSSC_original", pdir=bsscDir, length=3)
    bssc.set_param_at("-f", 1)
    bssc.set_param_at(workingFile + ".par", 2)
    bssc.set_param_at(1, 3)
    bssc.run(0)

    bssc_result = workingFile + ".paup"
    while not os.path.exists(bssc_result):
        bssc.run()


    input_handle = open(bssc_result, "rU")
    sequences = AlignIO.read(input_handle, "nexus")
    input_handle.close()

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


    output_handle = open(workingFile + "_seqgen.phylip", "w")
    ref_handle = open(workingFile + ".cons", "w")

    output_handle.write("7 1200\n")
    for i, seq in enumerate(sequences):
        if i == 0:
            ref_handle.write(">%s\n%s\n" % ("Ref", seq.seq))
        output_handle.write("%s %s\n" % (seq.id, seq.seq))
    ref_handle.close()

    output_handle.write("1\n")
    Phylo.write(tree, output_handle, "newick")
    output_handle.close()

    output_handle = open(srp_tree_file, "w")
    Phylo.write(tree, output_handle, "newick")
    output_handle.close()


#    print os.path.exists(bssc_result)
#    print os.remove(bssc_result)
#    print os.path.exists(bssc_result)

def runSeqGen(workingFile, srp_hap_file):

    seqgen_infile = workingFile + "_seqgen.phylip"
    seqgen = runExtProg(seqgenDir + "./seq-gen", pdir=seqgenDir, length=3)
    seqgen.set_param_at("-mHKY", 1)
    seqgen.set_param_at("-t2", 2)
    seqgen.set_param_at("-k1", 3)
    seqgen.set_stdin(seqgen_infile)
    seqgen.run(0)


    temp_handle = open(workingFile + "_seqgen_out.phylip", "w")
    temp_handle.write(seqgen.output)
    temp_handle.close()

    SeqIO.convert(workingFile + "_seqgen_out.phylip", "phylip", workingFile + ".fasta", "fasta")
    shutil.copy(workingFile + ".fasta", srp_hap_file)

# ./seq-gen -mHKY -t2 -k1  < H7.fasta > zz.fasta



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
    metaSim.run(0)

# ./MetaSim cmd --454 -f 250 -t 25 -r300 -c -d /home/sw167/Postdoc/Project_Srp/simulateData/H6a /home/sw167/Postdoc/Project_Srp/simulateData/H6a/H6a.fasta


def runShrimp(workingFile, srp_file):
    shrimp_infile = workingFile + "-454.fna"
    reference_file = workingFile + ".cons"
    sam_outfile = workingFile + ".sam"
    srp_temp_file = workingFile + "_temp.fasta"


    shrimp = runExtProg(shrimpDir + "gmapper-ls", length=3)
    shrimp.set_param_at(shrimp_infile, 1)
    shrimp.set_param_at(reference_file, 2)
    shrimp.set_param_at("-h 30%", 3)
#    shrimp.set_param_at("-P", 3)
    shrimp.run(0)

    outfile_handle = open(sam_outfile, "w")
    outfile_handle.write(shrimp.output)
#    print shrimp.output


    sam2fasta = runExtProg(softwareDir + "./sam2fasta.py", length=3)
    sam2fasta.set_param_at(reference_file, 1)
    sam2fasta.set_param_at(sam_outfile, 2)
    sam2fasta.set_param_at(srp_temp_file, 3)
    sam2fasta.run()

    temp_handle = open(srp_temp_file, "rU")
    srp_handle = open(srp_file, "w")

    last_line = 0
    for i, line in enumerate(temp_handle):
        if i > 1:
            srp_handle.write(line)
        last_line = i
#    print last_line / 2
    temp_handle.close()
    srp_handle.close()


def main():

    prefix = "H7"
    workingDir = dataDir + prefix + "/"
    workingFile = workingDir + prefix
#    base_bssc_infile = workingFile + ".par"

    for index in range(100):
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


        runBSSC(workingFile, srp_tree_file)
        runSeqGen(workingFile, srp_hap_file)
        runMetaSim(workingFile, workingDir)
        runShrimp(workingFile, srp_file)

if __name__ == '__main__':
    main()





