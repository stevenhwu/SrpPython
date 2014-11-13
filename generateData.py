from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import AlignInfo
from Bio.Cluster import distancematrix
from Bio.Phylo.BaseTree import Tree, Clade
from cStringIO import StringIO
from run_ext_prog import runExtProg
import Bio
import numpy as np
import os
import re
import shutil

MATCH_TREE_NODE = re.compile("(\d+\.)(\d+)(\.\d+)")

projectDir = "/home/steven/Postdoc/Project_Srp/"

dataDir = projectDir + "simulateData/"
softwareDir = projectDir + "software/"
metaSimDir = softwareDir + "metaSim/"
artDir = softwareDir + ""
seqgenDir = softwareDir + "SeqGen/"
shrimpDir = softwareDir + "SHRiMP/bin/"
bsscDir = softwareDir + "BSSC/"


def runBSSC(workingFile, srp_tree_file, debug):

    bssc = runExtProg(bsscDir + "BSSC_original", pdir=bsscDir, length=3)
    bssc.set_param_at("-f", 1)
    bssc.set_param_at(workingFile + "_BSSC.par", 2)
    bssc.set_param_at(1, 3)

    bssc_paup_result = workingFile + "_BSSC.paup"
#     bssc_tree_result = workingFile + "_true_trees.trees"

    try:
        os.remove(bssc_paup_result)
#         os.remove(workingFile + "_0.pau")
    except OSError:
        pass

    while not os.path.exists(bssc_paup_result):
        bssc.run()

#     input_handle = open(bssc_tree_result, "rU")
    input_handle = open(bssc_paup_result, "rU")
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


    input_handle = open(bssc_paup_result, "rU")
    sequences = AlignIO.read(input_handle, "nexus")
    input_handle.close()
    seq = sequences[0]

    ref_handle = open(workingFile + ".cons", "w")
    ref_handle.write(">%s\n%s\n" % ("Ref", seq.seq))
    ref_handle.close()

    output_handle = open(workingFile + "_seqgen.phylip", "w")
    output_handle.write("1 1200\n")
    output_handle.write("%s %s\n" % ("ancestor", seq.seq))
    output_handle.write("1\n")
    Phylo.write(tree, output_handle, "newick")
    output_handle.close()

    output_handle = open(srp_tree_file, "w")
    Phylo.write(tree, output_handle, "newick")
    output_handle.close()


#    print os.path.exists(bssc_paup_result)
#    print os.remove(bssc_paup_result)
#    print os.path.exists(bssc_paup_result)

def runSeqGen(workingFile, srp_hap_file, srp_tree_file, debug):

    seqgen_infile = workingFile + "_seqgen.phylip"
    seqgen = runExtProg(seqgenDir + "./seq-gen", pdir=seqgenDir, length=3)
    seqgen.set_param_at("-mHKY", 1)
    seqgen.set_param_at("-t2", 2)
    seqgen.set_param_at("-k1", 3)
#     seqgen.set_param_at("-d0.1", 4)
    # seqgen.set_param_at("-s0.00001", 4)
    seqgen.set_stdin(seqgen_infile)

    all_unique = False
    repeat = 0
    while not all_unique:
        if repeat == 100:
            runBSSC(workingFile, srp_tree_file, debug)
            print "==========rerun BSSC========"
            repeat = 0
        repeat += 1
#        print repeat
        seqgen.run(0)
        all_unique = check_unique_sequences(seqgen)


    temp_handle = open(workingFile + "_seqgen_out.phylip", "w")
    temp_handle.write(seqgen.output)
    temp_handle.close()

    SeqIO.convert(workingFile + "_seqgen_out.phylip", "phylip", workingFile + ".fasta", "fasta")
    shutil.copy(workingFile + ".fasta", srp_hap_file)

# ./seq-gen -mHKY -t2 -k1  < H7.fasta > zz        runBSSC(workingFile, srp_tree_file).fasta


def check_unique_sequences(seqgen):
    # create our hash table to add the sequences
    sequences = {}
    data = AlignIO.read(StringIO(seqgen.output), "phylip")

    for seq_record in data:
        sequence = str(seq_record.seq).upper()
        if sequence not in sequences:
            sequences[sequence] = seq_record.id
        else:
            return False

    return True


def runMetaSim(workingFile, workingDir, debug):
    metaSim_infile = workingFile + ".fasta"

    if os.path.exists(workingFile + "-454.fna"):
        os.remove(workingFile + "-454.fna")


    metaSim = runExtProg(metaSimDir + "./MetaSim", pdir=metaSimDir, length=9)
    metaSim.set_param_at("cmd", 1)
    metaSim.set_param_at("--454", 2)
    metaSim.set_param_at("-f250", 3)
    metaSim.set_param_at("-t25", 4)
    metaSim.set_param_at("-r2000", 5)
    metaSim.set_param_at("-c", 6)
    metaSim.set_param_at("-d", 7)
    metaSim.set_param_at(workingDir, 8)
    metaSim.set_param_at(metaSim_infile, 9)
    metaSim.run(debug)

# ./MetaSim cmd --454 -f 250 -t 25 -r300 -c -d /home/sw167/Postdoc/Project_Srp/simulateData/H6a /home/sw167/Postdoc/Project_Srp/simulateData/H6a/H6a.fasta


def runART(workingFile, art_output_prefix, error_free_model=True, debug=False):
    art_infile = workingFile + ".fasta"

    art = runExtProg(artDir + "./art_illumina", pdir=artDir)
    art.set_param_at("-sam", 1)
    art.add_switch("-i %s" % art_infile)
    art.add_switch("-l 75")
    art.add_switch("-f 20")
    art.add_param("-o %s%s" % (workingFile, art_output_prefix))
    if error_free_model:
        art.add_switch("-ef")
#     art.add_param("-na")
#    print art.get_extract_switch()
#     art.set_param_at("-c", 6)
#     art.set_param_at("-d", 7)
#     art.set_param_at(workingDir, 8)
#     art.set_param_at(art_infile, 9)
    art.run(debug)

# art_illumina -sam -i reference.fa -l 50 -f 10 -o single_dat


def runSam2Fasta(workingFile, result_file_prefix, samfile_prefix, srp_read_prefix, error_free_model=True, debug=False):

    reference_file = workingFile + ".cons"
    srp_temp_file = workingFile + "_temp.fasta"
    sam_outfile = workingFile + samfile_prefix + ".sam"

#     sam2fasta = runExtProg(softwareDir + "./sam2fasta.py", length=3)
    sam2fasta = runExtProg("./sam2fasta_mod.py", length=3)
    sam2fasta.set_param_at(reference_file, 1)
    sam2fasta.set_param_at(sam_outfile, 2)
    sam2fasta.set_param_at(srp_temp_file, 3)
    sam2fasta.run(debug)


    shutil.copy(sam_outfile, result_file_prefix + samfile_prefix + ".sam")
    shutil.copy(srp_temp_file, result_file_prefix + srp_read_prefix + samfile_prefix + ".fasta")

    if error_free_model:
        sam_outfile = workingFile + samfile_prefix + "_errFree.sam"
        sam2fasta.set_param_at(sam_outfile, 2)
        sam2fasta.run(debug)

        shutil.copy(sam_outfile, result_file_prefix + samfile_prefix + "_errFree.sam")
        shutil.copy(srp_temp_file, result_file_prefix + srp_read_prefix + samfile_prefix + "_errFree.fasta")



def runShrimp(workingFile, shrimp_infile_prefix, shrimp_prefix, debug):
    reference_file = workingFile + ".cons"

    shrimp_infile = workingFile + shrimp_infile_prefix
    sam_outfile = workingFile + shrimp_prefix + ".sam"

    shrimp = runExtProg(shrimpDir + "gmapper-ls", length=3)
    shrimp.set_param_at(shrimp_infile, 1)
    shrimp.set_param_at(reference_file, 2)
    shrimp.set_param_at("-h 25%", 3)
    shrimp.add_switch("--qv-offset 33")  # # for ART output
#     shrimp.add_switch("--ignore-qvs")
#    shrimp.set_param_at("-P", 3)
    shrimp.run(debug)

    outfile_handle = open(sam_outfile, "w")
    outfile_handle.write(shrimp.output)

#
# def runSam2FastaShrimp(workingFile, samfile_prefix, srp_file, debug=False):
#
#     reference_file = workingFile + ".cons"
#     srp_temp_file = workingFile + "_temp.fasta"
#     sam_outfile = workingFile + samfile_prefix + ".sam"
#
# #     sam2fasta = runExtProg(softwareDir + "./sam2fasta.py", length=3)
#     sam2fasta = runExtProg("./sam2fasta_mod.py", length=3)
#     sam2fasta.set_param_at(reference_file, 1)
#     sam2fasta.set_param_at(sam_outfile, 2)
#     sam2fasta.set_param_at(srp_temp_file, 3)
#     sam2fasta.run(1)
#
#     temp_handle = open(srp_temp_file, "rU")
#     srp_handle = open(srp_file, "w")
# #     print reference_file, sam_outfile, srp_temp_file
#
#     last_line = 0
#     for i, line in enumerate(temp_handle):
#         if i > 1:
#             srp_handle.write(line)
#         last_line = i
# #    print last_line / 2
#     temp_handle.close()
#     srp_handle.close()



def main():

    prefix = "H10"
    working_directory = os.path.join(dataDir, prefix , "")
    working_file_prefix = os.path.join(working_directory, prefix)
#    base_bssc_infile = working_file_prefix + ".par"

    srp_read_prefix = "_ShortRead"
    art_prefix = "_ART"
    shrimp_prefix = "_Shrimp"

    print working_directory, working_file_prefix

    for index in range(0, 10):
        str_index = str(index)

#        if not os.path.exists(working_directory):
#            os.mkdir(working_directory)
#        shutil.copy(base_bssc_infile, working_file_prefix + ".par")
        prefix_index = prefix + "_" + str_index

        result_dir = os.path.join(working_directory, prefix_index, "")
        result_file_prefix = os.path.join(result_dir, prefix_index)
        print index, result_dir, result_file_prefix

        if not os.path.exists(result_dir):
            os.mkdir(result_dir)

        srp_tree_file = result_file_prefix + "_FullTree.tree"
        srp_hap_file = result_file_prefix + "_FullHaplotype.fasta"

#         srp_read_file_prefix = result_file_prefix + srp_read_prefix
        debug = 0

        error_free_model = True
        runBSSC(working_file_prefix, srp_tree_file, debug)
        runSeqGen(working_file_prefix, srp_hap_file, srp_tree_file, debug)
#
        runART(working_file_prefix, art_prefix, error_free_model, debug)
        runSam2Fasta(working_file_prefix, result_file_prefix, art_prefix, srp_read_prefix, error_free_model, debug)

        shrimp_infile = art_prefix + ".fq"
        runShrimp(working_file_prefix, shrimp_infile, shrimp_prefix, debug)
        error_free_model = False
        runSam2Fasta(working_file_prefix, result_file_prefix, shrimp_prefix, srp_read_prefix, error_free_model, debug)

##### old steps

#         runMetaSim(working_file_prefix, working_directory, debug)
#         shrimp_infile = "-454.fna"  # metaSim
#         runShrimp(working_file_prefix, shrimp_infile, "_metaSim", debug)

#         runBSSC(working_file_prefix, srp_tree_file, debug)
#         runSeqGen(working_file_prefix, srp_hap_file, srp_tree_file, debug)
#         runMetaSim(working_file_prefix, working_directory, debug)
#         runShrimp(working_file_prefix, srp_file, debug)



def mainFromAlignment():

    prefix = "H7"
    workingDir = dataDir + prefix + "/"
    workingFile = workingDir + prefix
#    base_bssc_infile = workingFile + ".par"

#    for index in [12, 18, 19, 23, 26, 33, 36, 41, 43, 54, 60, 67, 68, 70, 73, 75, 89, 93]:  #
    for index in range(59, 60):
        str_index = str(index)
#        workingDir = dataDir + prefix + "/" + prefix + "_" + str_index + "/"
#        workingFile = workingDir + prefix + "_" + str_index

#        if not os.path.exists(workingDir):
#            os.mkdir(workingDir)
#        shutil.copy(base_bssc_infile, workingFile + ".par")
        prefix_index = prefix + "_" + str_index
        resultDir = workingDir + prefix_index + "/"
#         if not os.path.exists(resultDir):
#             os.mkdir(resultDir)
        workingDir = "/home/sw167/Postdoc/Project_Srp/simulateData/H7/H7_" + str_index + "/"
        workingFile = workingDir + "H7_" + str_index + "_Srp_fullHaplotype"
        srp_file = workingDir + "H7_" + str_index + "_Srp50x.fasta"
#         srp_tree_file = resultDir + prefix_index + "_Srp.tree"
#         srp_hap_file = resultDir + prefix_index + "_Srp_fullHaplotype.fasta"

#         print index, resultDir

        debug = 1
#         runBSSC(workingFile, srp_tree_file, debug)
#         runSeqGen(workingFile, srp_hap_file, srp_tree_file, debug)
        input_handle = open(workingFile + ".fasta", "rU")
        sequences = AlignIO.read(input_handle, "fasta")
        input_handle.close()
        seq = sequences[0]

        ref_handle = open(workingFile + ".cons", "w")
        ref_handle.write(">%s\n%s\n" % ("Ref", seq.seq))
        ref_handle.close()

        runMetaSim(workingFile, workingDir, debug)
        runShrimp(workingFile, srp_file, debug)

if __name__ == '__main__':
    main()
#     mainFromAlignment()
#     workingFile = "/home/sw167/Postdoc/Project_Srp/test/H1_100bp"
#     srp_file = workingFile + "out.fasta"
#     runShrimp(workingFile, srp_file, 1);





