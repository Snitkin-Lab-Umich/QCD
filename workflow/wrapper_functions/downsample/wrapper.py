__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from snakemake.shell import shell
import os
import sys

genome_size = snakemake.params.get("gsize", "")

file = {snakemake.input.r1}
file2 = {snakemake.input.r2}
outr1={snakemake.output.outr1}
outr2={snakemake.output.outr2}

def downsample_reads(file, file2, out1, out2, genome_size):
    file = file.pop()
    file2 = file2.pop()
    out1 = out1.pop()
    out2 = out2.pop()

    # Extract basic fastq reads stats with seqtk

    gsize = genome_size

    print("Using Genome Size: %s to calculate coverage" % gsize)
    
    seqtk_check = "seqtk fqchk -q3 %s > %s_fastqchk.txt" % (file, file)

    print(seqtk_check)

    try:
        os.system(seqtk_check)
    except sp.CalledProcessError:
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % file, 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    print('Original Covarage Depth: %s x' % ori_coverage_depth)

    # proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    # (nproc, err) = proc.communicate()
    # nproc = nproc.strip()

    if ori_coverage_depth >= 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        # r1_sub = "/tmp/%s" % os.path.basename(file)
        r1_sub = out1

        # Downsample using seqtk
        try:
            #print("Generating seqtk Downsampling command")
            print("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file, factor, r1_sub))

            seqtk_downsample = "seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (
                file, factor, r1_sub)
            os.system(seqtk_downsample)
            #call(seqtk_downsample, logger)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if file2:
            # r2_sub = "/tmp/%s" % os.path.basename(file2)
            r2_sub = out2
            try:
                print("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                os.system("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                #call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        os.system("cp %s %s" % (file, out1))
        os.system("cp %s %s" % (file2, out2))

downsample_reads(file, file2, outr1, outr2, genome_size)