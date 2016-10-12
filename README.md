MoTeX: A word-based HPC tool for single and structured MoTif eXtraction
=====

<b>Description</b>: MoTeX is an accurate and efficient tool for single and structured MoTif eXtraction. It comes in three flavors: the standard CPU version; the OpenMP-based version; and the MPI-based version. It includes a tool that implements measures for assesing the statistical significance of the reported motifs. 

<b>Installation</b>: To compile MoTeX-II, please follow the instructions given in file INSTALL.

<b>Usage</b>:
```
m    m       mmmmmmm        m    m
##  ##  mmm     #     mmm    #  #
# ## # #" "#    #    #"  #    ##
# "" # #   #    #    #""""   m""m
#    # "#m#"    #    "#mm"  m"  "m

 Usage: motexCPU|motexOMP|motexMPI <options>
 Standard (Mandatory):
  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'
                                      for protein  sequences. You may use `USR'
                                      for user-defined alphabet; edit the file
                                      motexdefs.h accordingly.
  -i, --input-file          <str>     (Multi)FASTA input filename.
  -o, --output-file         <str>     MoTeX output filename.
  -d, --distance            <int>     The  distance  used  for extracting  the
                                      motifs. It can be  either 0 (for Hamming
                                      distance) or 1 (for edit distance).
  -k, --motifs-length       <int>     The length for motifs.
  -e, --errors              <int>     Limit the  max number  of errors to this
                                      value.
  -q, --quorum              <int>     The quorum is the minimum percentage (%)
                                      of sequences in which a motif must occur.

 Optional:
  -Q, --max-quorum          <int>     The maximum percentage (%) of sequences
                                      in which a motif can occur (default: 100).
  -n, --num-of-occ          <int>     The minimum  number of  occurrences of a
                                      reported  motif in any  of the sequences
                                      (default: 1).
  -N, --max-num-of-occ      <int>     The maximum  number of  occurrences of a
                                      reported  motif in any  of the sequences
                                      (default: 10000).
  -s, --structured-motifs   <str>     Input filename  for the structure of the
                                      boxes in the case of structured motifs.
  -S, --SMILE-out-file      <str>     SMILE-like output filename to be used by
                                      SMILE.
  -b, --background-in-file  <str>     MoTeX background filename for statistical
                                      evaluation passed as input.
  -t, --threads             <int>     Number of threads to be used by the OMP
                                      version (default: 4).
  -L, --long-sequences      <int>     If the number of input sequences is less
                                      than  the number of  processors  used by
                                      the MPI version, this should be set to 1
                                      (default: 0); useful  for a few (or one)
                                      very long sequence(s), e.g. a chromosome.
  -u, --un-out-file         <str>     Output filename for foreground motifs not
                                      matched exactly with any background motif
                                      in the file passed with the `-b' option.
  -I, --un-in-file          <str>     Input filename of the aforementioned file
                                      with the unmatched  motifs. These  motifs
                                      will be approximately searched  as motifs
                                      in the file passed with the `-i' option.
  -U, --SMILE-un-out-file   <str>     SMILE-like output filename for foreground
                                      motifs  not  matched  exactly  with  any
                                      background motif in the file passed with
                                      the `-b' option.
```
<b>Example</b>: For typical runs, see file EXAMPLES.

<b>Citation</b>:

```
S. P. Pissis, "MoTeX-II: structured MoTif eXtraction from large-scale datasets", BMC Bioinformatics, vol. 15, 2014, pp. 235.
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2012 Solon P. Pissis.
