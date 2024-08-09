for cfile in Genomes/*.fna ; do echo $cfile && /media/thomashitch/MyDATA/Projects/MIMIC/Gapseq/gapseq/gapseq doall $cfile && gzip $cfile; done
