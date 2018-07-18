gc_plotter
=========

## Simple tool to map GC bins on chromosome ideograms with existing tools

Inspired by (and code taken from) [CMplot](https://github.com/YinLiLin/R-CMplot)
<br>
The Unix commands for creating GC content tracks also have been taken from [here](https://wiki.bits.vib.be/index.php/Create_a_GC_content_track)

---

### Run

# You must have bedtools, samtools and basic unix utils installed

Simply change the script to executable

```bash
chmod 755 gc.sh
```

Then run the script with 2 arguments, first is the reference genome fasta and the second is the desired width to scan for GC percenatge. 
```bash
./gc /home/reference.fasta 1000
```
The following files will be created in the directory which has the reference genome

".sizes"       file which specifies all the chromosomes/scaffolds and the lengths in basepairs
".bed"         a bed file
".fai"         fasta index from samtools
".igv"         igv viewer file with positions, bins, and gc content. Output file needed for gc_plotter 

---

### Plots

Two simple examples are the plasmodium species *P. falciparum* and *P. knowlesi*

*P. falsiparum* is well-known as a very [AT rich](https://genomevolution.org/wiki/index.php/Plasmodia_comparative_genomics) and GC poor genome. While *P. knowlesi* is still AT rich, but much less rich than *falciparum*. We can plot these and see if we get the expected results.

The genomes were downloaded from [ncbi](https://www.ncbi.nlm.nih.gov/genome/?term=txid1245013[Organism:noexp]) and used as is, except for a renaming. 

```bash
./gc.sh /home/plasmodium_knowlesi_genomic.fna 10000
```
<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/plas1.jpg">
<img src="plas1.jpg" height="800px" width="800px">
</a>
</p>

---
Same for the second genome:

```bash
./gc.sh /home/plasodium_falciparum_genomic.fna 10000
```
<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/plas.jpg">
<img src="plas.jpg" height="800px" width="800px">
</a>
</p>

As we can see, there is lots of red for *falicaprum* indicating that yes, it seems to be well within the 20's (23%), and *knowlesi* is around 40% as expected. 

---

### Hints

The legend can be freely changed, simply replace "topright" or "bottomright" on the last line in gcplot.r to another location or coordinates. The colors can easily be changed within the Rscript. Finally, the range for legend and colors can also be modified. 

If grey boxes appear in the plot, the bin is probably too high (around 50k is seems to fail). It will require modification, as I usually work with fungal genomes of about 15-100Mb, so it is optimized for this range. 
