gc_plotter
=========

## Simple tool to map GC bins on chromosome ideograms with existing tools

Inspired by (and code directly taken from) [CMplot](https://github.com/YinLiLin/R-CMplot)
<br>
The Unix commands for creating GC content tracks also have been taken from [here](https://wiki.bits.vib.be/index.php/Create_a_GC_content_track)

GC_PLOTTER came from the idea that ascomycetes and basidiomycetes fungal genomes contain unusual changes called Repeat-Induced Point mutations (RIP). Tools like [OcculterCut](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4943192/#evw121-B11) can identify RIP areas, but visualization is lacking. In addition, GC_PLOTTER does essentially the same thing as the vanilla OcculterCut run, with a finer level of detail (genome regions by any window size of the users choosing, with the location on the genome), in 5 lines of bash. 

Plotting is done as part of the R package (OcculterCut needs a separate GNUplot call). As long as R, bedtools, and samtools are installed it should work on any unix-like system, as the only other dependencies are basic unix tools.

---

### Install Dependencies

#### samtools

Most likely one will already have [samtools](http://www.htslib.org/) installed if you are doing almost any bioinformatics. If not, please install from [here](http://www.htslib.org/download/). 

#### R programming language

R is a statistical programming language. Please download and install from [here](http://archive.linux.duke.edu/cran/) .


### bedtools

Finally, we have [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
. The Bedtools package makes the windows and counts nucleotides occurrences in the fasta file. 

---

### Run

#### You must have bedtools, samtools, R, and basic unix utils installed

Simply change the script to executable

```bash
chmod 755 gc.sh
```

#### Linux
Change the script from #!/bin/bash to /usr/bin/bash or appropriate shell. 

Then run the script with 3 arguments for the "linear" form, first is the reference genome fasta, the second is the desired width to scan for GC percenatge, "d" is for density plot (linear plot), and "jpg" for the file type. The filetype can be "jpg", "pdf", "tiff", or "png" (default).

```bash
./gc /home/reference.fasta 1000 d jpg
```
The following files will be created in the directory which has the reference genome:

".sizes" &nbsp;&nbsp;&nbsp; file which specifies all the chromosomes/scaffolds and the lengths in basepairs <br>
".bed"   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  bed file <br>
".fai"   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  fasta index from samtools <br>
".igv"   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; igv viewer file with positions, bins, and gc content. Output file needed for gc_plotter<br>

---

### Plots

Two simple examples are the plasmodium species *P. falciparum* and *P. knowlesi*

*P. falsiparum* is well-known as a very [AT rich](https://genomevolution.org/wiki/index.php/Plasmodia_comparative_genomics) and GC poor genome. While *P. knowlesi* is still AT rich, but much less rich than *falciparum*. We can plot these and see if we get the expected results.

The genomes were downloaded from [ncbi](https://www.ncbi.nlm.nih.gov/genome/?term=txid1245013[Organism:noexp]) and used as is, except for a renaming. 

```bash
./gc.sh /home/plasmodium_knowlesi_genomic.fna 10000 d jpg
```
<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/figures/plas1.jpg">
<img src="figures/plas1.jpg" height="800px" width="800px">
</a>
</p>

---
Same for the second genome:

```bash
./gc.sh /home/plasodium_falciparum_genomic.fna 10000 d jpg
```
<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/figures/plas.jpg">
<img src="figures/plas.jpg" height="800px" width="800px">
</a>
</p>

As we can see, there is lots of red for *falciparum* indicating that yes, it seems to be well within the 20's (23%), and *knowlesi* is around 40% as expected. 

---

We can also try a slightly larger genome to see if it scales:

```bash
./gc.sh /home/arabidopsis_thaliana_genomic.fna 10000 d jpg
```

<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/figures/arabidopsis.jpg">
<img src="figures/arabidopsis.jpg" height="800px" width="800px">
</a>
</p>

---

### Circular Plot

I have now added circular plots. It has an extended command on the terminal taking a hard-masked fasta. If you have a soft-masked fasta ("atcg" instead of "ATCG") then please run the "mask.py" script to coerce to hardmasked. I suppose it might be possible to use either, this feature will be added later.

```python
python mask.py soft-masked.fasta
```

The full command line for the circular plot:
```bash
./gc.sh /home/knufia.genome.fasta 10000 c jpg /home/Desktop/knufia.hard_masked.fasta
```
This will then generate something like this:

<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/figures/Knufia_polished.jpg">
<img src="figures/Knufia_polished.jpg" height="600px" width="800px">
</a>
</p>

The figure above is purposefully enlarged to show the small contigs. Normally, the spacing will not be so explicit on the legend. 

---

Finally, there is a histogram function to see what the bins are graphically. If you would like to use the "histoplot" please install ggplot2 from R. (install.packages(""ggplot2")). 

```bash
/path/to/gc.sh exophiala.genome.fasta 1000 h green black
```
<p align="center">
<a href="https://raw.githubusercontent.com/TheRincon/gc_plotter/figures/Exo_1000.png">
<img src="figures/Exo_1000.png" height="600px" width="800px">
</a>
</p>

This also includes a table with the data for use in other contexts.

---

### Hints

The legend can be freely changed in the linear plot, simply replace "topright" or "bottomright" on the last line in gcplot.r to another location or coordinates. The colors can easily be changed within the Rscript under col (or col2 in case of the circular plot). 

If grey boxes appear in the plot, the bin is probably too high (around 50k or over is seems to fail). Set the bin or "window" size lower. Usually 1 order of magnitude lower will work, i.e. 100000 => 10000. It will require modification, as I usually work with fungal genomes of about 15-100Mb, so it is optimized for this range. 

The title must be changed by hand (for now) in the plot Rscript. It will be under paste("..."), just search for this string in the file. 

---

### To Do

1. Make the legend more "publication worthy", as it looks very simple now. 
2. Add "title" option. User should be able to pass an argument. 
3. Gene density plots as another plot option/function?
4. Add a scale for size in b/n the plots on circular.
5. mkdir tmp_$DATE and rm after use, unless the user wants to keep the files with "k" option.
6. Proper options flags for bash?
7. Hard and softmasking plot options. 
