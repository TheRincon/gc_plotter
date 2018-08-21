import sys

def makeN(fasta):
	with open(fasta, "rt") as f:
		with open("hard_masked.fasta", "wt") as fout:
			for line in f:
				if line.startswith(">"):
					fout.write(line)
					print "hey"
				else:
					for i in line:
						if (i == "a" or i == "t" or i == "c" or i == "g"):
							i = "N"
							fout.write(i)
						else:
							fout.write(i)

if __name__ == "__main__":
	fasta = sys.argv[1]
	makeN(fasta)