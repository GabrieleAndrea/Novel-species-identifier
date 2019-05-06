# Novel-species-identifier

                                         Version 1.0, May 2019

Novel species identifier - METAnnotatorX custom-made script
Copyright (C) 2019 Lugli Gabriele Andrea

Novel species identifier is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Novel species identifier is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

1. What could Novel species identifier do?:

The script identifies those putative contigs that belongs to novel species of a
selected genus starting from a shotgun metagenomic sample analyzed through the
METAnnotatorX pipeline.

2. System requirements:

Novel species identifier should run on all Unix platforms, although it has not 
tested in all platforms.

3. Installation:

The script does not need any installation. To be executed, the script needs a 
functional copy of METAnnotatorX pipeline.

4. Software requirements and dependencies:

Blastn programs (specifically blastn and makeblastdb programs)

5. Databases:

The user has to build three custom databases based on the interested bacterial genus
sequences.

	A. Collection of all protein sequences:
bin/./prerapsearch –f T –d /folder/sequences.fasta –n /folder/protein

	B. Collection of all nucleotide sequences:
makeblastdb -in /folder/sequences.fasta -dbtype nucl -out /folder/nuc_genomes

	C. Collection of all nucleotide genes' sequences:
makeblastdb -in /folder/sequences.fasta -dbtype nucl -out /folder/nuc_genes

6. Usage:

Place the script novel_species_identifier.sh in the same folder of the METAnnotatorX
pipeline, e.g., alongside with folder "output", "parameters_folder", "software" and
the file "parameters". In other words, after a complete run of the pipeline. The
assembly phase and classification at genus level of each contig are mandatory.

At line 77 to 79 the user must enter the path of each databases.

Usage option:

-ic, --idcontigs	maximum identity for contigs selection"

-ig, --idgenes		maximum identity for genes selection"

-g, --genus		selected genus"

example: "./novel_species_identifier.sh -ic 94 -ig 90 -g Bifidobacterium"

7. Results:

In the folder novel_genus are located the results of each analysis step.
