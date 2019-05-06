#!/usr/bin/env bash

# Novel species identifier - METAnnotatorX custom-made script
# Copyright (C) 2019 Lugli Gabriele Andrea

# Novel species identifier is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Novel species identifier is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#############################################################################
#                                                                           #
#   METAnnotatorX results has to be in the same folder of this script.      #
#   Expecially output, software, parameters_folder and parameters file.     #
#   Assembly and classification at genus level of each contig are needed.   #
#                                                                           #
#############################################################################

display_help() {
	echo "---------------------------------------------------------------------------------"
	echo ""
	echo "	Usage: $0 [option..]" >&2
	echo
	echo "	-ic, --idcontigs	maximum identity for contigs selection"
	echo "	-ig, --idgenes		maximum identity for genes selection"
	echo "	-g, --genus		selected genus"
	echo 
	echo "	example: $0 -ic 94 -ig 90 -g Bifidobacterium"
	echo 
	echo "---------------------------------------------------------------------------------"
	exit 1
}

case "$1" in 
	-h | --help )
		display_help
		exit 0
		;;
esac 

case "$1" in
	-ic | --idcontigs )
		identcontigs=($2)
		;;
esac 

case "$3" in
	-ig | --idgenes)
		identgenes=($4)
		;;
esac 

case "$5" in
	-g | --genus)
		selgenus=($6)
		;;
esac


if [[ $# -eq 0 ]] ; then
    display_help
    exit 0
fi

#variables
PROJECT=$(cat parameters_folder/project_name) ;
THREADNUM=$(cat parameters_folder/threads) ;
BACLENGTH=$(cat parameters_folder/bacteria_contigs_minimum_length) ;
DBPROT=$(echo '/folder/protein') ; #Collection of all protein sequences:		bin/./prerapsearch –f T –d /folder/sequences.fasta –n /folder/protein
DBCONT=$(echo '/folder/nuc_genomes') ; #Collection of all nucleotide sequences:		makeblastdb -in /folder/sequences.fasta -dbtype nucl -out /folder/nuc_genomes
DBGENE=$(echo '/folder/nuc_genes') ; #Collection of all nucleotide genes' sequences:	makeblastdb -in /folder/sequences.fasta -dbtype nucl -out /folder/nuc_genes

#selection of contigs
mkdir temp ;
mkdir novel_${selgenus} ;
perl software/Fasta1line.pl output/${PROJECT}_assembly/contigs.fasta temp/${PROJECT}_contigs_1line.fasta ;
grep "$selgenus" output/assembled_contigs_taxonomy/${PROJECT}_genera.txt > temp/${PROJECT}_just_${selgenus}.txt ;
gawk -F"\t" '{print $2}' temp/${PROJECT}_just_${selgenus}.txt > temp/${PROJECT}_${selgenus}_genera ;
grep -A 1 -f temp/${PROJECT}_${selgenus}_genera temp/${PROJECT}_contigs_1line.fasta > novel_${selgenus}/METAnnotatorX_${selgenus}_contigs.fasta ;
sed -i 's/--//g' novel_${selgenus}/METAnnotatorX_${selgenus}_contigs.fasta ;

#ORFs prediction and rapsearch
software/./prodigal.linux -f gff -a temp/aaORFs.fasta -i novel_${selgenus}/METAnnotatorX_${selgenus}_contigs.fasta -o temp/ORFs.gff -p meta ;
software/./rapsearch -q temp/aaORFs.fasta -d $DBPROT -o temp/bacteria.rapsearch -s f -e 1e-100 -l $BACLENGTH -z $THREADNUM -b 1 -v 1 -a T -p T ;

#CLEANER PHASE 1
cp temp/bacteria.rapsearch.m8 temp/bacteria_clean.rapsearch.m8
sed -i "s|\[\[1\]\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\[3\]\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\[h\]\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(3aS\,4S\,7aS)-7a-methyl-1\,5-dioxo-octahydro-1H-inden-4-yl\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(5-\ phosphoribosylamino)methylideneamino\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(5-phosphoribosylamino)\ methylideneamino\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(5-phosphoribosylamino)methylideneamino\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(R)-acetoin\ forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[(S)-acetoin\ forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\^PR\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[10\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[13\.6\ kda\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[14\.8\ kda\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[18\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[2+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[20\ kda\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[21\ kda\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[2Fe-2S\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[2fe2s\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-hydroxylauroyl\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-hydroxy-myristory\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-hydroxymyristoyl\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-hydroxy-phenyl\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-methyl-2-oxobutanoate\ dehydrogenase\ (acetyl-transferring)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[46\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[4Fe-2S-2O\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[4Fe-4S\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[5[']-phosphate\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[6\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[9[']\,10[']\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[AA\ 29-492\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acceptor\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetolactate\ synthase\ pyruvate\ dehydrogenase\ (Cytochrome)\ glyoxylate\ carboligase\ phosphonopyruvate\ decarboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetylCoA\ carboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetyl-CoA\ carboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetyl-CoA-carboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetyl-transferring\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ACP\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl\ carrier\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-\ carrier-protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acylating\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier-(ACP)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier-protein\ (ACP)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier-protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ADP/GDP-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ADP-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ADP-ribose\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Alcaligin\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Alcaligin-like\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[alpha\ type\ 1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[alpha\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[alternative\ form\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[amino\ acid\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ammonia\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[apurinic\ or\ apyrimidinic\ site\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[arbutin-salicin-cellobiose\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Arg8\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[asymmetrical\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ATP\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[B\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[B12\ dependent\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[bacillibactin\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[bacterial\ glycogen\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[benzylsuccinate\ synthase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[BtrI\ acyl-carrier\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Burkholderia\ phage\ BcepGomr\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Burkholderia\ phage\ KL1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[C\.\ elegans\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[C1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[carboxylating\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[carrier\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[CCR4-NOT\ transcription\ complex\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cd00048\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[citrate\ (pro-3S)-lyase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cl00173\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cl03829\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[COG1132\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\ diacylglycerol\ kinase\ catalytic\ domain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\ methyltransferase\,\ helicase\ and\ papain-like\ proteinase\ domains\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Alpha\ crystallin\ A\ chain\,\ short\ form\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Aspartic\ protease\;\ Endonuclease\;\ Reverse\ transcriptase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Endonuclease\ PI-MtuHIP\ (Mtu\ DnaB\ intein)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\:\ gamma-glutamyltransferase\ heavy\ chain\;\ gamma-glutamyltransferase\ light\ chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Gamma-glutamyltranspeptidase\ large\ chain\;\ Gamma-glutamyltranspeptidase\ small\ chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\:\ n-acetyl-gamma-glutamyl-phosphate\ reductase\ (n-acetyl-glutamate\ semialdehyde\ dehydrogenase)\ (nagsa\ dehydrogenase)\;\ acetylglutamate\ kinase\ (nag\ kinase)\ (agk)\ (n-acetyl-l-glutamate\ 5-phosphotransferase)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Neuregulin-4\ (NRG-4)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\:\ phosphatidylserine\ decarboxylase\ 1\ beta\ chain\;\ phosphatidylserine\ decarboxylase\ 1\ alpha\ chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\:\ phosphatidylserine\ decarboxylase\ beta\ chain\;\ phosphatidylserine\ decarboxylase\ alpha\ chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ Protein\ B\*\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[contains\:\ s-adenosylmethionine\ decarboxylase\ alpha\ chain\;\ s\ adenosylmethionine\ decarboxylase\ beta\ chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ thymosin\ alpha-1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[copper\ containing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[copper-containing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Cronobacter\ phage\ vB\_CsaM\_GAP32\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Cu-Zn\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cyclizing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[CysO\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cytidine(C)-cytidine(C)-adenosine\ (A)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cytochrome\ c\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[cytochrome\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[D(-)-3-hydroxyalkanoate\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[D(-)-3-hydroxybutyrate\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[D/E\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[deacetylating\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[decarboxylating\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[dimethylamine--corrinoid\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[dinitrogen\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[dinitrogen\ reductase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[dinitrogenase\ reductase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[DNA\ replication\,\ recombination\ and\ repair\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ 9g\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ IME08\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ ime09\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ JSE\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ RB14\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ RB32\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ RB51\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ T4\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ vB\_EcoM\_ACG-C40\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Enterobacteria\ phage\ vB\_EcoM-FV3\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ER\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Erwinia\ phage\ vB\_EamM-Y2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ 2\ JES-2013\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ e11/2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ phiEB49\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ rv5\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ vB\_EcoM\_FFH2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Escherichia\ phage\ wV7\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[FAD\,\ quinone\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[FMN\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fe/Mn\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[FeFe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fe-Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ferredoxin\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fe-S\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fe-Zn\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[flavin-containing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[flavocytochrome\ C\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[formate-C-acetyltransferase\ 1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[formate-C-acetyltransferase\ 2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[formate-C-acetyltransferase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[fragment\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fructose-bisphosphate\ aldolase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[G1/S-specific\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GcvH\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GDP-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GFA\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[glutamate--ammonia-ligase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[glutamate-ammonia-ligase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[glutamine-hydrolyzing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GPI-anchored\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GSEE\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[GTP\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[hemophilus\ phage\ HP1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[heparan\ sulfate\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[histone\ H3\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[hydroxy(phenyl)methyl\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[hydroxylating\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Hydroxymethylglutaryl-CoA\ reductase(NADPH)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[I/L\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[II\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[III\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\ phosphoribosylaminoimidazolecarboxamide\ formyltransferase\ (ec\ 2\.1\.2\.3)\;\ IMP\ cyclohydrolase\ (ec\ 3\.5\.4\.10)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ 4-alpha\ glucanotransferase\ (ec\ 2\.4\.1\.25)\ (oligo-1\,4-1\,4-glucantransferase)\ amylo-alpha-1\,6-glucosidase\ (ec\ 3\.2\.1\.33)\ (amylo-1\,6-glucosidase\ (dextrin\ 6-alpha-d-glucosidase)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ 4-hydroxy-2-oxoglutarate\ aldolase\;\ 2-dehydro-3-deoxy-phosphogluconate\ aldolase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ 8-oxoguanine\ DNA\ glycosylase\;\ DNA-(apurinic\ or\ apyrimidinic\ site)\ lyase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ alkaline\ phosphodiesterase\ (ec\ 3\.1\.4\.1)\;\ nucleotide\ pyrophosphatase\ (ec\ 3\.6\.1\.9)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Biotin\ carboxylase\;\ Biotin\ carboxyl\ carrier\ protein\ (BCCP)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Biotin\ carboxylase\;\ Biotin\ carboxyl\ carrier\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ cystathionine\ beta-lyase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ dihydroneopterin\ aldolase\ (ec\ 4\.1\.2\.25)\;\ dihydro-6-hydroxymethylpterin-pyrophosphokinase\ (ec\ 2\.7\.6\.3)\;\ dihydropteroate\ synthetase\ (ec\ 2\.5\.1\.15)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ dihydroneopterinaldolase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ dimethylallyltranstransferase\ (ec\ 2\.5\.1\.1)\ geranyltranstransferase\ (ec\ 2\.5\.1\.10)\;\ farnesyltranstransferas\ (ec\ 2\.5\.1\.29)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ glutamine\ amidotransferase\;\ cyclase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Glutamine\ amidotransferase\;\ Indole-3-glycerol\ phosphate\ synthase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Histidinol-phosphatase\ (EC\ 3\.1\.3\.15)\;\ Imidazoleglycerol-phosphate\ dehydratase\ (EC\ 4\.2\.1\.19)\ (IGPD)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ l-serin\ dehydratase\ (ec\ 4\.3\.1\.17)\ (l-serine\ deaminase)\;\ l-threonin\ dehydratase\ (ec\ 4\.3\.1\.19)\ (l-threonine\ deaminase)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ l-serine\ dehydratase\ (ec\ 4\.3\.1\.17)\ (l-serine\ deaminase)\;\ l-threonine\ dehydratase\ (ec\ 4\.3\.1\.19)\ (l-threonine\ deaminase)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Methylphosphotriester-DNA--protein-cysteine\ S-methyltransferase\;\ Methylated-DNA--protein-cysteine\ methyltransferase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Phosphopantothenoylcysteine\ decarboxylase\ (CoaC)\;\ Phosphopantothenate--cysteine\ ligase\ (Phosphopantothenoylcysteine\ synthase)\ (CoaB)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ phosphoribosyl-amp\ cyclohydrolase\ (ec\ 3\.5\.4\.19)\;\ phosphoribosyl-atp\ pyrophosphohydrolase\ (ec\ 3\.6\.1\.31)\;\ histidinol\ dehydrogenase\ (ec\ 1\.1\.1\.23)\ (hdh)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ Phosphoribosyl-AMP\ cyclohydrolase\;\ Phosphoribosyl-ATP\ pyrophosphatase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ prephenate\ dehydratase\;\ arogenate\ dehydratase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ thiamine-phosphate\ pyrophosphorylase\ (ec\ 2\.5\.1\.3)\;\ hydroxyethylthiazole\ kinase\ (ec\ 2\.7\.1\.50)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Includes\:\ threonine-phosphate\ decarboxylase\;\ Cobyric\ acid\ synthase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ UDP-glucose\ 4-epimerase\ (galactowaldenase)\;\ aldose\ 1-epimerase\ (mutarotase)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[isoleucine\ degradation\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[isomerizing\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[lipoamide\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[lipopolysaccharide\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[LSU\ ribosomal\ protein\ L11P\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[menaquinone\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[methionine\ synthase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[missing\ catalytic\ residues\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Mn/Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Mn\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[myosin\ heavy-chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[myosin\ light-chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Myosin-light-chain\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(+)\,\ L-lysine-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(+)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(P)(+)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(P)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(P)+\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(P)+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD(P)H\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD+\,\ L-lysine-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD+\,L-lysine-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NAD+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADH\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP(+)\,\ L-glutamate-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP(+)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP+\,\ L-glutamate-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP+\,\ L-lysine\ forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADP+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADPH\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[negative\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ni\,Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ni/Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ni\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NiFe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ni-Fe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NiFeSe\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NMDA\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[N-oxide-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NU+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[nuclear\ export\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[nucleotide\ sugar/triose\ phosphate\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[OFD1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[oxidative\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[P(3HO)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[P\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Pectobacterium\ phage\ phiTE\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[peptidoglycan\ hydrolase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Peptidoglycan\ synthetase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[phd\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[PI(3)P\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[PIN+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[positive\ chemotaxis\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Precursor\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[pro-3S\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Protein\ ADP-ribosylarginine\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[protein-PII\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[PSI+\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[pt\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[putative\ sequencing\ error\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[pyochelin\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[pyrroloquinoline-quinone\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Pyruvate\ dehydrogenase\ (acetyl-transferring)\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[pyruvate\,\ phosphate\ dikinase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[quinone\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Rad50-Mre11-Xrs2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ribulose\ forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ribulose-bisphosphate\ carboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ribulose-bisphosphate-carboxylase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ribulose-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[RNA-polymerase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Salmonella\ phage\ Vi06\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[second\ part\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[see\ also\ xtr\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[serine\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[seryl-trnaser\ selenium\ transferase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Shigella\ phage\ Shfl2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Skp1-protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[SSC1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[starch\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Streptococcus\ thermophilus\ bacteriophage\ Sfi11\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[sugar\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Tau\ protein\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[TIGR03118\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[TIR\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[transaldolase\ and\ glucose-6-phosphate\ isomerase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[tri-or\ tetra-\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[truncated\ ORF\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ubiquinone\ cytochrome\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ubiquinone\,\ cytochrome\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[ubiquinone\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[UDP-forming\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[UI\:99417512\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[V\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[3-methyl-2-oxobutanoate\ dehydrogenase\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[citrate\ -lyase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[citrate\ \[pro-3s\]-lyase\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Pyruvate\ dehydrogenase\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Pyruvate\ dehydrogenase\ \]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[NADPH\,\ B-specific\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[nadph\,\ b-specific\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[VPS32\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Colombo\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Giza\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ban5\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ageratum\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Quivican\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Y4536\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Y322\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[IT\:Sic2/2\:04\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Lembang\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Kochi\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Togo\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Nigeria\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Madagascar\:Morondova\:2001\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Venezuela\:Zulia\:2004\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Uganda\:Kampala\:2008\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[China\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Hainan\ 8\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Jamaica\:St\.\ Elizabeth\:2004\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India/Ahmedabad/2014\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[YN323\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[IT\:Sic2/2\:04\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Pakistan\:Lahore1\:2004\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Panipat\:Papaya\:2008\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[CM\:Lys1sp2\:09\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Sonipat\ EL14A\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Munthal\ EL37\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Sonipat\:EL10\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:\ Faizabad\:\ Cow\ Pea\:2012\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[91\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Guatemala\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Ama\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Gurdaspur\:Okra\:2013\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Amadalavalasa\:Hibiscus\:2007\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Bahraich\:2007\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Japan\:Fukui\:2001\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Mexico\:Yucatan\:2004\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Tomato\:Aragua\:2003\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Fz1\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Hoa\ Binh\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India/Ahmedabad/2014\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Coimbator\:OYCO1\:2005\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Dharwad\ OYDWR2\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[India\:Aurangabad\:OY164\:2006\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[2004\:New\ Delhi\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Hn2\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[G52\]||g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[\_\_|\\[|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acetyl-CoA-carboxylaseligase\ \[Pyrenophora\ tritici-repentis\ Pt-1C-BFP\]|\[Pyrenophora\ tritici-repentis\ Pt-1C-BFP\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[acyl-carrier\ protein\ \[Phaeodactylum\ tricornutum\ CCAP\ 1055/1\]|\[Phaeodactylum\ tricornutum\ CCAP\ 1055/1\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[bovine\ \[Plasmodium\ yoelii\ yoelii\ 17XNL\]|\[Plasmodium\ yoelii\ yoelii\ 17XNL\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Contains\:\ HCF\ N-terminal\ chain\ \[Schistosoma\ mansoni\]|\[Schistosoma\ mansoni\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[deaminase\ reductase\ \[Burkholderia\ pseudomallei\]|\[Burkholderia\ pseudomallei\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[deaminase\ reductase\ \[Burkholderia\ pseudomallei\]|\[Burkholderia\ pseudomallei\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[deaminase\,\ reductase\ \[Burkholderia\ pseudomallei\ K96243\]|\[Burkholderia\ pseudomallei\ K96243\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[deaminase\,\ reductase\ \[Burkholderia\ pseudomallei\ K96243\]|\[Burkholderia\ pseudomallei\ K96243\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[EC\:\ \[Bacillus\]|\[Bacillus\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[EC\:\ \[Bacillus\]|\[Bacillus\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[gm18\ methyltransferase)\ (trna\ methylase\ 3)\ \[Candida\ dubliniensis\ CD36\]|\[Candida\ dubliniensis\ CD36\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[human\ \[Cryptococcus\ neoformans\ var\.\ neoformans\ JEC21\]|\[Cryptococcus\ neoformans\ var\.\ neoformans\ JEC21\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\ B\ component\ of\ phosphotransferase\ system\ and\ transcriptional\ \[Clostridiales\]|\[Clostridiales\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\ B\ component\ of\ phosphotransferase\ system\ and\ transcriptional\ \[Clostridiales\]|\[Clostridiales\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ 2-hydroxyhepta-2\,4-diene-1\,7-dioate\ isomerase\ \[Salmonella\ enterica\]|\[Salmonella\ enterica\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[includes\:\ 2-hydroxyhepta-2\,4-diene-1\,7-dioate\ isomerase\ \[Salmonella\ enterica\]|\[Salmonella\ enterica\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Oryza\ sativa\ (ISS)\ \[Ostreococcus\ tauri\]|\[Ostreococcus\ tauri\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Oryza\ sativa\ (ISS)\,\ partial\ \[Ostreococcus\ tauri\]|\[Ostreococcus\ tauri\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[Prec\ (IC)\ \[Ostreococcus\ tauri\]|\[Ostreococcus\ tauri\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[rpteom\ \[Elusimicrobium\ minutum\]|\[Elusimicrobium\ minutum\]|g" temp/bacteria_clean.rapsearch.m8 ;
sed -i "s|\[rpteom\ \[Elusimicrobium\ minutum\]|\[Elusimicrobium\ minutum\]|g" temp/bacteria_clean.rapsearch.m8 ;

#CLEANER PHASE 2
grep "MULTISPECIES" temp/bacteria_clean.rapsearch.m8 > temp/bacteria_multispecies.rapsearch.m8 ;
grep -v "MULTISPECIES" temp/bacteria_clean.rapsearch.m8 > temp/bacteria_no_multispecies.rapsearch.m8 ;
sed -i 's/\]/\ sp.\]/g' temp/bacteria_multispecies.rapsearch.m8 ;
cat temp/bacteria_multispecies.rapsearch.m8 temp/bacteria_no_multispecies.rapsearch.m8 > temp/bacteria_modded.rapsearch.m8 ;

grep "\[\[" temp/bacteria_modded.rapsearch.m8 > temp/bacteria_doublesquare.rapsearch.m8 ;
grep -v "\[\[" temp/bacteria_modded.rapsearch.m8 > temp/bacteria_no_doublesquare.rapsearch.m8 ;
sed -i 's/\[\[/\[/g' temp/bacteria_doublesquare.rapsearch.m8 ;
sed -i 's/\]//' temp/bacteria_doublesquare.rapsearch.m8 ;
cat temp/bacteria_doublesquare.rapsearch.m8 temp/bacteria_no_doublesquare.rapsearch.m8 > temp/bacteria_corrected.rapsearch.m8 ;
sed 's/\[/\[Bacteria\_/g' temp/bacteria_corrected.rapsearch.m8 > temp/bacteria_final.rapsearch.m8

cat temp/bacteria_final.rapsearch.m8 | sort -u -k1,1 > temp/temp.rapsearch.m8 ;
gawk -F'\t' '$4 > 75 {print ;}' temp/temp.rapsearch.m8 > temp/temp_rapsearch_2.m8 ;
gawk '!/MULTISPECIES/' temp/temp_rapsearch_2.m8 > novel_${selgenus}/1_aaORFs.rapsearch ;
sed -i -e 's|\/|-|g' -e 's|\:|-|g' novel_${selgenus}/1_aaORFs.rapsearch ;

#split records into files
grep ">" temp/aaORFs.fasta > temp/temp_output.txt &&
cat temp/temp_output.txt | cut -d "#" -f 1  > temp/temp_output2.txt &&
awk '{gsub(">" , ""); print}' temp/temp_output2.txt > temp/temp_output_3.txt &&
awk -F '_' '{out="orfs_for_node_"$2} out!=prev{close(prev)} {print >> out; prev=out}' temp/temp_output_3.txt &&

#IDENTIFICATION
grep ">" novel_${selgenus}/METAnnotatorX_${selgenus}_contigs.fasta > temp/contigs_name ;
cat temp/contigs_name | gawk '{print $1}' > temp/contigs_name_edit ;
sed -i 's/>//g' temp/contigs_name_edit ;

while read CNAME
do
	let count=$count+1
	grep "^${CNAME}_" novel_${selgenus}/1_aaORFs.rapsearch > temp/contig_${count}_NCBI
	sed 's/\t//g' temp/contig_${count}_NCBI > temp/contig_${count}_NCBI_edit
	sed -i 's/\[/\t/g' temp/contig_${count}_NCBI_edit

	cat temp/contig_${count}_NCBI_edit | gawk -F"\t" '{print $2}' > temp/contig_${count}_tail
	cut -d " " -f 2- temp/contig_${count}_tail > temp/contig_${count}_tail_species
	sed -i 's/\]/\t/g' temp/contig_${count}_tail_species
	cat temp/contig_${count}_tail_species | gawk -F"\t" '{print $1}' > temp/contig_${count}_tail_species_edit

	cat temp/contig_${count}_tail | gawk '{print $1}' > temp/contig_${count}_tail_genus
	sed -i 's/\]//g' temp/contig_${count}_tail_genus
	cat temp/contig_${count}_tail_genus | gawk '{print $1}' > temp/contig_${count}_tail_genus_edit

	paste temp/contig_${count}_tail_genus_edit temp/contig_${count}_tail_species_edit > temp/contig_${count}_species
	sed -i 's/\t/\ /' temp/contig_${count}_species

#	nomber of genes placed as number of orfs
	grep "^>${CNAME}_" temp/aaORFs.fasta > temp/contig_${count}_number_of_genes
	wc -l temp/contig_${count}_number_of_genes > temp/contig_${count}_rows
	cat temp/contig_${count}_rows | gawk '{print $1}' > temp/contig_${count}_number_rows
	export NUMORFS=$( cat temp/contig_${count}_number_rows )

	#IDENTIFICATION GENERA
	cat temp/contig_${count}_species | gawk '{print $1}' > temp/contig_${count}_genera
	sort temp/contig_${count}_genera | uniq -c > temp/contig_${count}_genera_unique
	sort -k1,1nr -k2,2 temp/contig_${count}_genera_unique > temp/contig_${count}_genera_unique_sorted
	cat temp/contig_${count}_genera_unique_sorted | gawk '{print $1}' > temp/contig_${count}_genera_numbers
	export HIGHERNUM1=$(sed -n '1p' temp/contig_${count}_genera_numbers)
	NUMRAPP1=$( echo "scale=3; $HIGHERNUM1 / $NUMORFS" | bc -l )
	if [[ "$HIGHERNUM1" == "$NUMORFS" ]] || [[ $(bc <<< "$NUMRAPP1 > 0.1999") -eq 1 && "$NUMORFS" > 3 ]]  ; then
		cat temp/contig_${count}_genera_unique_sorted | gawk '{print $2}' > temp/contig_${count}_genera_names
		sed -n '1p' temp/contig_${count}_genera_names > temp/contig_${count}_genera_name
		echo "$CNAME" > temp/num_${count}_genera
		paste temp/contig_${count}_genera_name temp/num_${count}_genera > temp/final_contig_${count}_name.genera
	fi

	#IDENTIFICATION SPECIES
	sort temp/contig_${count}_species | cut -d " " -f 1,2 > temp/contig_${count}_species_trimmed
	uniq -c temp/contig_${count}_species_trimmed > temp/contig_${count}_species_unique
	sort -k1,1nr -k2,2 temp/contig_${count}_species_unique > temp/contig_${count}_species_unique_sorted
	cat temp/contig_${count}_species_unique_sorted | gawk '{print $1}' > temp/contig_${count}_species_numbers
	export HIGHERNUM2=$( sed -n '1p' temp/contig_${count}_species_numbers )
	NUMRAPP2=$( echo "scale=3; $HIGHERNUM2 / $NUMORFS" | bc  )
	if [[ "$HIGHERNUM2" == "$NUMORFS" ]] || [[ $(bc <<< "$NUMRAPP2 > 0.4999") -eq 1 && "$NUMORFS" > 3 ]]  ; then
		sed -i 's/^\ //g' temp/contig_${count}_species_unique_sorted ;
		sed -i 's/\ /\t/' temp/contig_${count}_species_unique_sorted ;
		cat temp/contig_${count}_species_unique_sorted | gawk -F"\t" '{print $2}' > temp/contig_${count}_species_names
		sed -n '1p' temp/contig_${count}_species_names > temp/contig_${count}_species_name
		echo "$CNAME" > temp/num_${count}_species
		paste temp/contig_${count}_species_name temp/num_${count}_species > temp/final_contig_${count}_name.species
	fi

	rm -rf contig_${count}_*
done < temp/contigs_name_edit

#EXTRAPOLATION SPECIES AND GENERA
cat temp/*.species | gawk -F"\t" '{print $2,"\t",$1}' > temp/${PROJECT}_species.txt ;
cat temp/*.genera | gawk -F"\t" '{print $2,"\t",$1}' > temp/${PROJECT}_genera.txt ;

rm -rf orfs_for_node_*

#creation of report file
awk '
BEGIN { OFS="\t" }
(NR==FNR) || ($1 in vals) {
    vals[$1][ARGIND] = $2
}
END {
    printf "%s%s", "Query", OFS
    for (fileNr=1; fileNr<=ARGIND; fileNr++) {
        printf "%s%s", ARGV[fileNr], (fileNr<ARGIND ? OFS : ORS)
    }
    for (key in vals) {
        printf "%s%s", key, OFS
        for (fileNr=1; fileNr<=ARGIND; fileNr++) {
            val = (fileNr in vals[key] ? vals[key][fileNr] : "NO MATCH")
            printf "%s%s", val, (fileNr<ARGIND ? OFS : ORS)
        }
    }
}' < temp/${PROJECT}_just_${selgenus}.txt temp/${PROJECT}_{genera.txt,species.txt} > novel_${selgenus}/${PROJECT}_1st_step_report.txt

#selection of unidentified contigs
sed 's/\ /\t/g' temp/${PROJECT}_species.txt > temp/${PROJECT}_species_tab.txt ;
gawk -F"\t" '{print $1}' temp/${PROJECT}_species_tab.txt > temp/${PROJECT}_just_species.txt ;
grep -v -f temp/${PROJECT}_just_species.txt temp/${PROJECT}_${selgenus}_genera > temp/${PROJECT}_interesting_contigs ;
grep -A 1 -f temp/${PROJECT}_interesting_contigs temp/${PROJECT}_contigs_1line.fasta > novel_${selgenus}/${PROJECT}_1st_step_unidentified_contigs.fasta ;
sed -i -e 's/--//g' novel_${selgenus}/${PROJECT}_1st_step_unidentified_contigs.fasta ;

#SELECTION based of contigs and genes identity values
blastn -query novel_${selgenus}/${PROJECT}_1st_step_unidentified_contigs.fasta -db $DBCONT -evalue 1e-20 -gapopen 1 -gapextend 1 -max_target_seqs 1 -num_threads $THREADNUM -outfmt "6 qseqid pident qcovs evalue stitle" -out temp/${PROJECT}.blastn ;

gawk '!seen[$1]++' temp/${PROJECT}.blastn > novel_${selgenus}/2_contigs.blastn ;
gawk -F"\t" -v iden="$identcontigs" '$2 < iden {print ;}' novel_${selgenus}/2_contigs.blastn > temp/${PROJECT}_interesting_results ;
gawk -F"\t" -v iden="$identcontigs" '$2 > iden {print ;}' novel_${selgenus}/2_contigs.blastn > temp/${PROJECT}_not_interesting_results ;
gawk -F"\t" '{print $1}' temp/${PROJECT}_not_interesting_results > temp/${PROJECT}_not_interesting_contigs ;
grep -v -f temp/${PROJECT}_not_interesting_contigs temp/${PROJECT}_interesting_contigs > temp/${PROJECT}_2nd_step_unidentified_contigs.fasta ;
grep -A 1 -f temp/${PROJECT}_2nd_step_unidentified_contigs.fasta temp/${PROJECT}_contigs_1line.fasta > novel_${selgenus}/${PROJECT}_2nd_step_unidentified_contigs.fasta ;
#gawk -F"\t" '{print $1}' temp/${PROJECT}_interesting_results > temp/${PROJECT}_interesting_contigs ;
#grep -A 1 -f temp/${PROJECT}_interesting_contigs temp/${PROJECT}_contigs_1line.fasta > novel_${selgenus}/${PROJECT}_2nd_step_unidentified_contigs.fasta ;
sed -i 's/--//g' novel_${selgenus}/${PROJECT}_2nd_step_unidentified_contigs.fasta ;

software/prodigal.linux -f gff -i novel_${selgenus}/${PROJECT}_2nd_step_unidentified_contigs.fasta -d temp/${PROJECT}_ORFs.fasta -o temp/${PROJECT}_temp.gff -p meta  ;
blastn -query temp/${PROJECT}_ORFs.fasta -db $DBGENE -out novel_${selgenus}/3_genes.blastn -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid pident qcovs stitle" -num_threads $THREADNUM ;
gawk '!seen[$1]++' novel_${selgenus}/3_genes.blastn > temp/3_gene_unique_blast ;
gawk -F"\t" -v ident="$identgenes" '$2 < ident {print ;}' temp/3_gene_unique_blast > temp/CDS_genomic_int_unique.txt ;
cat temp/CDS_genomic_int_unique.txt | gawk -F"\t" '{print $1}' > temp/CDS_genomic_list ;
sed -i 's/\_/\t/g' temp/CDS_genomic_list ;
awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' temp/CDS_genomic_list > temp/CDS_genomic_list_contigs ;
awk '!seen[$1]++' temp/CDS_genomic_list_contigs > temp/CDS_genomic_list_contigs_unique

	cat temp/3_gene_unique_blast | gawk -F"\t" '{print $1}' > temp/CDS_total_list ;
	sed -i 's/\_/\t/g' temp/CDS_total_list ;
	awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' temp/CDS_total_list > temp/CDS_total_list_contigs ;
	awk '!seen[$1]++' temp/CDS_total_list_contigs > temp/CDS_total_list_contigs_unique ;
	grep -v -f temp/CDS_total_list_contigs_unique temp/${PROJECT}_2nd_step_unidentified_contigs.fasta > temp/CDS_total_list_contigs_lost ;
	cat temp/CDS_genomic_list_contigs_unique temp/CDS_total_list_contigs_lost > temp/CDS_genomic_list_contigs_unique_and_lost

grep -A 1 -f temp/CDS_genomic_list_contigs_unique_and_lost temp/${PROJECT}_contigs_1line.fasta > novel_${selgenus}/${PROJECT}_3rd_step_unidentified_contigs.fasta ;
sed -i 's/--//g' novel_${selgenus}/${PROJECT}_3rd_step_unidentified_contigs.fasta ;

#rm -rf temp/* ;

echo ""
echo "  Attached results are reported into the novel_${selgenus} directory "
echo ""

exit 1  
