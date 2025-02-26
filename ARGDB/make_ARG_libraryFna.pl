#!/usr/bin/perl
use strict;
use Getopt::Long;
use File::Basename;

my $PROG = basename $0;

sub display_version {
  print STDERR " 
  ###############= VERSION =################
    2022/07/13
    $PROG v0.0.2 (dev), 
    hudeneil (hudeneil\@asiapathgenomics.com)
   ";
  exit 0;
}


sub version_log {
  print STDERR " 
  ###############= VERSION =################
    2022/07/13
    $PROG v0.0.2 (dev), 
    hudeneil (hudeneil\@asiapathgenomics.com)

Version log:
2022-03-04 v0.0.1 (dev) create
2022-07-13 v0.0.2 (dev) 

   ";
  exit 0;
}




sub usage {
    my $error_message = $_[1];
    my $exit_code = @_ ? shift : 64;
    print STDERR <<EOF;

This is summary all infomation from bwa, samtools coverage & idxstats

Usage: $PROG [options] 
Options:
Input/Output:
 column data input: 
  
  --inputFna,-i  STR      ARG homolog fasta
  --ARO_index  STR        ARO index file (download from CARD) [aro_index.tsv]
  
  --card_prevalence       ARG prevalence data [card_prevalence.txt]
  
  --output,-o  STR        output result path [result]

  --help                  Print this message
  --version               Print version information
  --version_log           Print version log

example:


perl make_ARG_libraryFna.pl \
--inputFna /home/shareFiles/APGDB/ARGDB/unzip/CARD_homolog_Data/all_nucleotide_homolog_fasta.fasta \
--ARO_index /home/shareFiles/APGDB/ARGDB/unzip/CARD_homolog_Data/aro_index.tsv \
--card_prevalence /home/shareFiles/APGDB/ARGDB/unzip/prevalence/card_prevalence.txt \
-o /home/shareFiles/APGDB/ARGDB/2022Sep02/library


$error_message
EOF
    exit $exit_code;
}

## Option Input Required variables
my ($inputFna, $output);
my $ARO_index;
my $card_prevalence;

GetOptions(
    "h" => \&display_help,      ## -h
    "v" => \&display_version,   ## -v
    "o=s" => \$output,          ## -o
    "i=s" => \$inputFna,        ## -i

    "help" => \&display_help,                   ## --help
    "version" => \&display_version,             ## --version
    "version_log" => \&version_log,             ## --version_log
    "inputFna=s" => \$inputFna,                 ## --inputcov
    "ARO_index=s" => \$ARO_index,               ## --ARO_index
    "card_prevalence=s" => \$card_prevalence,   ## --ARO_index
    "output=s" => \$output,                     ## --output
);

my $errorTime = 0;

if ( ! defined $inputFna) { print STDERR "ERROR: No --inputFna input!!!!!\n"; $errorTime++;}
if ( ! defined $ARO_index) { print STDERR "ERROR: No --ARO_index input!!!!!\n"; $errorTime++;}
if ( ! defined $card_prevalence) { print STDERR "ERROR: No --card_prevalence input!!!!!\n"; $errorTime++;}
if ( ! defined $output) { $output = "result";}

if ( ! -f $inputFna) { print STDERR "ERROR: Input No reads hit table not exists!!!!!\n"; $errorTime++;}
if ( ! -f $ARO_index) { print STDERR "ERROR: --ARO_index file not exists!!!!!\n"; $errorTime++;}
if ( ! -f $card_prevalence) { print STDERR "ERROR: --card_prevalence file not exists!!!!!\n"; $errorTime++;}

sub display_help { usage(0);}


my ($message, $errors);
if ($errorTime > 0)
{
    if ($errorTime > 1){$errors = "errors"} else{$errors = "error"}
    $message = "
There are $errorTime $errors at your input
Please check your input carefully";
    usage(1, $message);
}

`
if [ ! -d "$output" ] ;then
    mkdir -p $output
fi
`;





## Global tmp variables
my $header;
my @tmp;
my ($accessionID, $ARO_ids, $name);


## 確認是可以刪除的ARO id
my @blackList=qw/3003820 3004056 3005069/;
my %blackList;
foreach $ARO_ids (@blackList)
{
    $blackList{$ARO_ids} = 1;
}

#### aro_index.tsv content
#### download from https://card.mcmaster.ca/latest/data
# ARO Accession	CVTERM ID	Model Sequence ID	Model ID	Model Name	ARO Name	Protein Accession	DNA Accession	AMR Gene Family	Drug Class	Resistance Mechanism	CARD Short Name
# ARO:3005099	43314	6143	3831	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)	23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)	AAB60941.1	AF002716.1	Erm 23S ribosomal RNA methyltransferase	lincosamide antibiotic;macrolide antibiotic;streptogramin antibiotic	antibiotic target alteration	Spyo_ErmA_MLSb
# ARO:3000317	36456	51	874	AAC(1)	AAC(1)-I	ADH03009.1	HM036080.1	AAC(1)	aminoglycoside antibiotic	antibiotic inactivation	AAC(1)-I
# ARO:3002523	38923	8144	1781	AAC(2')-Ia	AAC(2')-Ia	AAA03550.1	L06156.2	AAC(2')	aminoglycoside antibiotic	antibiotic inactivation	AAC(2')-Ia
# ARO:3002524	38924	85	746	AAC(2')-Ib	AAC(2')-Ib	AAC44793.1	U41471.1	AAC(2')	aminoglycoside antibiotic	antibiotic inactivation	AAC(2')-Ib
# ARO:3002525	38925	4719	1246	AAC(2')-Ic	AAC(2')-Ic	CCP42991.1	AL123456.3	AAC(2')	aminoglycoside antibiotic	antibiotic inactivation	AAC(2')-Ic

## Columns:
#0 ARO Accession	
#1 CVTERM ID	
#2 Model Sequence ID	
#3 Model ID	
#4 Model Name	
#5 ARO Name	
#6 Protein Accession	
#7 DNA Accession	
#8 AMR Gene Family	
#9 Drug Class	
#10 Resistance Mechanism	
#11 CARD Short Name

my $ARO_Accession;
my $AMR_Gene_Family;
my $Resistance_Mechanism;


my %AMR_Gene_Family;
my %Resistance_Mechanism;
my %AMRGFid;
my %id2Name;
my %GFid2ARO_ids;

my $AMRGFid;
my $count = 0;


open(IN, "$ARO_index")||die "$!";

$header =<IN>;
while(<IN>)
{
    chomp;
    @tmp = split "\t", $_;
    $ARO_Accession = $tmp[0];
    $AMR_Gene_Family = $tmp[8];
    $Resistance_Mechanism = $tmp[10];
    
    if ($AMR_Gene_Family=~/AAC\(|ANT\(|APH\(/) {$AMR_Gene_Family = "AAC, ANT, APH gene family"}
    if ($AMR_Gene_Family=~/glycopeptide resistance gene cluster|[Vv]an/) {$AMR_Gene_Family = "glycopeptide resistance gene cluster (van)"}
    
    
    $AMR_Gene_Family{$ARO_Accession} = $AMR_Gene_Family;
    $Resistance_Mechanism{$ARO_Accession} = $Resistance_Mechanism;
    
    ## 命名 gene family 流水號 
    if ( ! exists ($AMRGFid{$AMR_Gene_Family}))
    {
        $count++;
        $AMRGFid = "GF_$count";
        $AMRGFid{$AMR_Gene_Family} = $AMRGFid;
        $id2Name{$AMRGFid} = $AMR_Gene_Family;
    }
    
    $AMRGFid = $AMRGFid{$AMR_Gene_Family};
    $GFid2ARO_ids{$AMRGFid}="$GFid2ARO_ids{$AMRGFid},$ARO_Accession";
    #
}
close IN;





### card_prevalence.txt
### download from https://card.mcmaster.ca/latest/variants
# ARO Accession	Name	Model ID	Model Type	Pathogen	NCBI Plasmid	NCBI WGS	NCBI Chromosome	NCBI Genomic Island	Criteria	ARO Categories
# ARO:3002501	PDC-4	1	protein homolog model	Pseudomonas aeruginosa	0	0.02	0	0	perfect_strict	antibiotic inactivation; carbapenem; cephalosporin; monobactam
# ARO:3002999	CblA-1	2	protein homolog model	Phocaeicola dorei	0	1.14	0	0	perfect	antibiotic inactivation; cephalosporin
# ARO:3002999	CblA-1	2	protein homolog model	Phocaeicola dorei	0	1.14	0	0	perfect_strict	antibiotic inactivation; cephalosporin
# ARO:3003390	Escherichia coli ompF with mutation conferring resistance to beta-lactam antibiotics	3	protein variant model	Escherichia coli	0	0.01	0	0	perfect_strict	carbapenem; cephalosporin; cephamycin; monobactam; penam; penem; reduced permeability to antibiotic
# ARO:3001109	SHV-52	4	protein homolog model	Klebsiella pneumoniae	0	0.11	0	0	perfect	antibiotic inactivation; carbapenem; cephalosporin; penam
# ARO:3001109	SHV-52	4	protein homolog model	Klebsiella pneumoniae	0	0.16	0	0	perfect_strict	antibiotic inactivation; carbapenem; cephalosporin; penam
# ARO:3002867	dfrF	5	protein homolog model	Erysipelatoclostridium ramosum	0	19.23	0	0	perfect_strict	antibiotic target replacement; diaminopyrimidine antibiotic
# ARO:3002867	dfrF	5	protein homolog model	Streptococcus dysgalactiae	0	0	2.56	0	perfect_strict	antibiotic target replacement; diaminopyrimidine antibiotic

#0 ARO Accession	
#1 Name	
#2 Model ID	
#3 Model Type	
#4 Pathogen	
#5 NCBI Plasmid	
#6 NCBI WGS	
#7 NCBI Chromosome	
#8 NCBI Genomic Island	Criteria	
#9 ARO Categories

my $Pathogen;
my %Gene_prevalence;

open(IN, "$card_prevalence")||die "$!";

$header =<IN>;
while(<IN>)
{
    chomp;
    @tmp =split "\t", $_;
    $ARO_Accession = $tmp[0];
    $Pathogen = $tmp[4];
    
    if (! exists($Gene_prevalence{$ARO_Accession})){$Gene_prevalence{$ARO_Accession} = $Pathogen}
    else { $Gene_prevalence{$ARO_Accession} = "$Gene_prevalence{$ARO_Accession},$Pathogen" }
    #
}
close IN;






open(IN, "$inputFna")||die "$!";
# open(OUTheader, ">$output/header.txt")||die "$!";
# open(OUTfna, ">$output/library.fna")||die "$!";



my $family_groups;
my %family_groups;
my %No_family_groups;
my %all_ARO_Accession;
my %ARO2fasta;
my %ARO2header;
my %ARO2Name;


while(<IN>)
{
    chomp;
    
# >gb|HQ845196.1|+|0-861|ARO:3001109|SHV-52 [Klebsiella pneumoniae] 
# >gb|AF028812.1|+|392-887|ARO:3002867|dfrF [Enterococcus faecalis] 
# >gb|JX017365.1|+|244-1120|ARO:3001989|CTX-M-130 [Escherichia coli] 
# >gb|JN967644.1|+|0-813|ARO:3002356|NDM-6 [Escherichia coli]
    
    if(/>/)
    {
        if(/gb\|(.+)\|\+/){$accessionID = $1;}
        if(/ARO:(\d+)\|(.+) \[(.+)\]/)
        {
            $ARO_ids = "$1"; 
            $name = $2;
            $Pathogen = $3;
            $header =">ARO:$ARO_ids,$accessionID $name";
            # print OUTfna "$header\n";
            # print OUTheader "$header\n";
            # print OUTtaxidname "$ARO_ids\t$name\n";
            
            $ARO_Accession = "ARO:$ARO_ids";
            
            if ($Gene_prevalence{$ARO_Accession}=~/$Pathogen/){}
            else
            {
                if (! exists($Gene_prevalence{$ARO_Accession})){$Gene_prevalence{$ARO_Accession} = $Pathogen}
                else { $Gene_prevalence{$ARO_Accession} = "$Gene_prevalence{$ARO_Accession},$Pathogen" }
            }
            
            $ARO2header{$ARO_ids} = "$header";
            $ARO2Name{$ARO_ids} = "$name";
            $all_ARO_Accession{$ARO_ids} = 1;
            
            
            ## 確認是否為特定gene family
            if ( exists($AMR_Gene_Family{$ARO_Accession}) )
            {
                $AMR_Gene_Family = $AMR_Gene_Family{$ARO_Accession};
                $AMRGFid = $AMRGFid{$AMR_Gene_Family};
                
                if ( ! exists($family_groups{$AMRGFid}))
                {
                    $family_groups{$AMRGFid} = "$ARO_ids";
                }
                else
                {
                    $family_groups{$AMRGFid} = "$family_groups{$AMRGFid},$ARO_ids";
                }
            }
            else
            {
                $No_family_groups{$ARO_Accession} = 1;
            }
        }
    }
    else
    {
        $ARO2fasta{$ARO_ids} = $_;
        # print OUTfna "$_\n";
    }
}

close IN ;
# close OUTheader ;


open(OUTtaxidname, ">$output/library.taxid-name.txt")||die "$!";
open(OUTtaxid, ">$output/library.taxid.txt")||die "$!";
open(OUTfna, ">$output/library.fna")||die "$!";


print OUTtaxidname "0000000\tCARD-homolog_model\n";

foreach $ARO_ids(sort keys %ARO2fasta)
{
    $ARO_Accession = "ARO:$ARO_ids";
    
    my $prevalence_frequence = (split ",", $Gene_prevalence{$ARO_Accession} ) ;
    
    if (! exists($blackList{$ARO_ids}))
    {
        print OUTfna "$ARO2header{$ARO_ids}\n";
        print OUTfna "$ARO2fasta{$ARO_ids}\n";
        print OUTtaxidname "$ARO_ids\t$ARO2Name{$ARO_ids}\t$prevalence_frequence\t$Gene_prevalence{$ARO_Accession}\n";
    }
}
close OUTfna ;

print "Fasta output : $output/library.fna\n";

#### 重新定義無 gene family group 分類
my $count=0;
$AMRGFid = "0000000";
foreach $ARO_ids (sort keys %No_family_groups)
{
    $ARO_Accession = "ARO:$ARO_ids";
    if ($count==0)
    {
        $family_groups{$AMRGFid} = "$ARO_ids";
    }
    else 
    {
        $family_groups{$AMRGFid} = "$family_groups{$AMRGFid},$ARO_ids";
    }
    $count++;
}



foreach $AMRGFid (sort keys %family_groups)
{
    print OUTtaxidname "$AMRGFid\t$id2Name{$AMRGFid}\n";
    print OUTtaxid "$AMRGFid\t$family_groups{$AMRGFid}\n";
}

close OUTtaxidname;
close OUTtaxid;

print "library output : $output/library.taxid.fna\n";
print "library output : $output/library.taxid-name.fna\n";

