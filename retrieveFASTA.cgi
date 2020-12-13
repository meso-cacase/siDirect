#!/usr/bin/perl

# retrieveFASTA.cgi?accession=XXX でFASTAをtxt形式で返す。
#
# 2009-10-18.y-naito
# 2010-11-10.y-naito（NCBIのURL変更）

use warnings ;
use strict ;
use WWW::Mechanize ;

my %query = get_query_parameters() ;  # CGIが受け取ったデータの処理
$query{'accession'} and (my $accession = $query{'accession'}) =~ s/^\s*(.*?)\s*$/$1/ ;  # 前後の空白文字を除去

unless ($accession and $accession =~ /^[\w\.]+$/){
	my $sampleseq =
'>sample sequence
ggctgccaag aacctgcagg aggcagaaga atggtacaaa tccaagtttg ctgacctctc
tgaggctgcc aaccggaaca atgacgccct gcgccaggca aagcaggagt ccactgagta
ccggagacag gtgcagtccc tcacctgtga agtggatgcc cttaaaggaa ccaatgagtc
cctggaacgc cagatgcgtg aaatggaaga gaactttgcc gttgaagctg ctaactacca
agacactatt ggccgcctgc aggatgagat tcagaatatg aaggaggaaa tggctcgtca
ccttcgtgaa taccaagacc tgctcaatgt taagatggcc cttgacattg agattgccac
ctacaggaag ctgctggaag gcgaggagag caggatttct ctgcctcttc caaacttttc
ctccctgaac ctgagggaaa ctaatctgga ttcactccct ctggttgata cccactcaaa
aaggacactt ctgattaaga cggttgaaac tagagatgga caggttatca acgaaacttc
tcagcatcac gatgaccttg aataaaaatt gcacacactc agtgcagcaa tatattacca' ;
	print_fasta_text($sampleseq) ;
	exit ;
}

my $mech = WWW::Mechanize->new ;
#$mech->get("https://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?dopt=fasta&sendto=t&val=$query{'accession'}") ;
$mech->get("https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=$query{'accession'}&dopt=fasta&sendto=t") ;
my $fasta = $mech->content ;
$fasta =~ s/\s+\z// ;
print_fasta_text($fasta) ;

exit ;

# ====================
sub get_query_parameters {  # CGIが受け取ったデータの処理
my $buffer = '' ;
if (defined $ENV{'REQUEST_METHOD'} and $ENV{'REQUEST_METHOD'} eq 'POST' and defined $ENV{'CONTENT_LENGTH'}){
	read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'}) ;
} elsif (defined $ENV{'QUERY_STRING'}){
	$buffer = $ENV{'QUERY_STRING'} ;
}
my %query ;
my @query = split /&/, $buffer ;
foreach (@query){
	my ($name,$value) = split /=/ ;
	if (defined $name and defined $value){
		$value =~ tr/+/ / ;
		$value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack('C', hex($1))/eg ;
		$name =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack('C', hex($1))/eg ;
		$query{$name} = $value ;
	}
}
return %query ;
} ;
# ====================
sub print_fasta_text {
my $fasta = $_[0] || 'not found.' ;
print "Content-type: text/plain; charset=utf-8\n\n$fasta\n" ;
exit ;
} ;
# ====================
