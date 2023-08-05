#!/usr/bin/perl

# --------------------------------
# siDirect web server, DETAIL page
# --------------------------------
# This page displays a detailed list of off-target sites.
#
# GGGenome https://GGGenome.dbcls.jp/ を利用してオフターゲットリストを表示する方法に変更（2023-07-01）
#
# 2009-04-28 Yuki Naito (@meso_cacase)
# 2023-07-10 Yuki Naito (@meso_cacase) Version 2.1

use warnings ;
use strict ;
use LWP::Simple ;

my $siDirect_top_url = 'http://siDirect2.RNAi.jp/' ;

my $timestamp = timestamp() ;
my %query = get_query_parameters() ;  # CGIが受け取ったデータの処理

#- ▼ 検索パラメータのセット
my $seq = flatsequence($query{'seq'}) || '' ;
#my $seq = flatsequence($query{'seq'}) || 'GGGTCCGGTTGCAATGCAA' ;  # テスト用にデフォルト値を設定
my $maxmismatch = $query{'MaxMismatch'} || 3 ;  # 指定なき場合は3ミスマッチまで
my $seqlength = $query{'TargetSize'} || length($seq) ;
my $maxoutputsize = $query{'MaxOutputSize'} || -1 ;
my $strand = 
	(not $query{'strand'}) ? 'both' :
	($query{'strand'} eq 'plus') ? 'plus' :
	($query{'strand'} eq 'minus') ? 'minus' :
		'both' ;  # plus, minus 以外は both にセット
#- ▲ 検索パラメータのセット

#- ▼ GGGenomeにアクセスして結果を取得
my $baseurl = 'https://gggenome.dbcls.jp/' ;
# 検索対象db
my $db = $query{'spe'} ;
# k
my $k = $maxmismatch ;
# strand_ggg
my $strand_ggg = '' ;
if ($strand eq 'both'){
	$strand_ggg = '' ;
} elsif ($strand eq 'plus'){
	$strand_ggg = '+' ;
} elsif ($strand eq 'minus'){
	$strand_ggg = '-' ;
}

# nogap
my $nogap = 'nogap' ;

# GGGenomeに対する検索
my $url = "${baseurl}/${db}/${k}/${strand_ggg}/${nogap}/${seq}.txt" ;
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime ;
if ($min == 59 && $sec >= 50){
	# 毎時00分のデーモン再起動を避ける
	sleep(20) ;
}
my $rawdata = get($url) ;
#- ▲ GGGenomeにアクセスして結果を取得

my @hits = parse_ggg_txt_for_si($rawdata) ;

@hits =
	($strand eq 'plus') ? grep {(split /\t/, $_)[1] eq 'plus'} @hits :
	($strand eq 'minus') ? grep {(split /\t/, $_)[1] eq 'minus'} @hits :
	@hits ;

my $offtargetlist_table = table_siDirectCore(@hits) ;
print_result_html($timestamp,$offtargetlist_table) ;

exit ;

# ====================
sub timestamp {
my ($sec,$min,$hour,$mday,$mon,$year) = (0,0,0,1,0,0) ;  # undefを返さないようにデフォルト値を代入しておく。
($sec,$min,$hour,$mday,$mon,$year) = localtime ;
my $timestamp = ($year+1900)*1e10 + ($mon+1)*1e8 + $mday*1e6 + $hour*1e4 + $min*100 + $sec ;
if ($timestamp =~ /(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/){
	return "$1-$2-$3 $4:$5:$6" ;
} else {
	return $timestamp ;
}
} ;
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
sub flatsequence {  # 塩基配列の整形
if (defined $_[0]){
	my $seq = $_[0] ;
	$seq =~ s/[^ATGCUNRYMKSWHBVD-]//gi ;
	return $seq ;
} else {
	return '' ;
}
} ;
# ====================
sub parse_siDirectCore {
my $rawdata = $_[0] || '' ;
my @hits = grep /^null\t/, split(/\n/,$rawdata) ;
@hits = sort {
	(split /\t/, $b)[3] cmp (split /\t/, $a)[3] ||  # 第1キー：match数
	(split /\t/, $b)[1] cmp (split /\t/, $a)[1] ||  # 第2キー：plus/minus
	(split /\t/, $a)[11] cmp (split /\t/, $b)[11]  # 第3キー：pattern
} @hits ;
return @hits ;
} ;
# ====================
sub table_siDirectCore {
my @hits = @_ ;
my $offtargetlist_table = "<table>\n" ;
foreach (@hits){
	my (undef,$strand,$mis,$match,undef,undef,undef,undef,$definition,$sseq,$tseq,$pattern) = split /\t/ ;
	my $match_bgcolor =
		($mis == 0) ? 'mis0' :
		($mis == 1) ? 'mis1' :
		($mis == 2) ? 'mis2' :
			'mis3' ;
	my $strand_bgcolor = ($strand eq 'plus') ? 'gs' : 'ps' ;
	$strand =~ s/plus/+/ ;
	$strand =~ s/minus/-/ ;
	$sseq = mismatch_red($sseq,$pattern) ;
	$definition = parse_definition_ggg($definition) ;
	$offtargetlist_table .= "<tr><td class=$match_bgcolor>$match\n" ;
	$offtargetlist_table .= "	<td class=$match_bgcolor>$strand\n" ;
	$offtargetlist_table .= "	<td><pre><span class=$strand_bgcolor>$tseq</span><br>$pattern<br>$sseq</pre>\n" ;
	$offtargetlist_table .= "	<td>$definition\n" ;
	$offtargetlist_table .= "</tr>\n" ;
}
$offtargetlist_table .= "</table>" ;
return $offtargetlist_table ;
} ;
#================================================================
# parse_ggg_txt_for_si
#================================================================
# gggenomeのtxt形式結果をパースして、siDirect用の結果に整形
# タブ区切り
# undef
# $strand	:	plusないしminus
#	$mis	：	mismatch数（背景色の設定に必要）
# $match	：	match数
# undef
# undef
# undef
# undef
#	$definition
# $sseq	:	similar seq（検索された配列）
# $tseq	：	target seq（クエリ配列）
# $pattern	:	一致不一致のパターン（|とXで標記、アライメントのミスマッチ文字を赤にするために必要）
# 結果はタブ区切りテキストの配列で下記のソートを行って返す
# 第1キー：match数の数値降順
# 第2キー：strand、plusが優先（文字列で降順ソート）
# 第3キー：patternで文字列昇順
#================================================================
sub parse_ggg_txt_for_si (){
	my @lines = split("\n", shift) ;
	my %strandDist2cnt = () ;
	my %strandDistSymbolAlign2cnt = () ;
	my %matchStrandAlignSymbol2def = () ;
	my %matchStrandAlignSymbol2pat = () ;
	my %matchStrandAlignSymbol2mis = () ;
	my %matchStrandAlignSymbol2sseq = () ;
	my %matchStrandAlignSymbol2tseq = () ;

	my %strandHash = () ;
	$strandHash{'+'} = "plus" ;
	$strandHash{'-'} = "minus" ;

	my $metaCnt = 0 ;
	my $lineCnt = 0 ;
	foreach my $res (@lines){
		# コメント行は読み飛ばし
		if ($res =~ /^#/){
			if ($res =~ /^# count:\t(\d+)/){
				$metaCnt += $1 ;
			}
			next ;
		}
		$lineCnt++ ;
		my @data = split("\t", $res) ;
		# シンボル抽出
		my $symbol = "" ;
		if ($data[0] =~ /^.+ \(([^ ]+)\)/){
				$symbol = $1 ;
		}
		# distance毎に、シンボルとターゲットのアライメントの組み合わせでユニーク化してカウントをとっていく
		my $strand = $strandHash{$data[1]} ;
		my $align = $data[8] ;  # target(similar sequence)のアライメント配列
		my $match = $data[11] ;
		my $dist = $data[12] ;  # nogapなのでdistance = mismatch
		my $pattern = $data[9] ;
		$pattern =~ s/ /X/g ;
		my $patternAlign = $pattern . $align ;
		my $def = $data[0] ;

		if ($matchStrandAlignSymbol2def{$match}{$strand}{$patternAlign}{$symbol}){
			$matchStrandAlignSymbol2def{$match}{$strand}{$patternAlign}{$symbol} .= "\n" ;
		}
		$matchStrandAlignSymbol2def{$match}{$strand}{$patternAlign}{$symbol} .= $def ;
		$matchStrandAlignSymbol2pat{$match}{$strand}{$patternAlign}{$symbol} = $pattern ;
		$matchStrandAlignSymbol2mis{$match}{$strand}{$patternAlign}{$symbol} = $dist ;
		$matchStrandAlignSymbol2sseq{$match}{$strand}{$patternAlign}{$symbol} = $align ;
		$matchStrandAlignSymbol2tseq{$match}{$strand}{$patternAlign}{$symbol} = $data[7] ;
	}
	# 念のため件数の一致も確認？
	if ($metaCnt != $lineCnt){
		warn("meta count ($metaCnt) does not match with results ($lineCnt)") ;
	}
	# 結果を配列にまとめる
	my @resAry = () ;
	foreach my $match (sort {$b <=> $a} keys %matchStrandAlignSymbol2def){
		foreach my $strand (sort {$b cmp $a} keys %{$matchStrandAlignSymbol2def{$match}}){
			foreach my $patternAlign (sort {$a cmp $b} keys %{$matchStrandAlignSymbol2def{$match}{$strand}}){  # 同じミスマッチ位置でも配列が異なる場合には区別するのでパターンとアライメントで一意に決める
				foreach my $symbol (sort {$a cmp $b} keys %{$matchStrandAlignSymbol2def{$match}{$strand}{$patternAlign}}){
					push @resAry, join(
						"\t",
						'',
						$strand,
						$matchStrandAlignSymbol2mis{$match}{$strand}{$patternAlign}{$symbol},
						$match,
						'',
						'',
						'',
						'',
						$symbol . "\n" . $matchStrandAlignSymbol2def{$match}{$strand}{$patternAlign}{$symbol},
						$matchStrandAlignSymbol2sseq{$match}{$strand}{$patternAlign}{$symbol},
						$matchStrandAlignSymbol2tseq{$match}{$strand}{$patternAlign}{$symbol},
						$matchStrandAlignSymbol2pat{$match}{$strand}{$patternAlign}{$symbol}
					) ;
				}
			}
		}
	}
	return @resAry ;
} ;
# ====================
sub mismatch_red {
my ($sseq,$pattern,undef) = @_ ;
if (length($sseq) eq length($pattern)){
	my @sseq = split //, $sseq ;
	my @pattern = split //, $pattern ;
	foreach (1..length($pattern)){
		if ($pattern[$_-1] eq 'X'){
			$sseq[$_-1] = "<font color=red>$sseq[$_-1]</font>" ;
		}
	}
	$sseq = join "", @sseq ;
}
return $sseq ;
} ;
# ====================
sub parse_definition {
my $def = $_[0] ;
my $out =
	($def =~ s|^skip/sz=(1)/?||) ? "<font color=#999966>unaligned-skip($1)</font>" :
	($def =~ s|^exon/sg=\d+/eg=\d+/sz=(\d+)/?||) ? "<font color=#999966>exon($1)</font>" :
	($def =~ s|^intron/sg=\d+/eg=\d+/sz=(\d+)/?||) ? "<font color=#999966>exon-exon junction($1)</font>" :
	($def =~ s|^notalign/sz=(1)/?||) ? "<font color=#999966>unaligned-notalign($1)</font>" :
	"[parse_error]" ;
while (
	$def =~ s/^(?:sq=)?gi\|\d+\|ref\|(.*?)\|\s*(.*?)(\Z|\/sq=)// and $out .= "<br>" . "<a target=\"_blank\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$1\">$1</a> | $2" or
	$def =~ s/^(?:sq=)?gnl\|UG\|(.*?)(?:\s+\/.*)?\/gb=(.*?)\s.*?\/gi=(.*?)\s.*?(\Z|\/sq=)// and $out .= "<br>" . "<a target=\"_blank\" href=\"https://www.ncbi.nlm.nih.gov/nuccore/$3\">$2</a> | $1"
){}
if ($def){
	$out .= "<br>" . $def ;
}
return $out ;
} ;
# ====================
sub parse_definition_ggg {
my @defAry = split("\n", $_[0]) ;
my $out = "<font color=#999966>$defAry[0]</font>" ;
for (my $i = 1 ; $i < @defAry ; $i++){
	$defAry[$i] =~ s/^(\S+) /<a target="_blank" href="https:\/\/www.ncbi.nlm.nih.gov\/nuccore\/$1">$1<\/a> | / ;
	$out .= "<br />" . $defAry[$i] ;
}

return $out ;
} ;
# ====================
sub print_result_html {
my $timestamp = $_[0] || '' ;
my $text = $_[1] || '' ;
print 'Content-type: text/html; charset=utf-8

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html lang=en>

<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="Content-Style-Type" content="text/css">
<meta name="author" content="Yuki Naito et al.">
<title>siDirect</title>
<style type="text/css">
<!--
/* font */
	p,table,div,h1,h2,h3 { font-family:verdana,arial,helvetica,sans-serif }
	p,table,div { font-size:7pt }
	pre,textarea { font-family:courier,monospace; font-size:8pt }
/* hyperlink */
	a:link,a:visited {
		text-decoration:none;
		color:#004080;
	}
	a:hover {
		text-decoration:none;
		color:red;
	}
/* siDirect */
	td { border:1px solid gray; vertical-align:baseline }
	.mis0 { text-align:center; font-weight:bold; background-color:#AAAAAA }
	.mis1 { text-align:center; font-weight:bold; background-color:#CCCCCC }
	.mis2 { text-align:center; font-weight:bold; background-color:#DDDDDD }
	.mis3 { text-align:center; font-weight:bold; background-color:#EEEEEE }
	.gs { background-color:#BBFFBB }
	.ps { background-color:#BBBBFF }
-->
</style>
</head>

<body>

<a href="/"><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt">
<font size=5>siDirect </font><font size=4>version 2.1<small>&beta;</small> </font>result page.
<a target="_blank" href="doc/"><img src="help.png" alt="Help" width=36 height=15 border=0></a>
</div>

<hr><!-- __________________________________________________ -->

<p><font color=gray>' . $timestamp . ',  siDirect v2.1&beta;</font></p>

<h3>Similar Sequences</h3>

' . $text . '

</body>
</html>
' ;
exit ;
} ;
# ====================
