#!/usr/bin/perl

# 2008-08-26.y-naito
# 2008-08-27.y-naito
# 2008-08-28.y-naito
# 2008-09-01.y-naito
# 2008-09-02.y-naito
# 2008-09-17.y-naito
# 2008-09-24.y-naito（マウス追加）
# 2008-12-08.y-naito（@siRNAlistの構造をハッシュへのリファレンスに変更）
# 2008-12-09.y-naito（オプションを追加）
# 2008-12-10.y-naito（オプションを追加）
# 2008-12-11.y-naito（オプションを追加）
# 2008-12-12.y-naito（オプションを追加）
# 2008-12-13.y-naito（Excel出力を追加）
# 2008-12-22.y-naito（タイムスタンプ出力）
# 2008-12-26.y-naito（precomputed siRNAにない場合の対応）
# 2008-12-29.y-naito（Ui-Tei + Reynolds ^ Amarzguiouiのcombinationを実装）
# 2009-01-20.y-naito（ラットNRDBを追加）
# 2009-02-14.y-naito（combinedルール変更）
# 2009-03-16.y-naito（siRNAが1個も設計できなかったときに別画面を表示）
# 2009-04-14.y-naito（sense側のseed Tmを考慮）
# 2009-06-30.y-naito（off-target検索をguide 2-20, passenger 4-22に変更）
# 2009-07-01.y-naito（ミスマッチ耐性ではなくseed Tmで色をつける）
# 2009-08-20.y-naito（関数を整理）
# 2009-08-20.y-naito（23ntの範囲にNがあるものは候補からはずす）
# 2009-10-19.y-naito（クエリのオプション一覧を表示）
# 2009-10-19.y-naito（GC含量で0未満や100より大の数値を置き換える）
# 2009-10-19.y-naito（TmとGC含量の入力に0をプラスして数値にする：-0→0、0.30→0.3）
# 2009-10-19.y-naito（RNAオリゴ配列を表示）
# 2010-05-12.y-naito（DBが落ちている場合にspe=noneの結果を表示）
# 2010-05-12.y-naito（配列名を結果ページに表示する）

use warnings ;
use strict ;

my $siDirect_top_url = 'http://siDirect2.RNAi.jp/' ;

my $timestamp = timestamp() ;
my %query = get_query_parameters() ;  # CGIが受け取ったデータの処理
my %option = set_options(%query) ;  # 検索オプションを決定

unless ($query{'yourSeq'}){
	print_redirect_html('Please input nucleotide sequence.',$siDirect_top_url) ;
	exit ;
}

($option{'targetname'},my $targetsequence) = readFASTA($query{'yourSeq'}) ;
unless ($targetsequence){
	print_redirect_html('Please input nucleotide sequence.',$siDirect_top_url) ;
	exit ;
}

# ▼ Tm計算のための準備
my ($rna_dh_ref,$rna_ds_ref) = get_nearest_neighbor_param() ;
my %rna_dh = %{$rna_dh_ref} ;
my %rna_ds = %{$rna_ds_ref} ;
# ▲ Tm計算のための準備

# ▼ DBに接続
use DBI ;
my $dbh = DBI->connect(
	'DBI:mysql:sidirect2009','root','designsi',
	# { RaiseError => 1, AutoCommit => 0 }
	{ AutoCommit => 0 }
) or ($option{'hs'},$option{'mm'},$option{'rn'}) = (0,0,0) ;  # DBが落ちている場合はspe=noneの結果を表示（2010-05-12変更）
# ▲ DBに接続

if (my @siRNAlist = select_sirna($targetsequence,\%option)){
	my $siRNA_table_html = print_siRNA_table(\@siRNAlist,\%option) ;
	my $siRNA_table_txt = print_siRNA_txt(\@siRNAlist,\%option) ;
	my $graphical_view_html = print_graphical_view($targetsequence,undef,undef,\@siRNAlist) ;
	my $query_summary_html = print_Qsummary_html($targetsequence,\%option) ;
	print_result_html($siRNA_table_html,$graphical_view_html,$siRNA_table_txt,$timestamp,$query_summary_html) ;
} else {
	my $query_summary_html = print_Qsummary_html($targetsequence,\%option) ;
	print_no_siRNA_html($timestamp,$query_summary_html) ;  # siRNAがひとつも設計できなかった場合にトラブルシューティングを案内。
}

$dbh->disconnect ;
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
sub set_options {
my %query = @_ ;
# ▼ 効くsiRNAの選択オプション
# 優先順位：
# (combined rule指定) ＞ (show all siRNAs) ＞ (individual rule指定)、3つのうち優先度が高いもののみ有効。
# combined rule指定のなかでは、(UorRorA) ＞ (UorRA) ＞ (URA)、3つのうち優先度が高いもののみ有効。
$query{'uitei'} and $option{'uitei'} = 1 ;
$query{'reynolds'} and $option{'reynolds'} = 1 ;
$query{'amarzguioui'} and $option{'amarzguioui'} = 1 ;
$query{'UorRorA'} and ($option{'uitei'},$option{'reynolds'},$option{'amarzguioui'},$option{'UorRorA'}) = (1,1,1,1) ;
$query{'UorRA'} and ($option{'uitei'},$option{'reynolds'},$option{'amarzguioui'},$option{'UorRA'}) = (1,1,1,1) ;
$query{'URA'} and ($option{'uitei'},$option{'reynolds'},$option{'amarzguioui'},$option{'URA'}) = (1,1,1,1) ;
# ALL: デバッグ用、効かないものも含めて全siRNAを表示。
$query{'ALL'} and ($option{'uitei'},$option{'reynolds'},$option{'amarzguioui'},$option{'ALL'}) = (1,1,1,1) ;
# ▲ 効くsiRNAの選択オプション
# ▼ 19mer/3mm検索のprecomputed DB照会オプション
if ($query{'spe'} and $query{'spe'} eq 'hs'){
	$option{'hs'} = 1 ;  # ヒトNRDBに対する3mm検索
} elsif ($query{'spe'} and $query{'spe'} eq 'mm'){
	$option{'mm'} = 1 ;  # マウスNRDBに対する3mm検索
} elsif ($query{'spe'} and $query{'spe'} eq 'rn'){
	$option{'rn'} = 1 ;  # ラットNRDBに対する3mm検索
}
$query{'hidenonspe'} and $option{'hidenonspe'} = 1 ;
$query{'hitcount'} and $option{'hitcount'} = 1 ;
# ▲ 19mer/3mm検索のprecomputed DB照会オプション
# ▼ seed部分のTmによる選択オプション
if ($query{'seedTm'}){
	$option{'seedTm'} = 1 ;
}
if (defined $query{'seedTmMax'} and $query{'seedTmMax'} =~ /^\s*(-?\s*?\d+(\.\d*)?)\s*$/){
	($option{'seedTmMax'} = $1) =~ s/\s//g ;  # マイナスと数字との間のスペースを除去
}  # 数値が入るseedTmMax、posStart、posEnd、percentGCMin、percentGCMaxは、ゼロが入力されることを考慮してdefinedで判定
# ▲ seed部分のTmによる選択オプション
# ▼ Target rangeによる選択オプション
$query{'pos'} and $option{'pos'} = 1 ;
if (defined $query{'posStart'} and $query{'posStart'} =~ /^\s*(\d+)\s*$/){
	$option{'posStart'} = $1 ;
}
if (defined $query{'posEnd'} and $query{'posEnd'} =~ /^\s*(\d+)\s*$/){
	$option{'posEnd'} = $1 ;
}
# ▲ Target rangeによる選択オプション
# ▼ Gの連続・Cの連続による選択オプション
$query{'consGC'} and $option{'consGC'} = 1 ;
if ($query{'consGCmax'} and $query{'consGCmax'} =~ /^[4-7]$/){
	$option{'consGCmax'} = $query{'consGCmax'} ;
}
# ▲ Gの連続・Cの連続による選択オプション
# ▼ Aの連続・Tの連続による選択オプション
$query{'consAT'} and $option{'consAT'} = 1 ;
if ($query{'consATmax'} and $query{'consATmax'} =~ /^[4-7]$/){
	$option{'consATmax'} = $query{'consATmax'} ;
}
# ▲ Aの連続・Tの連続による選択オプション
# ▼ GC含量による選択オプション
$query{'percentGC'} and $option{'percentGC'} = 1 ;
if (defined $query{'percentGCMin'} and $query{'percentGCMin'} =~ /^\s*(-?\s*?\d+(\.\d*)?)\s*$/){
	($option{'percentGCMin'} = $1) =~ s/\s//g ;  # マイナスと数字との間のスペースを除去
	$option{'percentGCMin'} =
		($option{'percentGCMin'} > 100) ? 100 :
		($option{'percentGCMin'} < 0) ? 0 :
			$option{'percentGCMin'} ;
}
if (defined $query{'percentGCMax'} and $query{'percentGCMax'} =~ /^\s*(-?\s*?\d+(\.\d*)?)\s*$/){
	($option{'percentGCMax'} = $1) =~ s/\s//g ;  # マイナスと数字との間のスペースを除去
	$option{'percentGCMax'} =
		($option{'percentGCMax'} > 100) ? 100 :
		($option{'percentGCMax'} < 0) ? 0 :
			$option{'percentGCMax'} ;
}
# ▲ GC含量による選択オプション
# ▼ Custom patternによる選択オプション
$query{'custom'} and $option{'custom'} = 1 ;
if ($query{'customPattern'} and $query{'customPattern'} =~ /^\s*([ATGCUNRYMKSWHBVD]{0,23})\s*$/i){  # 23文字以下の塩基構成文字だけからなる文字列
	$option{'customPattern'} = uc $query{'customPattern'} ;
}
# ▲ Custom patternによる選択オプション
# ▼ Exclude patternによる選択オプション
$query{'exclude'} and $option{'exclude'} = 1 ;
if ($query{'excludePattern'} and $query{'excludePattern'} =~ /^\s*([ATGCUNRYMKSWHBVD]{0,23})\s*$/i){  # 23文字以下の塩基構成文字だけからなる文字列
	$option{'excludePattern'} = uc $query{'excludePattern'} ;
}
# ▲ Exclude patternによる選択オプション
# ▼ すべての条件にマッチするsiRNAだけを残すオプション
$query{'hide'} and $option{'hide'} = 1 ;
# ▲ すべての条件にマッチするsiRNAだけを残すオプション
return %option ;
} ;
# ====================
sub readFASTA {  # 塩基配列の整形
if (defined $_[0]){
	(my $seq = $_[0]) =~ s/\A\s*(>[\ \t]*(.*?)[\ \t]*$)?//m ;
	my $name = $2 || '' ;
	$seq =~ s/>.*//sm ;  # 2個目以降の配列を削除
	$seq = flatsequence($seq) ;
	return ($name,$seq) ;
} else {
	return ('','') ;
}
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
sub print_Qsummary_html {
if (not defined $_[0]){
	print_error_html('ERROR : print_Qsummary_html() : target sequence not defined.') ;
	exit ;
}
my $targetseq = $_[0] ;
if ($targetseq =~ /[^ATGCUNRYMKSWHBVD-]/i){  # 塩基構成文字以外が含まれていないかチェック
	print_error_html('ERROR : print_Qsummary_html() : target sequence contains invalid characters.') ;
	exit ;
}
my %option ;
if (defined $_[1] and ref $_[1] eq 'HASH'){
	%option = %{$_[1]} ;
}
# ▼ オプションを抽出
my $targetname = $option{'targetname'} || '' ;
my $querylength = length($targetseq) ;
my $rule = '' ;
$rule .= $query{'uitei'} ? 'Ui-Tei, ' : '' ;
$rule .= $query{'reynolds'} ? 'Reynolds, ' : '' ;
$rule .= $query{'amarzguioui'} ? 'Amarzguioui, ' : '' ;
$rule =~ s/,\ *$// ;
$rule = 
	$option{'UorRorA'} ? 'Ui-Tei + Reynolds + Amarzguioui' :
	$option{'UorRA'} ? 'Ui-Tei + Reynolds &times; Amarzguioui' :
	$option{'URA'} ? 'Ui-Tei &times; Reynolds &times; Amarzguioui' :
	$option{'ALL'} ? 'Show all siRNAs' :
	$rule ? $rule :
		'(not selected)' ;
my $seedTmMax = ($option{'seedTm'} and defined $option{'seedTmMax'} and not $option{'seedTmMax'} eq '') ? $option{'seedTmMax'} + 0 . '&deg;C' : '(blank)' ;
my $posStart = ($option{'pos'} and $option{'posStart'}) ? $option{'posStart'} : '(blank)' ;
my $posEnd = ($option{'pos'} and $option{'posEnd'}) ? $option{'posEnd'} : '(blank)' ;
my $spe =
	$option{'hs'} ? 'Homo sapiens non-redundant database' :
	$option{'hs'} ? 'Mus musculus non-redundant database' :
	$option{'hs'} ? 'Rattus norvegicus non-redundant database' :
		'none' ;
my $consGCmax = ($option{'consGC'} and $option{'consGCmax'} and $option{'consGCmax'} =~ /^[4-7]$/) ? "$option{'consGCmax'} nt" : '(blank)' ;
my $consATmax = ($option{'consAT'} and $option{'consATmax'} and $option{'consATmax'} =~ /^[4-7]$/) ? "$option{'consATmax'} nt" : '(blank)' ;
my $percentGCMin = ($option{'percentGC'} and defined $option{'percentGCMin'} and not $option{'percentGCMin'} eq '') ? $option{'percentGCMin'} + 0 . '%' : '(blank)' ;
my $percentGCMax = ($option{'percentGC'} and defined $option{'percentGCMax'} and not $option{'percentGCMax'} eq '') ? $option{'percentGCMax'} + 0 . '%' : '(blank)' ;
my $customPattern = ($option{'custom'} and $option{'customPattern'}) ? $option{'customPattern'} : '(blank)' ;
my $excludePattern = ($option{'exclude'} and $option{'excludePattern'}) ? $option{'excludePattern'} : '(blank)' ;
# ▲ オプションを抽出
# ▼ HTMLを生成
my $html = "<p>\n" ;
$html .= "	<b>Query name:</b> $targetname<br>\n" ;  # 配列名を結果ページに表示する（2010-05-12変更）
$html .= "	<b>Query sequence:</b> $querylength bp<br>\n" ;
$html .= "	<b>Functional siRNA selection:</b> $rule<br>\n" ;
$option{'seedTm'} and $html .= "	<b>Seed-duplex stability - Max Tm:</b> $seedTmMax<br>\n" ;
$html .= "	<b>Specificity check:</b> $spe<br>\n" ;
$option{'pos'} and $html .= "	<b>Target range:</b> $posStart - $posEnd<br>\n" ;
$option{'consGC'} and $html .= "	<b>Avoid contiguous G's or C's:</b> $consGCmax<br>\n" ;
$option{'consAT'} and $html .= "	<b>Avoid contiguous A's or T's:</b> $consATmax<br>\n" ;
$option{'percentGC'} and $html .= "	<b>GC content:</b> $percentGCMin - $percentGCMax<br>\n" ;
$option{'custom'} and $html .= "	<b>Custom pattern:</b> $customPattern<br>\n" ;
$option{'exclude'} and $html .= "	<b>Exclude pattern:</b> $excludePattern<br>\n" ;
$html .= "</p>\n" ;
# ▲ HTMLを生成
return $html ;
} ;
# ====================
sub select_sirna {
# usage:
# my @siRNAlist = select_sirna($targetsequence,\%option) ;
# siRNAリストを出力する。
# @siRNAlistの要素は、ハッシュのリファレンスになっている。
# ex:
# 'startpos' => 1
# 'endpos' => 23
# 'si23' => 'agctcaagggccaaggcaagtcg'
# 'name' => '1-23'
# 'efficacy' => 'URA'
# 'seed_tm' => 10.0
# 'hide' => 0

if (not defined $_[0]){
	print_error_html('ERROR : select_sirna() : target sequence not defined.') ;
	exit ;
}
my $targetseq = $_[0] ;
if ($targetseq =~ /[^ATGCUNRYMKSWHBVD-]/i){  # 塩基構成文字以外が含まれていないかチェック
	print_error_html('ERROR : select_sirna() : target sequence contains invalid characters.') ;
	exit ;
}
my %option ;
if (defined $_[1] and ref $_[1] eq 'HASH'){
	%option = %{$_[1]} ;
}
# ▼ 効くsiRNAの選択
my @siRNAlist ;
foreach (1..length($targetseq)-22){
	my $si23 = substr($targetseq,$_ - 1,23) ;
	my $efficacy = '' ;
	if ($option{'uitei'} and uitei_chk($si23) > 0){
		$efficacy .= 'U' ;
	}
	if ($option{'reynolds'} and reynolds_chk($si23) > 0){
		$efficacy .= 'R' ;
	}
	if ($option{'amarzguioui'} and amar_chk($si23) > 0){
		$efficacy .= 'A' ;
	}
	if ($efficacy and $si23 =~ /^[atugcATUGC]{23}$/ or $option{'ALL'}){  # OH内にNなどを許容しない（2009-08-20変更）
		my $endpos = $_ + 22 ;
		#- @siRNAlist：ハッシュのリファレンスを格納したアレイ（2008-12-08変更）
		push @siRNAlist, {
			'startpos' => $_,
			'endpos' => $endpos,
			'si23' => $si23,
			'name' => "$_-$endpos",
			'efficacy' => $efficacy,
			'hide' => 0  # このsiRNAを結果画面に表示するか否か
		}
	}
}
# ▲ 効くsiRNAの選択
# ▼ 効くsiRNA・combined ruleによる選択
if ($option{'UorRorA'}){
	foreach (@siRNAlist){
		if (not $$_{'efficacy'} =~ /U|R|A/){
			$$_{'hide'} = 1 ;
		}
	}
}
if ($option{'UorRA'}){
	foreach (@siRNAlist){
		if (not $$_{'efficacy'} =~ /U|RA/){
			$$_{'hide'} = 1 ;
		}
	}
}
if ($option{'URA'}){
	foreach (@siRNAlist){
		if (not $$_{'efficacy'} eq 'URA'){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ 効くsiRNA・combined ruleによる選択
# ▼ 19mer/3mm検索のprecomputed DB照会
if ($option{'hs'}){
	foreach (@siRNAlist){
		#my ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3) = sidirect2_hs($$_{'si23'}) ;
		my ($p0,undef,$p1,undef,$p2,undef,$p3,undef) = sidirect2_hs(substr($$_{'si23'},1,19)) ;  # 2009-06-30変更
		my (undef,$m0,undef,$m1,undef,$m2,undef,$m3) = sidirect2_hs(substr($$_{'si23'},3,19)) ;  # 2009-06-30変更
		my $mt_plus = mt_plus($p0,$p1,$p2,$p3) ;  # 2009-06-30変更
		my $mt_minus = mt_minus($m0,$m1,$m2,$m3) ;  # 2009-06-30変更
		@$_{'p0','m0','p1','m1','p2','m2','p3','m3','mt_plus','mt_minus'} = ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,$mt_plus,$mt_minus) ;
		if ($option{'hidenonspe'} and not ($mt_plus and $mt_plus >= 2 and $mt_minus and $mt_minus >= 2)){  # mt=2までOKとする。2009-07-01変更
			$$_{'hide'} = 1 ;
		}
	}
} elsif ($option{'mm'}){
	foreach (@siRNAlist){
		#my ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3) = sidirect2_mm($$_{'si23'}) ;
		my ($p0,undef,$p1,undef,$p2,undef,$p3,undef) = sidirect2_mm(substr($$_{'si23'},1,19)) ;  # 2009-06-30変更
		my (undef,$m0,undef,$m1,undef,$m2,undef,$m3) = sidirect2_mm(substr($$_{'si23'},3,19)) ;  # 2009-06-30変更
		my $mt_plus = mt_plus($p0,$p1,$p2,$p3) ;  # 2009-06-30変更
		my $mt_minus = mt_minus($m0,$m1,$m2,$m3) ;  # 2009-06-30変更
		@$_{'p0','m0','p1','m1','p2','m2','p3','m3','mt_plus','mt_minus'} = ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,$mt_plus,$mt_minus) ;
		if ($option{'hidenonspe'} and not ($mt_plus and $mt_plus >= 2 and $mt_minus and $mt_minus >= 2)){  # mt=2までOKとする。2009-07-01変更
			$$_{'hide'} = 1 ;
		}
	}
} elsif ($option{'rn'}){
	foreach (@siRNAlist){
		#my ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3) = sidirect2_rn($$_{'si23'}) ;
		my ($p0,undef,$p1,undef,$p2,undef,$p3,undef) = sidirect2_rn(substr($$_{'si23'},1,19)) ;  # 2009-06-30変更
		my (undef,$m0,undef,$m1,undef,$m2,undef,$m3) = sidirect2_rn(substr($$_{'si23'},3,19)) ;  # 2009-06-30変更
		my $mt_plus = mt_plus($p0,$p1,$p2,$p3) ;  # 2009-06-30変更
		my $mt_minus = mt_minus($m0,$m1,$m2,$m3) ;  # 2009-06-30変更
		@$_{'p0','m0','p1','m1','p2','m2','p3','m3','mt_plus','mt_minus'} = ($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,$mt_plus,$mt_minus) ;
		if ($option{'hidenonspe'} and not ($mt_plus and $mt_plus >= 2 and $mt_minus and $mt_minus >= 2)){  # mt=2までOKとする。2009-07-01変更
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ 19mer/3mm検索のprecomputed DB照会
# ▼ seed部分のTmによる選択
if ($option{'seedTm'}){
	foreach (@siRNAlist){
		$$_{'seed_tm'} = seedTm($$_{'si23'}) ;
		$$_{'seed_tm_sense'} = seedTm(complementaryDNA($$_{'si23'})) ;
		# ▼ Tmで色分けするよう変更（2009-07-01）
		$$_{'color'} =
			($$_{'seed_tm'} < 10 and $$_{'seed_tm_sense'} < 10) ? 3 :
			($$_{'seed_tm'} < 15 and $$_{'seed_tm_sense'} < 15) ? 2 :
			($$_{'seed_tm'} < 21.5 and $$_{'seed_tm_sense'} < 21.5) ? 1 :
				0 ;
		# ▲ Tmで色分けするよう変更（2009-07-01）
		if (defined $option{'seedTmMax'} and not $option{'seedTmMax'} eq '' and
			($option{'seedTmMax'} < $$_{'seed_tm'} or $option{'seedTmMax'} < $$_{'seed_tm_sense'})  # 両側のseed Tmをチェック
		){  # 値がゼロの場合を考慮
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ seed部分のTmによる選択
# ▼ Target rangeによる選択
if ($option{'pos'}){
	foreach (@siRNAlist){
		$$_{'pos'} = ((defined $option{'posStart'} and $option{'posStart'} > $$_{'startpos'}) or
			(defined $option{'posEnd'} and $option{'posEnd'} < $$_{'endpos'})) ? 0 : 1 ;  # 下限・上限とも空欄を許容
		if ($option{'hide'} and $$_{'pos'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ Target rangeによる選択
# ▼ Gの連続・Cの連続による選択
if ($option{'consGC'} and $option{'consGCmax'} and $option{'consGCmax'} =~ /^[4-7]$/){
	foreach (@siRNAlist){
		$$_{'consGC'} = ($$_{'si23'} =~ /(G{$option{'consGCmax'}})|(C{$option{'consGCmax'}})/i) ? 0 : 1 ;
		if ($option{'hide'} and $$_{'consGC'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ Gの連続・Cの連続による選択
# ▼ Aの連続・Tの連続による選択
if ($option{'consAT'} and $option{'consATmax'} and $option{'consATmax'} =~ /^[4-7]$/){
	foreach (@siRNAlist){
		$$_{'consAT'} = ($$_{'si23'} =~ /(A{$option{'consATmax'}})|([TU]{$option{'consATmax'}})/i) ? 0 : 1 ;
		if ($option{'hide'} and $$_{'consAT'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ Aの連続・Tの連続による選択
# ▼ GC含量による選択
if ($option{'percentGC'}){
	foreach (@siRNAlist){
		$$_{'percentGC'} = ((defined $option{'percentGCMin'} and $option{'percentGCMin'} > gc_content($$_{'si23'})) or  # ゼロの可能性も考慮
		(defined $option{'percentGCMax'} and $option{'percentGCMax'} < gc_content($$_{'si23'}))) ? 0 : 1 ;  # 下限・上限とも空欄を許容
		if ($option{'hide'} and $$_{'percentGC'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ GC含量による選択
# ▼ Custom patternによる選択
if ($option{'custom'} and $option{'customPattern'} and my $customPattern_regexp = iub2regexp($option{'customPattern'})){
	# ↑変換後の正規表現（$customPattern_regexp）が空白でないこともチェック
	foreach (@siRNAlist){
		$$_{'custom'} = ($$_{'si23'} =~ /${customPattern_regexp}/i) ? 1 : 0 ;
		if ($option{'hide'} and $$_{'custom'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ Custom patternによる選択
# ▼ Exclude patternによる選択
if ($option{'exclude'} and $option{'excludePattern'} and my $excludePattern_regexp = iub2regexp($option{'excludePattern'})){
	# ↑変換後の正規表現（$excludePattern_regexp）が空白でないこともチェック
	foreach (@siRNAlist){
		$$_{'exclude'} = ($$_{'si23'} =~ /${excludePattern_regexp}/i) ? 0 : 1 ;  # Custom patternの判定を逆転しただけ
		if ($option{'hide'} and $$_{'exclude'} == 0){
			$$_{'hide'} = 1 ;
		}
	}
}
# ▲ Exclude patternによる選択（Custom patternの反転を逆転しただけ）
# ▼ 'hide' => 1 のsiRNAを除く
my @new_siRNAlist ;
foreach (@siRNAlist){
	if (not $$_{'hide'}){
		push @new_siRNAlist, $_ ;
	}
}
@siRNAlist = @new_siRNAlist ;
# ▲ 'hide' => 1 のsiRNAを除く
return @siRNAlist ;
} ;
# ====================
sub check_19mer {
# 19bpの規定配列を抽出して返す。失敗したら空文字列を返す。
if (defined $_[0] and my $seq = flatsequence($_[0])){
	my $si19 ;
	if (length $seq == 23 or length $seq == 21){
		$si19 = substr($seq,2,19) ;
	} elsif (length $seq == 19){
		$si19 = $seq ;
	}
	if (defined $si19 and $si19 =~ /^[atugcATUGC]{19}$/){
		return $si19 ;
	}
}
return '' ;
} ;
# ====================
sub uitei_chk {  # Ui-Tei Guidelineを満たすかチェック。(-1,0,1,2)
if (my $si19 = check_19mer($_[0])){
	my $at7 = (substr($si19,12,7) =~ tr/atuATU/atuATU/) ;
	if ($si19 =~ /^[gc][atugc]{17}[atu]$/i and  # 19bpにambiguousな文字を含まない
		$at7 >= 4 and
		not $si19 =~ /[gc]{10}/i  # ds部分で10 base以上の連続GC stretchを排除
	){
		if ($at7 >= 5){
			return 2 ;  # class Ia → 2
		} else {
			return 1 ;  # class Ib → 2
		}
	} else {
		return 0 ;
	}
} else {
	return -1 ;
}
} ;
# ====================
sub reynolds_chk {  # Reynoldsのルールを満たすかチェック。(-1,0,6〜8)
if (my $si19 = check_19mer($_[0])){
	my $pts = 0 ;
	# criteria I
	my $gc = ($si19 =~ tr/gcGC/gcGC/) ;
	if (7 <= $gc and $gc <= 10){ $pts ++ }
	# criteria II
	my $at5 = (substr($si19,14,5)  =~ tr/atuATU/atuATU/) ;
	$pts += $at5 ;
	# criteria IV + VII
	if ($si19 =~ /^.{18}a$/i){ $pts ++ }
	elsif ($si19 =~ /^.{18}[gc]$/i){ $pts -- }
	# criteria V
	if ($si19 =~ /^.{2}a.{16}$/i){ $pts ++ }
	# criteria VI
	if ($si19 =~ /^.{9}[tu].{9}$/i){ $pts ++ }
	# criteria VIII
	if ($si19 =~ /^.{12}g.{6}$/i){ $pts -- }
	# criteria III
	if ($pts == 5){
		#tm() ;  # 実装できない
	}
	# まとめ
	if ($pts >= 6){
		return $pts ;  # criteria IIIが実装不能なので除外、5以上を効くと判定。（6以上？）
	} else {
		return 0 ;
	}
} else {
	return -1 ;
}
} ;
# ====================
sub amar_chk {  # Amarzguiouiのルールを満たすかチェック。(-1,0,1)
if (my $si19 = check_19mer($_[0])){
	my $pts = 0 ;
	# delta_T3
	my $right3 = (substr($si19,16,3)  =~ tr/atuATU/atuATU/) ;
	my $left3 = (substr($si19,0,3)  =~ tr/atuATU/atuATU/) ;
	my $delta_T3 = $right3 - $left3 ;
	$pts += $delta_T3 ;
	# S1 + U1
	if ($si19 =~ /^[gc].{18}$/i){ $pts ++ }
	elsif ($si19 =~ /^[tu].{18}$/i){ $pts -- }
	# A6
	if ($si19 =~ /^.{5}a.{13}$/i){ $pts ++ }
	# G19 + W19
	if ($si19 =~ /^.{18}g$/i){ $pts -- }
	elsif ($si19 =~ /^.{18}[atu]$/i){ $pts ++ }
	# まとめ
	if ($pts >= 3){
		return $pts ;
	} else {
		return 0 ;
	}
} else {
	return -1 ;
}
} ;
# ====================
sub sidirect2_hs {
if (my $si19 = check_19mer($_[0])){
	if (my @sql_kekka = sql_select("
			SELECT
				p0,m0,p1,m1,p2,m2,p3,m3
			FROM
				sidirect2_hs
			WHERE
				targetseq = \"$si19\"
		")){
		my @offtargets = split /\t/, $sql_kekka[0] ;
		return @offtargets ;
	}
}
return ('','','','','','','','') ;
} ;
# ====================
sub sidirect2_mm {
if (my $si19 = check_19mer($_[0])){
	if (my @sql_kekka = sql_select("
			SELECT
				p0,m0,p1,m1,p2,m2,p3,m3
			FROM
				sidirect2_mm
			WHERE
				targetseq = \"$si19\"
		")){
		my @offtargets = split /\t/, $sql_kekka[0] ;
		return @offtargets ;
	}
}
return ('','','','','','','','') ;
} ;
# ====================
sub sidirect2_rn {
if (my $si19 = check_19mer($_[0])){
	if (my @sql_kekka = sql_select("
			SELECT
				p0,m0,p1,m1,p2,m2,p3,m3
			FROM
				sidirect2_rn
			WHERE
				targetseq = \"$si19\"
		")){
		my @offtargets = split /\t/, $sql_kekka[0] ;
		return @offtargets ;
	}
}
return ('','','','','','','','') ;
} ;
# ====================
sub sql_select {  # @result = select("SELECT * FROM refseq") でタブ区切りテキスト出力。
my $select = $_[0] || return ;
my $sth = $dbh->prepare($select) ;
$sth->execute ;
my @result ;
while (my @hit = $sth->fetchrow_array){
	my $hit = join "\t", @hit ;
	push @result, $hit ;
}
$sth->finish ;
return @result ;
} ;
# ====================
sub mt_plus {
my ($p0,$p1,$p2,$p3,undef) = @_ ;
my $mt_plus =
	(defined $p0 and defined $p1 and defined $p2 and defined $p3 and $p0 ne '' and $p1 ne '' and $p2 ne '' and $p3 ne '') ?
		($p0 <= 1) ?
			($p1 == 0) ?
				($p2 == 0) ?
					($p3 == 0) ? 4 : 3 : 2 : 1 : 0 : '' ;
return $mt_plus ;
} ;
# ====================
sub mt_minus {
my ($m0,$m1,$m2,$m3,undef) = @_ ;
my $mt_minus =
	(defined $m0 and defined $m1 and defined $m2 and defined $m3 and $m0 ne '' and $m1 ne '' and $m2 ne '' and $m3 ne '') ?
		($m0 == 0) ?
			($m1 == 0) ?
				($m2 == 0) ?
					($m3 == 0) ? 4 : 3 : 2 : 1 : 0 : '' ;
return $mt_minus ;
} ;
# ====================
sub get_nearest_neighbor_param {
# parameters from Freier et al. PNAS 83:9373 (1986).
my %rna_dh = (
'aa' => -6.6,
'uu' => -6.6,
'au' => -5.7,
'ua' => -8.1,
'ca' => -10.5,
'ug' => -10.5,
'cu' => -7.6,
'ag' => -7.6,
'ga' => -13.3,
'uc' => -13.3,
'gu' => -10.2,
'ac' => -10.2,
'cg' => -8.0,
'gc' => -14.2,
'gg' => -12.2,
'cc' => -12.2
) ;
my %rna_ds = (
'aa' => -18.4,
'uu' => -18.4,
'au' => -15.5,
'ua' => -22.6,
'ca' => -27.8,
'ug' => -27.8,
'cu' => -19.2,
'ag' => -19.2,
'ga' => -35.5,
'uc' => -35.5,
'gu' => -26.2,
'ac' => -26.2,
'cg' => -19.4,
'gc' => -34.9,
'gg' => -29.7,
'cc' => -29.7
) ;
return (\%rna_dh,\%rna_ds) ;
} ;
# ====================
sub seedTm {
if (my $si19 = check_19mer($_[0])){
	my $seed = substr($si19,11,7) ;
	$seed =~ tr/tT/uU/ ;
	(my $tm = tm_RNA($seed)) =~ s/(?<=\..).*// ;  # 小数第2位以下を切り捨て
	return $tm ;
} else {
	return '' ;
}
} ;
# ====================
sub tm_RNA {
my $oligoseq = lc ($_[0] || '') ;
if ($oligoseq =~ /^[augc]+$/){
#- 100nM, [Na+] = 100mM
return ( 1000 * deltaH_RNA($oligoseq) / ( -10.8 + deltaS_RNA($oligoseq) + 1.987 * log(0.0001/4) ) - 273.15 + 16.6 * log(0.1)/log(10) ) ;
#- 定数をゼロにしている
#return ( 1000 * deltaH_RNA($oligoseq) / ( deltaS_RNA($oligoseq) + 1.987 * log(0.0001) ) - 273.15 ) ;
#- 程先生の条件
#return ( 1000 * deltaH_RNA($oligoseq) / ( -10.8 + deltaS_RNA($oligoseq) + 1.987 * log(0.000001/4) ) - 273.15 + 16.6 * log(0.1)/log(10) ) ;
} else {
	return '' ;
}
} ;
# ====================
sub deltaH_RNA {
my $oligoseq = lc ($_[0] || '') ;  # [atgc]+ であることが保証されているのでチェックを省略。
my $sum_dH = 0 ;
foreach (0..length($oligoseq) - 2){
	my $dinucleotide = substr($oligoseq,$_,2) ;
	$sum_dH += $rna_dh{$dinucleotide} ;
}
return $sum_dH ;
} ;
# ====================
sub deltaS_RNA {
my $oligoseq = lc ($_[0] || '') ;  # [atgc]+ であることが保証されているのでチェックを省略。
my $sum_dS = 0 ;
foreach (0..length($oligoseq) - 2){
	my $dinucleotide = substr($oligoseq,$_,2) ;
	$sum_dS += $rna_ds{$dinucleotide} ;
}
return $sum_dS ;
} ;
# ====================
sub complementaryDNA {  # 相補鎖（DNA）を求める
my $seq = flatsequence($_[0]) ;
(my $comp = reverse $seq) =~ tr/GATUCRYMKSWHBVDNgatucrymkswhbvdn/CTAAGYRKMSWDVBHNctaagyrkmswdvbhn/ ;
return $comp ;
} ;
# ====================
sub gc_content {
if (my $si19 = check_19mer($_[0])){
	my $gc_count = ($si19 =~ tr/GCgc/GCgc/) ;
	my $gc_percent = $gc_count/length($si19)*100 ;
	return $gc_percent ;
} else {
	return '' ;
}
} ;
# ====================
sub iub2regexp {
if ($_[0] =~ /^[ATGCUNRYMKSWHBVD]+$/i){
	my $regexp = $_[0] ;
	$regexp =~ s/R/[AG]/gi ;
	$regexp =~ s/Y/[CTU]/gi ;
	$regexp =~ s/M/[AC]/gi ;
	$regexp =~ s/K/[GTU]/gi ;
	$regexp =~ s/S/[CG]/gi ;
	$regexp =~ s/W/[ATU]/gi ;
	$regexp =~ s/H/[ACTU]/gi ;
	$regexp =~ s/B/[CGTU]/gi ;
	$regexp =~ s/V/[ACG]/gi ;
	$regexp =~ s/D/[AGTU]/gi ;
	$regexp =~ s/N/[ACGTU]/gi ;
	return $regexp ;
} else {
	return '' ;
}
} ;
# ====================
sub RNAoligo_guide {
if (defined $_[0] and my $seq = complementaryDNA(flatsequence($_[0]))){
	if (length $seq == 23){
		(my $oligoseq = substr($seq,2,21)) =~ tr/tT/uU/ ;
		return uc $oligoseq ;
	}
}
return '' ;
} ;
# ====================
sub RNAoligo_passenger {
if (defined $_[0] and my $seq = flatsequence($_[0])){
	if (length $seq == 23){
		(my $oligoseq = substr($seq,2,21)) =~ tr/tT/uU/ ;
		return uc $oligoseq ;
	}
}
return '' ;
} ;
# ====================
sub print_siRNA_table {
# siRNAリストを引数に与える。
# @siRNAlistの要素は、ハッシュのリファレンス。

unless (defined $_[0] and ref $_[0] eq 'ARRAY'){
	print_error_html('ERROR : print_siRNA_table() : siRNA list not defined.') ;
	exit ;
}
my @siRNAlist = @{$_[0]} ;

my %option ;
if (defined $_[1] and ref $_[1] eq 'HASH'){
	%option = %{$_[1]} ;
}

my $html =
"<table cellspacing=0 cellpadding=2>

<tr>
	<td class=vw rowspan=3>target<br>position
	<td class=vw rowspan=3>target sequence<br>21nt target + 2nt overhang
	<td class=vw rowspan=3>RNA oligo sequences<br>21nt guide (5&prime;&rarr;3&prime;)<br>21nt passenger (5&prime;&rarr;3&prime;)
	<td class=ow rowspan=3>functional<br>siRNA<br>selection:<br>\n" ;

# 効くsiRNAの選択・結果表示
$option{'uitei'} and $html .= "\t\t<font class=u>U</font>i-Tei<br>\n" ;
$option{'reynolds'} and $html .= "\t\t<font class=r>R</font>eynolds<br>\n" ;
$option{'amarzguioui'} and $html .= "\t\t<font class=a>A</font>marzguioui<br>\n" ;
# seed部分のTmによる選択・結果表示
$option{'seedTm'} and $html .= "\t<td class=pw colspan=3 rowspan=2>seed-duplex<br>stabilty (Tm);\n" ;
# 19mer/3mm検索のprecomputed DB照会・結果表示
($option{'hs'} or $option{'mm'} or $option{'rn'}) and $html .=
	"\t<td class=gw colspan=2 rowspan=2>specificity check:<br>minimum number of<br>mismatches against<br>any off-targets;\n" ;
$option{'hitcount'} and ($option{'hs'} or $option{'mm'} or $option{'rn'}) and $html .=
	"\t<td class=gw colspan=8>specificity check:<br>number of off-target hits with<br>indicated mismatches(strand)\n" ;
# Target rangeによる選択・結果表示
$option{'pos'} and $html .= "\t<td class=vw rowspan=3>target<br>position\n" ;
# Gの連続・Cの連続による選択・結果表示
$option{'consGC'} and $option{'consGCmax'} and $html .= "\t<td class=vw rowspan=3>contiguous<br>G's or C's<br>constraint\n" ;
# Aの連続・Tの連続による選択・結果表示
$option{'consAT'} and $option{'consATmax'} and $html .= "\t<td class=vw rowspan=3>contiguous<br>A's or T's<br>constraint\n" ;
# GC含量による選択・結果表示
$option{'percentGC'} and $html .= "\t<td class=vw rowspan=3>GC<br>content\n" ;
# Custom patternによる選択・結果表示
$option{'custom'} and $option{'customPattern'} and $html .= "\t<td class=vw rowspan=3>custom<br>pattern\n" ;
# Exclude patternによる選択・結果表示
$option{'exclude'} and $option{'excludePattern'} and $html .= "\t<td class=vw rowspan=3>exclude<br>pattern\n" ;
$html .= "\t</tr>\n" ;

$html .= "<tr>\n" ;
# 19mer/3mm検索のprecomputed DB照会・結果表示
$option{'hitcount'} and ($option{'hs'} or $option{'mm'} or $option{'rn'}) and $html .=
	"\t<td class=gw colspan=4>guide (+)<td class=gw colspan=4>passenger (-)\n" ;
$html .= "\t</tr>\n" ;

$html .= "<tr>\n" ;
# seed部分のTmによる選択・結果表示
$option{'seedTm'} and $html .= "\t<td class=pw><td class=pw>guide<td class=pw>passenger\n" ;
# 19mer/3mm検索のprecomputed DB照会・結果表示
($option{'hs'} or $option{'mm'} or $option{'rn'}) and $html .=
	"\t<td class=gw>guide<td class=gw>passenger\n" ;
$option{'hitcount'} and ($option{'hs'} or $option{'mm'} or $option{'rn'}) and $html .=
	"\t<td class=gw>0(+)<td class=gw>1(+)<td class=gw>2(+)<td class=gw>3(+)<td class=gw>0(-)<td class=gw>1(-)<td class=gw>2(-)<td class=gw>3(-)\n" ;

$html .= "</tr>\n" ;

my ($v,$o,$g,$p) = ('w','w','w','w') ;
my ($v_nextline,$o_nextline,$g_nextline,$p_nextline) = ('v','o','g','p') ;

foreach (@siRNAlist){
	#- @siRNAlist：ハッシュのリファレンスを格納したアレイ（2008-12-08変更）
	%$_ = (
		'startpos' => '','endpos' => '','si23' => '','name' => '',
		'efficacy' => '','seed_tm' => '','seed_tm_sense' => '',
		'p0' => '','m0' => '','p1' => '','m1' => '','p2' => '','m2' => '','p3' => '','m3' => '',
		'pos' => '','consGC' => '','consAT' => '','percentGC' => '','custom' => '','exclude' => '',
		'color' => '','mt_plus' => '','mt_minus' => '',
	%$_) ;  # undefを返さないようにデフォルト値として空白文字を代入しておく。
	my ($startpos,$endpos,$si23,$name,$efficacy,$seed_tm,$seed_tm_sense,$p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,
		$pos,$consGC,$consAT,$percentGC,$custom,$exclude,$color,$mt_plus,$mt_minus) =
		@$_{'startpos','endpos','si23','name','efficacy','seed_tm','seed_tm_sense','p0','m0','p1','m1','p2','m2','p3','m3',
		'pos','consGC','consAT','percentGC','custom','exclude','color','mt_plus','mt_minus'} ;
	if ($efficacy){ # HTML化。U/R/Aに色をつける。
		$efficacy =~ s|(?<=.)(?=.)| |g ;
		$efficacy =~ s|U|<font class=u>U</font>| ;
		$efficacy =~ s|R|<font class=r>R</font>| ;
		$efficacy =~ s|A|<font class=a>A</font>| ;
	}
	foreach (($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3)){  # hit数が100を越えるものは>100と表示する。
		if ($_ and $_ > 100){
			$_ = '&gt;100' ;
		}
	}
	my $si19g = substr($si23,1,19) ;  # guide鎖off-target検索用
	my $si19p = substr($si23,3,19) ;  # passenger鎖off-target検索用
	my $rnaseq_gs = RNAoligo_guide($si23) ;  # guide鎖オリゴ配列計算
	my $rnaseq_ps = RNAoligo_passenger($si23) ;  # passenger鎖オリゴ配列計算
	$html .= "<tr><td class=$v>$name
	<td class=$v>" ;
	$html .= "$si23\n" ;
	$html .= "<td class=$v>$rnaseq_gs<br>$rnaseq_ps\n" ;
	$html .= "	<td class=$o>" . $efficacy . "\n" ;
	# ▼ 他のcriteriaをチェック
	# $color_html を定義
	my $color_html =
	(not $color) ? '' :  # 未定義または0の場合 → 白
	($color == 3) ? '<font class=color3>' :  # 濃い青
	($color == 2) ? '<font class=color2>' :  # 青
	($color == 1) ? '<font class=color1>' :  # うす青
		'' ;  # 上記以外の場合 → 白。
	my $color_html_closetag = $color_html ? '</font>' : '' ;
	$option{'seedTm'} and $html .= "\t<td class=$p>$color_html&nbsp;&nbsp;$color_html_closetag\n" ;
	$option{'seedTm'} and $html .= "\t<td class=$p>$seed_tm &deg;C\n" ;
	$option{'seedTm'} and $html .= "\t<td class=$p>$seed_tm_sense &deg;C\n" ;
	$html .=
		$option{'hs'} ? "\t<td class=$g>$mt_plus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19g . '&strand=plus&spe=hs">[detail]</a>' .
			"<td class=$g>$mt_minus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19p . '&strand=minus&spe=hs">[detail]</a>' . "\n" :
		$option{'mm'} ? "\t<td class=$g>$mt_plus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19g . '&strand=plus&spe=mm">[detail]</a>' .
			"<td class=$g>$mt_minus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19p . '&strand=minus&spe=mm">[detail]</a>' . "\n" :
		$option{'rn'} ? "\t<td class=$g>$mt_plus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19g . '&strand=plus&spe=rn">[detail]</a>' .
			"<td class=$g>$mt_minus " .
			'<a target="_blank" href="detail.cgi?seq=' . $si19p . '&strand=minus&spe=rn">[detail]</a>' . "\n" :
			'' ;
	$option{'hitcount'} and $html .= (defined $p0 and not $p0 eq '') ?
			"\t<td class=$g>$p0<td class=$g>$p1<td class=$g>$p2<td class=$g>$p3<td class=$g>$m0<td class=$g>$m1<td class=$g>$m2<td class=$g>$m3\n" :
		$option{'hs'} ? "\t<td class=$g colspan=8><font class=s>not available</font> " . "\n" :
#			'<a target="_blank" href="siDirectHs.cgi?val=' . $si23 . '">[show detail]</a>' . "\n" :
		$option{'mm'} ? "\t<td class=$g colspan=8><font class=s>not available</font> " . "\n" :
#			'<a target="_blank" href="siDirectMm.cgi?val=' . $si23 . '">[show detail]</a>' . "\n" :
		$option{'rn'} ? "\t<td class=$g colspan=8><font class=s>not available</font> " . "\n" :
#			'<a target="_blank" href="siDirectRn.cgi?val=' . $si23 . '">[show detail]</a>' . "\n" :
			'' ;
	$option{'pos'} and $html .= $pos ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;
	$option{'consGC'} and $option{'consGCmax'} and $html .= $consGC ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;  # &#10003; はwin IEで四角になる。
	$option{'consAT'} and $option{'consATmax'} and $html .= $consAT ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;
	$option{'percentGC'} and $html .= $percentGC ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;
	$option{'custom'} and $option{'customPattern'} and $html .= $custom ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;
	$option{'exclude'} and $option{'excludePattern'} and $html .= $exclude ? "\t<td class=$v>ok\n" : "\t<td class=$v>\n" ;
	# ▲ 他のcriteriaをチェック
	$html .= "</tr>\n" ;
	($v,$o,$g,$p,$v_nextline,$o_nextline,$g_nextline,$p_nextline) = ($v_nextline,$o_nextline,$g_nextline,$p_nextline,$v,$o,$g,$p) ;
}
$html .= "\n</table>\n" ;

return $html ;
} ;
# ====================
sub print_siRNA_txt {
# print_siRNA_tableをベースに改変
# siRNAリストを引数に与える。
# @siRNAlistの要素は、ハッシュのリファレンス。

unless (defined $_[0] and ref $_[0] eq 'ARRAY'){
	print_error_html('ERROR : print_siRNA_txt() : siRNA list not defined.') ;
	exit ;
}
my @siRNAlist = @{$_[0]} ;

my %option ;
if (defined $_[1] and ref $_[1] eq 'HASH'){
	%option = %{$_[1]} ;
}

my $txt = 'target position	target sequence	RNA oligo, guide	passenger	functional siRNA selection' ;

# seed部分のTmによる選択・結果表示
$option{'seedTm'} and $txt .= '	seed-duplex stabilty (Tm), guide	passenger' ;
# 19mer/3mm検索のprecomputed DB照会・結果表示
($option{'hs'} or $option{'mm'} or $option{'rn'}) and $txt .= '	min. number of mismatches against off-targets, guide	passenger' ;
$option{'hitcount'} and ($option{'hs'} or $option{'mm'} or $option{'rn'}) and $txt .= '	number of off-target hits, 0(+)	1(+)	2(+)	3(+)	0(-)	1(-)	2(-)	3(-)' ;
# Target rangeによる選択・結果表示
$option{'pos'} and $txt .= '	target position' ;
# Gの連続・Cの連続による選択・結果表示
$option{'consGC'} and $option{'consGCmax'} and $txt .= "	contiguous G's or C's constraint" ;
# Aの連続・Tの連続による選択・結果表示
$option{'consAT'} and $option{'consATmax'} and $txt .= "	contiguous A's or T's constraint" ;
# GC含量による選択・結果表示
$option{'percentGC'} and $txt .= '	GC content' ;
# Custom patternによる選択・結果表示
$option{'custom'} and $option{'customPattern'} and $txt .= '	custom pattern' ;
# Exclude patternによる選択・結果表示
$option{'exclude'} and $option{'excludePattern'} and $txt .= '	exclude pattern' ;

$txt .= "\n" ;

foreach (@siRNAlist){
	#- @siRNAlist：ハッシュのリファレンスを格納したアレイ（2008-12-08変更）
	%$_ = (
		'startpos' => '','endpos' => '','si23' => '','name' => '',
		'efficacy' => '','seed_tm' => '','seed_tm_sense' => '',
		'p0' => '','m0' => '','p1' => '','m1' => '','p2' => '','m2' => '','p3' => '','m3' => '',
		'pos' => '','consGC' => '','consAT' => '','percentGC' => '','custom' => '','exclude' => '',
		'mt_plus' => '','mt_minus' => '',
	%$_) ;  # undefを返さないようにデフォルト値として空白文字を代入しておく。
	my ($startpos,$endpos,$si23,$name,$efficacy,$seed_tm,$seed_tm_sense,$p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,
		$pos,$consGC,$consAT,$percentGC,$custom,$exclude,$mt_plus,$mt_minus) =
		@$_{'startpos','endpos','si23','name','efficacy','seed_tm','seed_tm_sense','p0','m0','p1','m1','p2','m2','p3','m3',
		'pos','consGC','consAT','percentGC','custom','exclude','mt_plus','mt_minus'} ;
	foreach (($p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3)){  # hit数が100を越えるものは>100と表示する。
		if ($_ and $_ > 100){
			$_ = '>100' ;
		}
	}
	my $rnaseq_gs = RNAoligo_guide($si23) ;  # guide鎖オリゴ配列計算
	my $rnaseq_ps = RNAoligo_passenger($si23) ;  # passenger鎖オリゴ配列計算
	$txt .= "$name	$si23	$rnaseq_gs	$rnaseq_ps	$efficacy" ;
	# ▼ 他のcriteriaをチェック
	$option{'seedTm'} and $txt .= "	$seed_tm	$seed_tm_sense" ;
	($option{'hs'} or $option{'mm'} or $option{'rn'}) and $txt .= "	$mt_plus	$mt_minus" ;
	$option{'hitcount'} and ($option{'hs'} or $option{'mm'} or $option{'rn'}) and $txt .= "	$p0	$p1	$p2	$p3	$m0	$m1	$m2	$m3" ;
	$option{'pos'} and $txt .= $pos ? "	1" : "	0" ;
	$option{'consGC'} and $option{'consGCmax'} and $txt .= $consGC ? "	1" : "	0" ;
	$option{'consAT'} and $option{'consATmax'} and $txt .= $consAT ? "	1" : "	0" ;
	$option{'percentGC'} and $txt .= $percentGC ? "	1" : "	0" ;
	$option{'custom'} and $option{'customPattern'} and $txt .= $custom ? "	1" : "	0" ;
	$option{'exclude'} and $option{'excludePattern'} and $txt .= $exclude ? "	1" : "	0" ;
	# ▲ 他のcriteriaをチェック
	$txt .= "\n" ;
}

return $txt ;
} ;
# ====================
sub print_graphical_view {
# usage:
# print_graphical_view($targetseq,$linewidth,$startpos,\@siRNAlist)
# @siRNAlistの要素は、ハッシュのリファレンス。

my $max_targetseq_length = 100000 ;  # 実験用に。10000→100000
my $default_linewidth = 100 ;
my $default_startpos = 1 ;

# ▼ 引数の処理
unless ($_[0]){
	print_error_html('ERROR : print_graphical_view() : target sequence not defined.') ;
	exit ;
}
my $targetseq = $_[0] ;
if (length $targetseq > $max_targetseq_length){  # 塩基配列の長さをチェック
	print_error_html('ERROR : print_graphical_view() : target sequence is too long.') ;
	exit ;
}
if ($targetseq =~ /[^ATGCUNRYMKSWHBVD-]/i){  # 塩基構成文字以外が含まれていないかチェック
	print_error_html('ERROR : print_graphical_view() : target sequence contains invalid characters.') ;
	exit ;
}

my $linewidth = $default_linewidth ;  # 100
if (defined $_[1] and $_[1] =~ /^\d+$/){
	$linewidth = $_[1] ;
}

my $startpos = $default_startpos ;  # 1
if (defined $_[2] and $_[2] =~ /^\d+$/){
	$startpos = $_[2] ;
}

my @siRNAlist ;
if (defined $_[3] and ref $_[3] eq 'ARRAY'){
	@siRNAlist = @{$_[3]} ;
}
# ▲ 引数の処理

my $html = '' ;

if ($targetseq =~ s/^(.{0,$linewidth})//){
	my $targetseq_thisline = $1 ;
	my $targetlength_thisline = length $targetseq_thisline ;

# ▼ 座標の表示
	$html .= "<tr>\n" ;
	my $pos = $startpos + 9 ;
	while ($pos < $startpos + $targetlength_thisline){
		$html .= "<td colspan=10 class=pos>$pos\n" ;
		$pos += 10 ;
	}
	$html .= "<td>&nbsp;\n" ;
	$html .= "</tr>\n\n" ;
# ▲ 座標の表示

# ▼ 塩基配列の表示
	$html .= "<tr>\n" ;
	(my $targetseq_html = $targetseq_thisline) =~ s/(?=.)/<td>/g ;
	$html .= $targetseq_html ;
	my $endpos = $startpos + $targetlength_thisline - 1 ;
	my $space_html = '' ;
	if ($targetlength_thisline < $linewidth){  # 最終行の処理
		my $space = $linewidth - $targetlength_thisline ;
		$space_html = " colspan=$space" ;
	}
	$html .= "\n<td class=posrange$space_html>&nbsp;&nbsp;$startpos-$endpos\n</tr>\n\n" ;
# ▲ 塩基配列の表示

# ▼ 現在の行に表示するsiRNAの選別
	my ($show_siRNAlist_ref,$next_siRNAlist_ref) = select_siRNA_thisline(\@siRNAlist,$targetlength_thisline) ;
	my @show_siRNAlist = @{$show_siRNAlist_ref} ;
	my @next_siRNAlist = @{$next_siRNAlist_ref} ;
# ▲ 現在の行に表示するsiRNAの選別

# ▼ siRNAの表示
	$html .= print_siRNA_thisline(@show_siRNAlist) ;
# ▲ siRNAの表示

# ▼ 空白行
	$html .= "\n<tr><td>&nbsp;</tr><!-- __________________________________________________ -->\n\n" ;
# ▲ 空白行

# ▼ $targetseq が残っている場合、再帰
	my $html2 = '' ;
	if ($targetseq){
		$html2 = print_graphical_view($targetseq,$linewidth,$startpos+$targetlength_thisline,\@next_siRNAlist) ;
	}
# ▲ $targetseq が残っている場合、再帰

	return $html . $html2 ;
}
} ;
# ====================
sub select_siRNA_thisline {  # 現在の行に表示するsiRNAの選別
# usage:
# my ($show_siRNAlist_ref,$next_siRNAlist_ref) = select_siRNA_thisline(\@siRNAlist,$targetlength_thisline) ;

unless (defined $_[0] and ref $_[0] eq 'ARRAY'){
	print_error_html('ERROR : select_siRNA_thisline() : siRNA list not defined.') ;
	exit ;
}
if (not defined $_[1]){
	print_error_html('ERROR : select_siRNA_thisline() : targetlength_thisline not defined.') ;
	exit ;
}
my @siRNAlist = @{$_[0]} ;
my $targetlength_thisline = $_[1] ;
my @show_siRNAlist ;
my @next_siRNAlist ;
foreach (@siRNAlist){
	%$_ = (
		'startpos' => '','endpos' => '','si23' => '','name' => '',
		'efficacy' => '','seed_tm' => '','seed_tm_sense' => '',
		'p0' => '','m0' => '','p1' => '','m1' => '','p2' => '','m2' => '','p3' => '','m3' => '',
		'color' => '',
	%$_) ;  # undefを返さないようにデフォルト値として空白文字を代入しておく。
	my ($startpos,$endpos,$si23,$name,$efficacy,$seed_tm,$seed_tm_sense,$p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,$color) =
		@$_{'startpos','endpos','si23','name','efficacy','seed_tm','seed_tm_sense','p0','m0','p1','m1','p2','m2','p3','m3','color'} ;
	if (defined $startpos and $startpos =~ /^\d+$/ and 
		defined $endpos and $endpos =~ /^\d+$/ and 
		$startpos < $endpos
	){
		defined $si23 or $si23 = '' ;
		defined $name or $name = '*' ;
		defined $efficacy or $efficacy = '' ;
		defined $seed_tm or $seed_tm = '' ;
		defined $seed_tm_sense or $seed_tm_sense = '' ;
		defined $color or $color = '' ;
		# ▼ 現在の行に表示するsiRNAのリスト
		if ($startpos <= $targetlength_thisline){
			push @show_siRNAlist, {
				'startpos' => $startpos,
				'endpos' => $endpos,
				'si23' => $si23,
				'name' => $name,
				'efficacy' => $efficacy,
				'seed_tm' => $seed_tm,
				'seed_tm_sense' => $seed_tm_sense,
				'p0' => $p0,
				'm0' => $m0,
				'p1' => $p1,
				'm1' => $m1,
				'p2' => $p2,
				'm2' => $m2,
				'p3' => $p3,
				'm3' => $m3,
				'color' => $color,
			} ;
		}
		# ▲ 現在の行に表示するsiRNAのリスト
		# ▼ 次以降の行に表示するsiRNAのリスト
		if ($endpos > $targetlength_thisline) {
			my $new_startpos = $startpos - $targetlength_thisline ;
			if ($new_startpos < 1){
				$new_startpos = 1 ;
				$name = '&nbsp;' ;
				$efficacy = '' ;
			}
			my $new_endpos = $endpos - $targetlength_thisline ;
			push @next_siRNAlist, {
				'startpos' => $new_startpos,
				'endpos' => $new_endpos,
				'si23' => $si23,
				'name' => $name,
				'efficacy' => $efficacy,
				'seed_tm' => $seed_tm,
				'seed_tm_sense' => $seed_tm_sense,
				'p0' => $p0,
				'm0' => $m0,
				'p1' => $p1,
				'm1' => $m1,
				'p2' => $p2,
				'm2' => $m2,
				'p3' => $p3,
				'm3' => $m3,
				'color' => $color,
			} ;
		}
		# ▲ 次以降の行に表示するsiRNAのリスト
	}
}
return (\@show_siRNAlist,\@next_siRNAlist) ;
} ;
# ====================
sub print_siRNA_thisline {
# usage:
# $html .= print_siRNA_thisline(@show_siRNAlist) ;
my @show_siRNAlist = @_ ;
my $html = '' ;
foreach (@show_siRNAlist){
	%$_ = (
		'startpos' => '','endpos' => '','si23' => '','name' => '',
		'efficacy' => '','seed_tm' => '','seed_tm_sense' => '',
		'p0' => '','m0' => '','p1' => '','m1' => '','p2' => '','m2' => '','p3' => '','m3' => '',
		'color' => '',
	%$_) ;  # undefを返さないようにデフォルト値として空白文字を代入しておく。
	my ($startpos,$endpos,$si23,$name,$efficacy,$seed_tm,$seed_tm_sense,$p0,$m0,$p1,$m1,$p2,$m2,$p3,$m3,$color) =
		@$_{'startpos','endpos','si23','name','efficacy','seed_tm','seed_tm_sense','p0','m0','p1','m1','p2','m2','p3','m3','color'} ;
	# $space_html を定義
	my $space_html = '' ;
	if ($startpos > 1){
		my $space = $startpos - 1 ;
		$space_html = "<td colspan=$space>" ;
	}
	# $si_length を定義
	my $si_length = $endpos - $startpos + 1 ;
	# $color_html を定義
	my $color_html =
	(not $color) ? '' :  # 未定義または0の場合 → 白
	($color == 3) ? ' bgcolor="#6666FF"' :  # 濃い青
	($color == 2) ? ' bgcolor="#9999FF"' :  # 青
	($color == 1) ? ' bgcolor="#DDDDFF"' :  # うす青
		'' ;  # 上記以外の場合 → 白。
	# $name_html を定義
	my $name_html = $name || '*' ;
	# $efficacy を整形
	if ($efficacy){ # HTML化。U/R/Aに色をつける。
		$efficacy =~ s|(?<=.)(?=.)| |g ;
		$efficacy =~ s|U|<font class=u>U</font>| ;
		$efficacy =~ s|R|<font class=r>R</font>| ;
		$efficacy =~ s|A|<font class=a>A</font>| ;
	}
	# 書き出し
	$html .= "<tr>$space_html<td class=si colspan=$si_length$color_html>\n" ;
	$html .= #$option{'hs'} ? '	<a target="_blank" href="siDirectHs.cgi?val=' . $si23 . '">' . "$name_html</a> $efficacy\n" :
				#$option{'mm'} ? '	<a target="_blank" href="siDirectMm.cgi?val=' . $si23 . '">' . "$name_html</a> $efficacy\n" :
				#$option{'rn'} ? '	<a target="_blank" href="siDirectRn.cgi?val=' . $si23 . '">' . "$name_html</a> $efficacy\n" :
										"	$name_html</a> $efficacy\n" ;
	# ▼ siRNA positionを示す罫線（new!）
	# 注意：siRNAlistがpositionでsortされている必要あり。
	my @others_startpos_list ;
	foreach (@show_siRNAlist){
		my $others_startpos = $$_{'startpos'} ;
		if ($others_startpos > $endpos + 1){
			push @others_startpos_list, $others_startpos ;
		}
	}
	@others_startpos_list = sort {$a <=> $b} @others_startpos_list ;
	my $current_pos = $endpos + 1 ;
	foreach (@others_startpos_list){
		my $space = $_ - $current_pos ;
		$current_pos = $_ ;
		$html .= "\t<td class=guide colspan=$space>\n" ;
	}
	# ▲ siRNA positionを示す罫線（new!）
	$html .= "\t</tr>\n" ;
}
return $html ;
} ;
# ====================
sub print_result_html {
my $siRNA_table_html = $_[0] || '' ;
my $graphical_view_html = $_[1] || '' ;
my $siRNA_table_txt = $_[2] || '' ;
my $timestamp = $_[3] || '' ;
my $query_summary_html = $_[4] || '' ;
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
	textarea { font-size:8pt }
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
	.pos {
		color:green;
		text-align:right;
	}
	.posrange {
		color:green;
		white-space:nowrap;
	}
	.si {
		border-top:1px solid gray;
		border-left:1px solid gray;
		border-bottom:2px solid black;
		border-right:2px solid black;
	}
	.guide { border-right:1px solid silver }
	.vw { background-color:#E6E6FA; border-left:1px solid white; white-space:nowrap; text-align:center }
	.ow { background-color:#FFF0E0; border-left:1px solid white; white-space:nowrap; text-align:center }
	.gw { background-color:#E6FAF2; border-left:1px solid white; white-space:nowrap; text-align:center }
	.pw { background-color:#FFEDED; border-left:1px solid white; white-space:nowrap; text-align:center }
	.v { background-color:#E6E6FA; border-left:1px solid #CCCCCC; white-space:nowrap; text-align:center }
	.o { background-color:#FFF0E0; border-left:1px solid #CCCCCC; white-space:nowrap; text-align:center }
	.g { background-color:#E6FAF2; border-left:1px solid #CCCCCC; white-space:nowrap; text-align:center }
	.p { background-color:#FFEDED; border-left:1px solid #CCCCCC; white-space:nowrap; text-align:center }
	.w { border-left:1px solid #CCCCCC; white-space:nowrap; text-align:center }
	.u { background-color:#FF99FF }
	.r { background-color:#66FFFF }
	.a { background-color:#99FF00 }
	.b { background-color:#6688FF }
	.c { background-color:#DDDDFF }
	.s { background-color:silver }
	.color1 { background-color:#DDDDFF; border:1px solid black }
	.color2 { background-color:#9999FF; border:1px solid black }
	.color3 { background-color:#6666FF; border:1px solid black }
-->
</style>
</head>

<body>

<a href="/"><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt">
<font size=5>siDirect </font><font size=4>version 2.0 </font>result page.
<a target="_blank" href="doc/"><img src="help.gif" alt="Help" border=0></a>
</div>

<hr><!-- __________________________________________________ -->

<p><font color=gray>' . $timestamp . ',  siDirect v.2.0</font></p>

<h3>Query</h3>

' . $query_summary_html .'
<h3>Effective siRNA candidates</h3>

<!--
<table style="border:1px solid gray; padding:5px; margin:1em">
<tr>
	<td style="border:1px dotted white"><font class=b><i>start-end</i></font>
	<td>Functional, off-target reduced siRNA
	</tr>
<tr>
	<td style="border:1px dotted gray"><i>start-end</i>
	<td>Functional siRNA
	</tr>
<tr>
	<td style="border:1px dotted white"><font class=s><i>start-end</i></font>
	<td>Functional siRNA (number of off-target hits not available)
	</tr>
</table>
-->

' . $siRNA_table_html . '

<h3>Graphical view of effective siRNA candidates</h3>

<table style="border:1px solid gray; padding:5px; margin:1em">
<tr>
	<td class=si bgcolor="#6666FF"><i>start-end</i>
	<td>Functional, off-target reduced siRNA (seed duplex Tm < 10 &deg;C)
	</tr>
<tr>
	<td class=si bgcolor="#9999FF"><i>start-end</i>
	<td>Functional, off-target reduced siRNA (seed duplex Tm < 15 &deg;C)
	</tr>
<tr>
	<td class=si bgcolor="#DDDDFF"><i>start-end</i>
	<td>Functional, off-target reduced siRNA (seed duplex Tm < 21.5 &deg;C)
	</tr>
<tr>
	<td class=si><i>start-end</i>
	<td>Functional siRNA
	</tr>
<!--
<tr>
	<td style="text-align:center"><font class=u>U</font>
	<td>Satisfies Ui-Tei algorithm
	</tr>
<tr>
	<td style="text-align:center"><font class=r>R</font>
	<td>Satisfies Reynolds algorithm
	</tr>
<tr>
	<td style="text-align:center"><font class=a>A</font>
	<td>Satisfies Amarzguioui algorithm
	</tr>
-->
</table>

<table cellspacing=0 cellpadding=0>

' . $graphical_view_html . '
</table>

<h3>Tab-delimited list <font size=2>(for data export)</font></h3>

<textarea readonly rows=15 cols=90>
[siDirect v.2.0 | ' . $timestamp . ']
' . $siRNA_table_txt . '</textarea>

</body>
</html>
' ;
exit ;
} ;
# ====================
sub print_error_html {
my $error_mesg = $_[0] || 'ERROR : print_error_html() : no mesg.' ;
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
	* { font-family:verdana,arial,helvetica,sans-serif }
	p,table,div { font-size:8pt }
/* hyperlink */
	a:link,a:visited {
		text-decoration:none;
		color:#004080;
	}
	a:hover {
		text-decoration:none;
		color:red;
	}
-->
</style>
</head>

<body>

<a href="/"><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt"><font size=5>siDirect </font><font size=4>version 2.0</font></div>

<hr><!-- __________________________________________________ -->

<p><font color=red>' . $error_mesg . '</font></p>

</body>
</html>
' ;
exit ;
} ;
# ====================
sub print_redirect_html {
# usage:
# print_redirect_html(msg,URL) ;
my $error_mesg = $_[0] || 'ERROR : print_redirect_html() : no mesg.' ;
unless ($_[1]){
	print_error_html('ERROR : print_redirect_html() : URL not defined.') ;
	exit ;
}
my $redirect_url = $_[1] ;
print 'Content-type: text/html; charset=utf-8

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html lang=en>

<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="Refresh" content="2;URL=' . $redirect_url . '">
<meta http-equiv="Content-Script-Type" content="text/javascript">
<meta http-equiv="Content-Style-Type" content="text/css">
<meta name="author" content="Yuki Naito et al.">
<title>siDirect</title>
<script type="text/javascript">
<!--
	function refresh() { location.href = "' . $redirect_url . '"; }
-->
</script>
<style type="text/css">
<!--
/* font */
	* { font-family:verdana,arial,helvetica,sans-serif }
	p,table,div { font-size:8pt }
/* hyperlink */
	a:link,a:visited {
		text-decoration:none;
		color:#004080;
	}
	a:hover {
		text-decoration:none;
		color:red;
	}
-->
</style>
</head>

<body onload="setTimeout(\'refresh()\', 2000);">

<a href="/"><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt"><font size=5>siDirect </font><font size=4>version 2.0</font></div>

<hr><!-- __________________________________________________ -->

<p><font color=red>' . $error_mesg . '</font></p>

<p>Return to: <a href="' . $redirect_url . '">
' . $redirect_url . '</a></p>

</body>
</html>
' ;
exit ;
} ;
# ====================
sub print_no_siRNA_html {
my $timestamp = $_[0] || '' ;
my $query_summary_html = $_[1] || '' ;
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
	textarea { font-size:8pt }
/* hyperlink */
	a:link,a:visited {
		text-decoration:none;
		color:#004080;
	}
	a:hover {
		text-decoration:none;
		color:red;
	}
-->
</style>
</head>

<body>

<a href="/"><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt">
<font size=5>siDirect </font><font size=4>version 2.0 </font>result page.
<a target="_blank" href="doc/"><img src="help.gif" alt="Help" border=0></a>
</div>

<hr><!-- __________________________________________________ -->

<p><font color=gray>' . $timestamp . ',  siDirect v.2.0</font></p>

<h3>Query</h3>

' . $query_summary_html .'
<h3><font color="#800000">No effective siRNA candidates were found.</font></h3>

<p style="font-size:9pt">
Try relaxing various parameters, including functional siRNA selection algorithms and seed-duplex stability parameters, etc.<br>
Giving the following parameters will help you find out why no acceptable siRNAs were found.
</p>

<img src="doc/nosiRNA.gif">

<h3>Tab-delimited list <font size=2>(for data export)</font></h3>

<textarea readonly rows=15 cols=90>
[siDirect v.2.0 | ' . $timestamp . ']
No effective siRNA candidates were found.
</textarea>

</body>
</html>
' ;
exit ;
} ;
# ====================
