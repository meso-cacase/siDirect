#!/usr/bin/perl

# -----------------------------
# siDirect web server, TOP page
# -----------------------------
#
# 2008-09-01 Yuki Naito (@meso_cacase)
# 2023-07-10 Yuki Naito (@meso_cacase) Version 2.1

use warnings ;
use strict ;
use LWP::Simple ;
use FindBin ;                     # perl 5.26以降はカレントディレクトリはlibrary pathに含まれないので追加する
use lib $FindBin::Bin ;
eval 'use DBlist ; 1' or          # データベースの正式名
	warn('ERROR : cannot load DBlist') ;

my $googleanalyticscode =
'<!-- Google Analytics tracking code -->
<script type="text/javascript">
var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
document.write(unescape("%3Cscript src=\'" + gaJsHost + "google-analytics.com/ga.js\' type=\'text/javascript\'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
try {
var pageTracker = _gat._getTracker("UA-9921668-1");
pageTracker._trackPageview();
} catch(err) {}
</script>

' ;

my %query = get_query_parameters() ;  # CGIが受け取ったデータの処理
(my $accession = $query{'accession'} || '') =~ s/(^\s+|\s+$)//g ;  # 前後の空白文字を除去
$query{'debug'} and my $debug = 1 ;

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
print_top_html('',$sampleseq,'',$debug) ;
	exit ;
}

my $fasta = get(
	'https://www.ncbi.nlm.nih.gov/sviewer/?report=fasta&retmode=text&val=' . $accession
) || "ERROR : cannot retrieve sequence." ;
$fasta =~ s/\s+\z// ;

my $speflag = ($fasta =~ /Homo sapiens/) ? 'Hs' :
					($fasta =~ /Mus musculus/) ? 'Mm' :
					($fasta =~ /Rattus norvegicus/) ? 'Rn' : '' ;
print_top_html($accession,$fasta,$speflag,$debug) ;

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
sub print_top_html {
my $accession = $_[0] || '' ;
my $fasta = $_[1] || '' ;
my $speflag = $_[2] || '' ;
my $debug = $_[3] || '' ;
if ($accession and not $fasta){
	$fasta = 'not found.' ;
}

# 2021-04-12 DB追加の度にソースを書き換えないで済むように変更
my $dbcnt = 0 ;	# 表示順用変数
my %dbinfo = () ;
my %species2defaultIndex = () ;
my $defaultIndexScript = "\tif (0) {}\n" ;
my $dbconf = $DBlist::dbconfig ;  # データベース名と正式名称のリスト
foreach (split /\n/, $DBlist::dbconfig){
	chomp ;
	my ($db, $desc, $species) = split /\t/ ;
	$dbinfo{$db} = {
		'no'   => ++$dbcnt,
		'desc' => $desc
	} ;
	unless (defined($species2defaultIndex{$species})){
		$species2defaultIndex{$species} = $dbcnt ;  # 配列取得時のデフォルトDB選択用
		# 2021-04-12 配列取得時のDB選択ロジックを動的に変更
		$defaultIndexScript .= <<EOS ;
	else if (res.indexOf("${species}") != -1) {
		document.getElementById("nrdbSpe").selectedIndex = $species2defaultIndex{$species} ;
		document.getElementById("hideLessSpecific").checked = true ;
	}
EOS
	}
}

my $nrdb = "" ;
my $fSelected = 0 ;
foreach my $dbname (sort {$dbinfo{$a}->{'no'} <=> $dbinfo{$b}->{'no'}} keys %dbinfo){
	if ($fSelected == 0 && length($speflag) == 0){
		# speflagに値がセットされていなければ最初のDBを選択状態にしておく
		$nrdb .= "<option value=${dbname} selected>${dbinfo{$dbname}->{'desc'}}</option>\n" ;
		$fSelected = 1 ;
	} elsif ($fSelected == 0 && ${dbname} =~ /^${speflag}/i){
		# speflagに値がセットされている場合、先頭がspeflagと一致したらデフォルト選択にしておく
		$nrdb .= "<option value=${dbname} selected>${dbinfo{$dbname}->{'desc'}}</option>\n" ;
		$fSelected = 1 ;
	} else {
		$nrdb .= "<option value=${dbname}>${dbinfo{$dbname}->{'desc'}}</option>\n" ;
	}
}
$nrdb .= <<EOF ;
</select><br>
&nbsp;&nbsp;&nbsp;&nbsp;
<input type=checkbox checked id=hideLessSpecific name=hidenonspe value=1>
	Hide less-specific siRNAs
EOF

# デバッグ用ページ。デバッグフラグを次のページに引き継ぐ。
my $debug_html1 = $debug ?
'<input type=hidden name=debug value=1>' : '' ;
# デバッグ用ページ。全siRNAを表示。
my $debug_html2 = $debug ?
'<p>
Debug:<br>
<input type=checkbox name=ALL value=1>
	Show all siRNAs
</p>
' : '' ;

print
'Content-type: text/html; charset=utf-8

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html lang=en>

<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="Content-Script-Type" content="text/javascript">
<meta http-equiv="Content-Style-Type" content="text/css">
<meta name="author" content="Yuki Naito et al.">
<title>siDirect</title>
<script type="text/javascript">
<!--
function retrieveSequence(){
	acc = document.getElementById("acc").value ;
	requestFile(\'\',\'GET\',\'./retrieveFASTA.cgi?accession=\'+acc,true) ;
}
function requestFile(data, method, fileName, async){
	var httpoj = createHttpRequest() ;
	httpoj.onreadystatechange = function(){
		if (httpoj.readyState==4){ 
			on_loaded(httpoj) ;
		}
	}
	httpoj.open(method, fileName, async) ;
	httpoj.send(data) ;
}
function createHttpRequest(){
	//for Win IE
	if(window.ActiveXObject){
		try {
			//for MSXML2
			return new ActiveXObject("Msxml2.XMLHTTP") ;
		} catch (e) {
			try {
				//for MSXML
				return new ActiveXObject("Microsoft.XMLHTTP") ;
			} catch (e2) {
				return null ;
			}
		}
	} else if(window.XMLHttpRequest){
		//other agents
		return new XMLHttpRequest() ;
	} else {
		return null ;
	}
}
function on_loaded(oj){
	res  = oj.responseText ;
	document.getElementById("useq").value = res ;
' . $defaultIndexScript . '
}

img1 = new Image(); img1.src = "doc/venn_all.svg" ;
img2 = new Image(); img2.src = "doc/venn_uitei.svg" ;
img3 = new Image(); img3.src = "doc/venn_reynolds.svg" ;
img4 = new Image(); img4.src = "doc/venn_amarzguioui.svg" ;
img5 = new Image(); img5.src = "doc/venn_combined.svg" ;
img6 = new Image(); img6.src = "doc/venn_combined2.svg" ;
img7 = new Image(); img7.src = "doc/venn_UorRorA.svg" ;
function showOption(optName){
	optMenu = document.getElementById(optName).style ;
	if (optMenu.display == \'none\') optMenu.display = "block" ; else optMenu.display = "none" ;
}
function more2hide(id){
	obj = document.getElementById(id) ;
	if (obj.innerHTML == \'[more options...]\') obj.innerHTML = "[hide options]" ; else obj.innerHTML = "[more options...]" ;
}
function chk_combined_rule(obj){
	if (obj.checked == true){
			document.getElementsByName("uitei"      )[0].checked  = true ;
			document.getElementsByName("uitei"      )[0].disabled = true ;
			document.getElementsByName("reynolds"   )[0].checked  = true ;
			document.getElementsByName("reynolds"   )[0].disabled = true ;
			document.getElementsByName("amarzguioui")[0].checked  = true ;
			document.getElementsByName("amarzguioui")[0].disabled = true ;
			document.getElementsByName("UorRorA"    )[0].checked  = false ;
			document.getElementsByName("UorRA"      )[0].checked  = false ;
			document.getElementsByName("URA"        )[0].checked  = false ;
			                                         obj.checked  = true ;
	} else {
			document.getElementsByName("uitei"      )[0].disabled = false ;
			document.getElementsByName("reynolds"   )[0].disabled = false ;
			document.getElementsByName("amarzguioui")[0].disabled = false ;
	}
}
//-->
</script>
<style type="text/css">
<!--
/* font */
	p,table,div,h1,h2,h3 { font-family:verdana,arial,helvetica,sans-serif }
	p,table,div,textarea { font-size:8pt }
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

<a href="."><img src="ocean.jpg" height=80 width="100%" border=0></a>
<div style="font-size:8pt">
<a style="color:black" href="."><font size=5>siDirect </font><font size=4>version 2.1<small>&beta;</small> </font></a>highly effective, target specific siRNA online design site.
<a target="_blank" href="doc/"><img src="help.png" alt="Help" width=36 height=15 border=0></a>
</div>

<hr><!-- __________________________________________________ -->

<script type="text/javascript">
<!--
var ua = navigator.userAgent.toLowerCase() ;
if (ua.indexOf("msie") != -1 && ua.indexOf("mac") != -1){
	document.write(\'<form name=getSequence method=GET action=".">\') ;
	document.write(\'' . $debug_html1 . '\') ;
	document.write(\'<p>Enter an accession number and retrieve sequence:<br>\') ;
	document.write(\'<input type=text name=accession size=30 maxlength=50 value="' . $accession . '">\') ;
	document.write(\'<input type=submit value="retrieve sequence">\') ;
	document.write(\'</p>\') ;
	document.write(\'</form>\') ;
} else {
	document.write(\'<p>Enter an accession number and retrieve sequence:<br>\') ;
	document.write(\'<input type=text id=acc name=accession size=30 maxlength=50 value="' . $accession . '" onkeypress="if(event.keyCode==13){retrieveSequence()}">\') ;
	document.write(\'<input type=button value="retrieve sequence" onclick="retrieveSequence()">\') ;
	document.write(\'<\/p>\') ;
}
//-->
</script>

<noscript>
<form name=getSequence method=GET action=".">
' . $debug_html1 . '
<p>Enter an accession number and retrieve sequence:<br>
<input type=text name=accession size=30 maxlength=50 value="' . $accession . '">
<input type=submit value="retrieve sequence">
</p>

</form>
</noscript>

<form name=siDirect method=POST action="design.cgi">
' . $debug_html1 . '
<p>or Paste in a nucleotide sequence:<br>
<textarea id=useq name=yourSeq rows=15 cols=80>
' . $fasta . '
</textarea></p>

<p><input type=submit value="design siRNA"></p>

<h3>Options:</h3>

<div style="background-color:#FFF0E0; width:90%; margin-bottom:1em">
<table cellpadding=0 cellspacing=0 width="100%">
<tr>
	<td>

<p><b>Functional siRNA selection</b> algorithm by: &nbsp;
<a target="_blank" href="doc/algorithms.html">
<img src="qicon.png" alt="?" width=15 height=15 border=0></a>
</p>
<p>
<input type=checkbox checked name=uitei value=1>
	<span onmouseover="venn.src=img2.src" onmouseout="venn.src=img1.src">
	Ui-Tei <i>et al.</i>, <i>Nucleic Acids Res</i> <b>32</b>, 936-948 (2004)
	<a target="_blank" href="https://doi.org/10.1093/nar/gkh247">Link</a></span><br>
<script type="text/javascript">
<!--
document.write(\'<a id=more1 href="javascript:showOption(\\\'functional_options\\\');showOption(\\\'venn\\\');more2hide(\\\'more1\\\')">[more options...]<\/a><br>\') ;
document.write(\'<div id=functional_options style="display:none">\') ;
//-->
</script>
<input type=checkbox name=reynolds value=1>
	<span onmouseover="venn.src=img3.src" onmouseout="venn.src=img1.src">
	Reynolds <i>et al.</i>, <i>Nat Biotechnol</i> <b>22</b>, 326-330 (2004)
	<a target="_blank" href="https://doi.org/10.1038/nbt936">Link</a></span><br>
<input type=checkbox name=amarzguioui value=1>
	<span onmouseover="venn.src=img4.src" onmouseout="venn.src=img1.src">
	Amarzguioui <i>et al.</i>, <i>BBRC</i> <b>316</b>, 1050-1058 (2004)
	<a target="_blank" href="https://doi.org/10.1016/j.bbrc.2004.02.157">Link</a></span>
</p>
<p>
use combined rule:<br>
<input type=checkbox name=UorRorA value=1 onclick="chk_combined_rule(this)">
	<span onmouseover="venn.src=img7.src" onmouseout="venn.src=img1.src">
	Ui-Tei + Reynolds + Amarzguioui</span><br>
<input type=checkbox name=UorRA value=1 onclick="chk_combined_rule(this)">
	<span onmouseover="venn.src=img5.src" onmouseout="venn.src=img1.src">
	Ui-Tei + Reynolds &times; Amarzguioui</span><br>
<input type=checkbox name=URA value=1 onclick="chk_combined_rule(this)">
	<span onmouseover="venn.src=img6.src" onmouseout="venn.src=img1.src">
	Ui-Tei &times; Reynolds &times; Amarzguioui</span>
<script type="text/javascript">
<!--
document.write(\'<\/div>\') ;
//-->
</script>
</p>
' . $debug_html2 . '
	</td>
	<td>
<script type="text/javascript">
<!--
document.write(\'<div id=venn style="display:none">\') ;
//-->
</script>
		<a target="_blank" href="doc/algorithms.html">
		<img src="doc/venn_all.svg" name=venn width=200 border=0></a>
<script type="text/javascript">
<!--
document.write(\'<\/div>\') ;
//-->
</script>
	</td>
</tr>
</table>
</div>

<div style="background-color:#FFEDED; width:90%; margin-bottom:1em">
<table cellpadding=0 cellspacing=0 width="100%">
<tr>
	<td>

<p><b>Minimization of seed-dependent off-target effects</b> &nbsp;
<a href="doc/seedTm.html" onclick="window.open(\'doc/seedTm.html\',\'\',\'width=530,height=700\') ; return false ;">
<img src="qicon.png" alt="?" width=15 height=15 border=0></a>
</p>
<p><input type=checkbox checked name=seedTm value=1>
	<b>Seed-duplex stability</b>: Max Tm
	<input type=text name=seedTmMax size=7 maxlength=12 value="21.5">&deg;C
	<br>&nbsp;&nbsp;&nbsp;&nbsp;
	(for reducing seed-dependent off-target effect)
	<br>&nbsp;&nbsp;&nbsp;&nbsp;
	Ui-Tei <i>et al.</i>, <i>Nucleic Acids Res</i> <b>36</b>, 7100-7109 (2008)
	<a target="_blank" href="https://doi.org/10.1093/nar/gkn902">Link</a><br>
</p>

	</td>
	<td>
		<a href="doc/seedTm.html" onclick="window.open(\'doc/seedTm.html\',\'\',\'width=530,height=700\') ; return false ;">
		<img src="doc/seedTm.svg" width=180 border=0></a>
	</td>
</tr>
</table>
</div>

<div style="background-color:#E6FAF2; width:90%; margin-bottom:1em">
<p><b>Specificity check</b>:
<select id=nrdbSpe name=spe>
	<option value=none>none</option>
' . $nrdb . '<br>
<script type="text/javascript">
<!--
document.write(\'<a id=more2 href="javascript:showOption(\\\'specificity_options\\\');more2hide(\\\'more2\\\')">[more options...]<\/a><br>\') ;
document.write(\'<div id=specificity_options style="display:none">\') ;
//-->
</script>
&nbsp;&nbsp;&nbsp;&nbsp;
<input type=checkbox checked name=hitcount value=1>
	Show number of off-target hits within three mismatches
<script type="text/javascript"> 
<!--
document.write(\'<\/div>\') ;
//-->
</script>
</p>
</div>

<div style="background-color:#E6E6FA; width:90%; margin-bottom:1em">
<p><b>Other options</b></p>
<p><input type=checkbox name=pos value=1>
	<b>Target range</b>: from
	<input type=text name=posStart size=7 maxlength=12 value=start> to
	<input type=text name=posEnd size=7 maxlength=12 value=end>
</p>

<p><input type=checkbox name=consGC value=1>
	<b>Avoid contiguous G\'s or C\'s</b>
	<select name=consGCmax>
	<option value=4 selected>4
	<option value=5>5
	<option value=6>6
	<option value=7>7
	</select>
	nt or more (for chemically synthesized siRNA)
	<br>
	<input type=checkbox name=consAT value=1>
	<b>Avoid contiguous A\'s or T\'s</b>
	<select name=consATmax>
	<option value=4 selected>4
	<option value=5>5
	<option value=6>6
	<option value=7>7
	</select>
	nt or more (for shRNA vectors with pol III promoter)
</p>

<p><input type=checkbox name=percentGC value=1>
	<b>GC content</b>: from
	<input type=text name=percentGCMin size=5 maxlength=5 value=0>% to
	<input type=text name=percentGCMax size=5 maxlength=5 value=100>%
</p>

<p><input type=checkbox name=custom value=1>
	<b>Custom pattern</b>:
	<input type=text name=customPattern size=40 maxlength=23 value="NNGNNNNNNNNNNNNNNNNNNNN">
	<a href="doc/pattern.html" onclick="window.open(\'doc/pattern.html\',\'\',\'width=530,height=700\') ; return false ;">
	<img src="qicon.png" alt="?" width=15 height=15 border=0></a><br>
	<input type=checkbox name=exclude value=1>
	<b>Exclude pattern</b>:
	<input type=text name=excludePattern size=40 maxlength=23 value="">
</p>

<hr><!-- __________________________________________________ -->

<p>&nbsp;&nbsp;&nbsp;&nbsp;<input type=checkbox name=hide value=1>
	Only show siRNAs that match all checked criteria 
</p>
</div>

</form>

<hr><!-- __________________________________________________ -->

<p>
<font color=gray>siDirect v2.1&beta; | Last modified on Aug 25, 2023.<br>
Previsous version: <a href="v2.0/">siDirect v2.0</a></font>
</p>

' . $googleanalyticscode . '</body>
</html>
' ;
exit ;
} ;
# ====================
