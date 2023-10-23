package DBlist ;

# オフターゲット検索対象データベース一覧
#
# フォーマット：
# DB	fullname	学名
#
# 生物種ごとのデフォルトDBはより前方に記載されているものが採用される

$dbconfig =
<<'--EOS--' ;
hs_refseq220	Human (Homo sapiens) transcript, RefSeq release 220 (Sep, 2023)	Homo sapiens
hs_refseq215	Human (Homo sapiens) transcript, RefSeq release 215 (Nov, 2022)	Homo sapiens
hs_refseq210	Human (Homo sapiens) transcript, RefSeq release 210 (Jan, 2022)	Homo sapiens
hs_refseq205	Human (Homo sapiens) transcript, RefSeq release 205 (Mar, 2021)	Homo sapiens
hs_refseq200	Human (Homo sapiens) transcript, RefSeq release 200 (May, 2020)	Homo sapiens
mm_refseq220	Mouse (Mus musculus) transcript, RefSeq release 220 (Sep, 2023)	Mus musculus
mm_refseq215	Mouse (Mus musculus) transcript, RefSeq release 215 (Nov, 2022)	Mus musculus
mm_refseq210	Mouse (Mus musculus) transcript, RefSeq release 210 (Jan, 2022)	Mus musculus
mm_refseq205	Mouse (Mus musculus) transcript, RefSeq release 205 (Mar, 2021)	Mus musculus
mm_refseq200	Mouse (Mus musculus) transcript, RefSeq release 200 (May, 2020)	Mus musculus
rn_refseq220	Rat (Rattus norvegicus) transcript, RefSeq release 220 (Sep, 2023)	Rattus norvegicus
rn_refseq215	Rat (Rattus norvegicus) transcript, RefSeq release 215 (Nov, 2022)	Rattus norvegicus
rn_refseq210	Rat (Rattus norvegicus) transcript, RefSeq release 210 (Jan, 2022)	Rattus norvegicus
rn_refseq205	Rat (Rattus norvegicus) transcript, RefSeq release 205 (Mar, 2021)	Rattus norvegicus
rn_refseq200	Rat (Rattus norvegicus) transcript, RefSeq release 200 (May, 2020)	Rattus norvegicus
--EOS--

return 1 ;
