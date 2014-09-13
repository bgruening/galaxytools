#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: prinseq-graphs
#   Date: 2012-12-22
#   Version: 0.6 graphs
#
#   Usage:
#      prinseq-graphs [options]
#
#      Try 'prinseq-graphs -h' for more information.
#
#    Purpose: PRINSEQ will help you to preprocess your genomic or metagenomic
#             sequence data in FASTA or FASTQ format. The graphs version allows
#             users of the lite version to generate graphs similar to the web
#             version.
#
#    Bugs: Please use http://sourceforge.net/tracker/?group_id=315449
#
#===============================================================================

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile); #for output files
use Fcntl qw(:flock SEEK_END); #for log file
use Cwd;
use JSON;
use Cairo;
use Statistics::PCA;
use MIME::Base64;
use File::Basename;
use Data::Dumper; ###

$| = 1; # Do not buffer output

my $PI = 4 * atan2(1, 1);
my $LOG62 = log(62);
my $DINUCODDS_VIR = [
                        [qw(1.086940308	0.98976932	1.034167044	0.880024041	1.070421277	0.990687084	0.890945575	1.069957074	0.92465631	0.803973303)],
                        [qw(1.101064857	0.986812783	1.038299155	0.896162618	1.081652847	0.976365237	0.867445186	1.06727283	0.94688543	0.768007295)],
                        [qw(1.071548411	0.912204166	1.196914981	0.80628184	1.294201511	1.148517794	0.269295791	1.033948026	0.895951033	0.623192149)],
                        [qw(1.090253719	0.907428629	1.203991784	0.786359294	1.281499107	1.145421568	0.235974709	1.033437274	0.899580091	0.631699771)],
                        [qw(1.075864745	1.003413074	1.01872902	0.897841689	0.980373171	1.05854979	0.934262259	1.052477953	0.88145851	0.889239724)],
                        [qw(1.101890467	1.030028291	1.019912674	0.84191395	1.0015174	1.069546264	0.900151602	0.996269395	0.889195343	0.904039022)],
                        [qw(1.152417359	0.855028574	0.91164793	1.017415486	1.114163672	1.128353311	0.846355573	0.916745489	1.206820475	0.811014651)],
                        [qw(1.142454218	0.8635465	0.923406967	1.026242747	1.134445058	1.131747833	0.79793368	0.920767641	1.179468556	0.799770057)],
                        [qw(1.124462747	0.873556143	0.945627041	1.013755408	1.159866153	1.096259526	0.757315047	0.972924919	1.105562567	0.772731886)],
                        [qw(1.143826972	0.866968779	0.995740249	0.945859278	1.109590621	1.089305083	0.76048874	0.971561388	1.157101408	0.792923027)],
                        [qw(1.131900141	0.82776996	0.996204924	0.999433455	1.024692372	1.071176333	0.921026216	1.088936699	1.054010776	0.773498892)],
                        [qw(1.042180476	0.930180412	1.019242897	0.98909997	1.006666828	1.046708539	0.959492164	1.011183418	1.055168776	0.937433818)],
                        [qw(1.086515695	0.985345815	0.930914307	0.969581792	1.043010232	1.087463712	0.939482285	0.990551965	0.954752469	0.893972874)],
                        [qw(1.096657826	0.950117614	0.936195529	0.965619788	1.114975275	1.077011195	0.843153131	0.989128406	1.043790912	0.840634731)],
                        [qw(1.158030995	0.935307365	0.874812261	1.056236525	1.117171274	0.937484692	1.057442372	0.970079538	1.174848738	0.725071711)],
                        [qw(1.15591506	0.93000227	0.883538923	1.0567652	1.095730954	0.944489906	1.074229471	0.983993745	1.156051409	0.726688465)],
                        [qw(1.205726473	0.924439339	1.049457756	0.805718412	0.975472778	1.07581991	0.726992211	1.075025787	0.8704929	0.726672843)],
                        [qw(1.188544681	0.95239611	1.049066985	0.790031334	1.038632598	1.056749787	0.665197397	1.057566244	0.862429061	0.708982398)],
                        [qw(1.063631482	0.925593715	1.014869316	0.944904401	1.119690731	1.325971834	0.273781451	0.943347677	1.06438014	0.920825904)],
                        [qw(1.077560287	0.911888545	1.044147857	0.927758054	1.058535939	1.296838544	0.421514996	0.945722451	1.128317986	0.926419928)],
                        [qw(1.163753415	0.989905668	0.893599328	0.955641844	1.176047687	0.941559156	0.950641089	0.959741692	1.100815282	0.72491925)],
                        [qw(1.139253929	0.946297517	0.922096125	1.024801537	1.205206793	0.968818717	0.915801342	0.971626058	1.107569276	0.627623404)]
                     ];
my $DINUCODDS_MIC = [
                        [qw(1.13127323	0.853587195	0.911041047	1.104520778	1.065586428	1.021434164	0.999734139	1.063684014	1.078035184	0.733596552)],
                        [qw(1.173267344	0.840539337	0.919534602	1.068050141	1.062394214	1.051999071	0.96770576	1.035511729	1.095600433	0.72328141)],
                        [qw(1.172939786	0.84567902	0.911836259	1.106288994	1.05351787	1.026143368	1.002308358	1.066319771	1.094918797	0.710733535)],
                        [qw(1.073527689	0.850290918	0.978455025	1.080882178	1.111174765	1.010754115	0.895668707	1.072980666	1.079304608	0.754057386)],
                        [qw(1.08807747	0.837444678	0.95824965	1.097310298	1.118897971	1.030863881	0.886827263	1.072349394	1.07406322	0.733440096)],
                        [qw(1.071685485	0.861055813	0.966566865	1.090268118	1.112945761	1.012538936	0.909535491	1.063745603	1.071156598	0.755770377)],
                        [qw(1.142698587	0.867936867	1.000612099	0.977934257	1.111801746	1.018318601	0.788556794	0.987763594	1.184649653	0.784776176)],
                        [qw(1.134560074	0.876651844	0.998190253	0.995723123	1.128448077	1.014172324	0.781776188	0.971020602	1.182411449	0.786449476)],
                        [qw(1.180029632	0.787899325	1.01316945	0.932268406	1.077837263	1.211699678	0.612128817	1.033036699	1.157314398	0.74940288)],
                        [qw(1.160925546	0.788308899	1.003702496	0.965371236	1.076051693	1.188304271	0.641536444	1.070331188	1.124067192	0.740126813)],
                        [qw(1.173873006	0.790118011	1.014718833	0.937979878	1.07453725	1.207167373	0.622279064	1.046150047	1.145627707	0.742212886)],
                        [qw(1.128383111	0.870541389	0.987269741	0.98353238	1.115643879	1.040107028	0.774505865	1.010896432	1.164757274	0.775254395)],
                        [qw(1.15297511	0.853883985	0.956393231	1.000027661	1.139915472	1.01355294	0.838843622	1.015553125	1.216219741	0.70447264)],
                        [qw(1.148264236	0.852123859	0.974568293	0.985455546	1.13192373	1.015879393	0.828987111	1.016820786	1.216647853	0.71634006)],
                        [qw(1.12933788	0.831777975	1.005434367	0.991081409	1.126146895	1.07421504	0.69343913	1.054032466	1.14809591	0.728541157)],
                        [qw(1.124157235	0.828112691	1.022348424	0.983822386	1.143028487	1.081830005	0.672594435	1.05685982	1.149537403	0.684432106)],
                        [qw(1.128029586	0.841853305	1.00983936	0.967179139	1.122524003	1.094555807	0.659238308	1.061578854	1.1243601	0.740148171)],
                        [qw(1.093521636	0.855071052	0.929160818	1.203773691	1.178257185	0.881341255	1.078305505	1.051988532	1.169143967	0.555057308)],
                        [qw(1.073737278	0.877396537	0.968017446	1.124155374	1.166244435	0.909044208	0.999147578	1.071098934	1.120156138	0.607444953)],
                        [qw(1.092150184	0.863407008	0.927040387	1.185387013	1.171670826	0.882276859	1.083058605	1.048379554	1.168635365	0.580337997)]
                    ];
my $DATA_VIR = [
                    [2,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [3,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [42,2,'Human (nasal)',[127/255, 127/255, 255/255,1]],
                    [43,2,'Human (nasal)',[127/255, 127/255, 255/255,1]],
                    [45,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [49,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [52,3,'Human (sputum)',[127/255, 127/255, 255/255,1]],
                    [54,3,'Human (sputum)',[127/255, 127/255, 255/255,1]],
                    [55,4,'Human (sputum, CF)',[127/255, 127/255, 255/255,1]],
                    [57,4,'Human (sputum, CF)',[127/255, 127/255, 255/255,1]],
                    [88,5,'Freshwater (Hot spring)',[127/255, 127/255, 255/255,1]],
                    [89,5,'Freshwater (Hot spring)',[127/255, 127/255, 255/255,1]],
                    [98,6,'Freshwater (Antartic lake)',[127/255, 127/255, 255/255,1]],
                    [99,6,'Freshwater (Antartic lake)',[127/255, 127/255, 255/255,1]],
                    [100,7,'Freshwater (reclaimed)',[127/255, 127/255, 255/255,1]],
                    [102,7,'Freshwater (reclaimed)',[127/255, 127/255, 255/255,1]],
                    [153,8,'Mouse (brain tissue)',[127/255, 127/255, 255/255,1]],
                    [154,8,'Mouse (brain tissue)',[127/255, 127/255, 255/255,1]],
                    [202,9,'Fish (gut)',[127/255, 127/255, 255/255,1]],
                    [206,9,'Fish (gut)',[127/255, 127/255, 255/255,1]],
                    [209,10,'Mosquito',[127/255, 127/255, 255/255,1]],
                    [211,10,'Mosquito',[127/255, 127/255, 255/255,1]],
		    ['U',0,'User input',[255/255, 127/255, 127/255,1]]
                ];
my $DATA_MIC = [
                    [17,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [20,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [22,1,'Human (fecal)',[127/255, 127/255, 255/255,1]],
                    [63,2,'Mouse (fecal)',[127/255, 127/255, 255/255,1]],
                    [65,2,'Mouse (fecal)',[127/255, 127/255, 255/255,1]],
                    [68,2,'Mouse (fecal)',[127/255, 127/255, 255/255,1]],
                    [93,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [95,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [109,4,'Marine (open ocean)',[127/255, 127/255, 255/255,1]],
                    [110,4,'Marine (open ocean)',[127/255, 127/255, 255/255,1]],
                    [111,4,'Marine (open ocean)',[127/255, 127/255, 255/255,1]],
                    [120,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [124,5,'Marine (estuary)',[127/255, 127/255, 255/255,1]],
                    [125,5,'Marine (estuary)',[127/255, 127/255, 255/255,1]],
                    [134,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [146,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [148,3,'Marine (coastal)',[127/255, 127/255, 255/255,1]],
                    [201,6,'Fish (gut)',[127/255, 127/255, 255/255,1]],
                    [203,7,'Fish (slime)',[127/255, 127/255, 255/255,1]],
                    [205,6,'Fish (gut)',[127/255, 127/255, 255/255,1]],
		    ['U',0,'User input',[255/255, 127/255, 127/255,1]]
                ];
my $BASE64_BASES = {A => 'iVBORw0KGgoAAAANSUhEUgAAAEkAAABJCAYAAABxcwvcAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAAAzVJREFUeNrsnMFxo0AQRWe7fJcyMBnYGawyMIe9a0JQJtbefDPOAB33
JmdgZyBlsIpgl9lCLkwJA/N7uhu0XTXlkstI8Oh+agbG355+/XDC8VaNu8htf1ZjI73DJPx59wCg
EN4phDQkNAsWGqCkIeUM7zFrSL7OBDS+VyObMyQrZWsSUlZnACfw5dwgcZ/5BZPfTEHyEwCvColL
2O24q/uuWUDKJ1TGKpCCsB8Sn4Dl1CGlbvxEBD51SCIlR4lL4VYAUnKB08SzSCSbUkFKLWxRgdMM
sii5wK1BOlksuRSQVoCwA9wjIPDVVCAhWVTWw1SZc0MK8lxHblvUP7fA569TCJyMZFET0qEa75ay
iRtSrDwDlLfG663CPohAQoRdtF4jXrrlFjgZKbU2lN/VeLFSclyQlkAzt6s95BiziVXgXJByFz/7
WH7x+6OFbOKCFCvL0wUffeUqFYFzQELu7/eVFAKJTeCkmEVDIARXvWqXHAoJEXbwzZ4BZJ/AM21I
iLCLESV50swmMlxqzZ6pnCqkDBD2a0dvlErguRYkiSw6x16zZyKlDy4FwDbjARE4AYBihf1Se0YS
EnRSaSJZpNozxUAKaRv7QNYR/KZSEXgMpI1CFjUhifdMMZBypUzgAB0lcIoAFDv72J6ijY0tuL1P
DckrZ5GrQSM90yYlpMxh9/cfq/GHaSBPq4xeVUBCWWQt/kMaEKNWFQyFJPVAlmRsuCF5N7/wnJCW
TvaBLKkYLHC60iwadWzEWbtzFXgfpNUMhT06CeiKS23wMVKPsNdXAKlX4HTlWTToWG8SQdoxXK3H
zA7E3r0JAr/vmqXogoSu3w87vFeA9AwK3I8pN+Rr/6gAKAQ669m5qoA6hJ0r7mxsoE/Hda4qoA6i
CzDttaJI0TMRc6mFKdqDIqS9w2YtLy4LowTC1o4tdzYR83VaaQASu8Dpwh/ERuzta+441H0am8Cp
1TwuJp5FSQROTB32yRgk9Om4TwI/Q8oc9g9XCmcv2LKJmIRtERL6LfexqoAYSo3r9nUKgb+D7+HP
kFBhW8wi1p6JHL4KujQMCRX4v1UFARJyu2infBky5KIXPYn+rwADAOL8qKxS08x7AAAAAElFTkSu
QmCC',
		    C => 'iVBORw0KGgoAAAANSUhEUgAAAEkAAABJCAYAAABxcwvcAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAAA7BJREFUeNrsnM1xqzAQxxUNDfBKwCWQCt7g+7vgEkgJ5pRDTnYJpgRz
eXeYVBBKCCU8SvAzM6sZxuMPaXclQaydYYKTGPBv/7tagdYvp9NJTO3Px6dwZPl5S2A/hdf3rD9v
1eT1nvuC/r7/vvr7SLizDGAUEzgmNr5nN3mt9ksAWNu6cNuQYoCyhX0bpmANoK4K9tlMWrrw0euH
8/YPPkTsQKkxnIv9nNKSZ79BQb5sy3kNkjnnfMMFzsFiUHNDVZVk9FyDTMguBowvGDS8QTpejDpz
tARAZT4gNRr1zZyswYCSrk84Azuahp58MkAqoR9NkjkG0m7BgG5V76yQcgtD/B6mFqvz9nJlW8Pf
uacdha6zI0P6B6YLbGH6UGv+b3tRbnCNpgdwDpuSOEr9cU61AXXUBOX9YlJWolOVS4MwyxnUs2L6
cAr2G1MhzAKJKu8K1DMw55UKYFHVlFMhYe//KKuZPH7v+CXxGCyQsNZbBjTNUzURUoyFlFEmhhAK
g3BjVDUVWEg5MV90DgvEy3vgppZi66ScGAKurTJMDxXAvXuPPMLGqUYy7T1A6mBLHxSlRg6MMPLT
hOTLWnBuNVELKS9GD5I2ttDzCalkSOJaiTsmKKkVP8wks4qE4xHNKyRKhd0HSCHcyCPb4LDC9g4p
DqFmL9yGZ4EUkrbhBDeYBSWJoKQAKViAFCAFSLOERKl1kqCkoKSgJFMl9QGSPUijpQHSE6rppypJ
tU5Y7Qig3IL1vZ5ydNJ403BcdzSuZBt71Rp4ncxJSbFHSNmN36melxMAK6iQhgWrSWf9wu6KylBL
byiQCo+hliIcqlTmFFLmaZSjOKfCQFIrNLDmuqUrIULqsHO3muhVl+UAxSl3F3lIDQlSHhMZ9XAQ
w9tKqOlAUs2/lBA4OAgz6jlIkDjUlFsEpTqOqGsXeiokqppUfmqYQy+BY1Lz3sPPJg0O1DPkDXSL
5xV1fjEAanVKHZM7kxtG72ObCjN4L9eAoLUQ36SVqwNFcdQ/GWzTUL6V+7aTn5zhqh0dpl/DUYLE
ueZm6lshhHDbEd4Lg8WnmAcBG7H8dZFGqQMDSfWa9QsG1NmGpOS6XiAoVC+vJMb164JCr8TWe9SH
kwOAqmcO6I1SEEvGON/MEI5KC5QWL9bH3KOaVjNSVQXXQ15XLi14TrW0+1r03kIKYGtrlRYvdM0h
dUPlvMI5WQeTyIFXW/Cqeu5VMIPpheUuTZdfobifjDTTXvxYcz5YXsBxtrD+vwADADoA0kx0ZQr1
AAAAAElFTkSuQmCC',
		    G => 'iVBORw0KGgoAAAANSUhEUgAAAEkAAABJCAYAAABxcwvcAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAAA6lJREFUeNrsnN2RqjAUgANjA9wS3CefsQQsAUtgS8AStAQsQUuQEuR5
nzYlSAkumTnZy/UKJpyTEANn5syqs0L8cn4hIbjf70zI19eaWZS40aT1Pm80Uvhe1ei59b6Ez8hl
tbr+vl5YhpLCa8xx4h54RqCZhCQsI9OwEkYEr2700OgRXqMlNARn3+gN/kbMrrTPXzS6dA2SHFzO
3BBhyd8wrtEhJTAYV+A8Sg7ji8eCJGbpQmHWhkWM7wrJwxqkCODk7L3kpDvmBWJW3sF6+qyfQRY0
YknvDqgNKjUByRdAUgqVYK4DKQJ/9gWQ/E0FJSQl6gNExIVdo9tGgw5dw/8cDJw/fhXIA8UGN8cW
ZA9ybPVaQ4vEjHDSapgI/qzBDRXjEBUgAeaj0U8EIAl5Dcepidwux7hbQTRTG3ApTmyRa6LOP+3q
M0OFLybIk1fwQ0pmRjhMQEVgTdkQSHsCQBti6+mzVE5gTVqQMmS6l4BqZkckKGymi3UhYQa8tQio
7Xo7gisaSpASZHrdWXCxvrqLI61JqcFNkW52HLmSPmrG0yOA5ezfGw2dxaSI8t9s+GXXjcFMppOp
bj21WgWhoHMyX90tSRCAuAOAZEws4XecdS6LPJOFik9qmq0rsqE6UEic1VyCxExBWiJcrRoh5Y8C
CeNqJfNUKCFVU4GEaUP4DGm2JDQkb63oEVKEyGz1lCCxGZJaMemKiKL2PpJeuiDNme0NLck7SNFU
INUzJLOQ2AzptSxnSLO7kaTyyGdQVJC8drmQsJOPpwJpDt4KkDCXYBPisYmbCgFSuSl3qxHuFk3B
krDWlE0FEiZ4p1OBdEZmuHgKkDjSmrIpQMJaU2Yg0zkJCXtPfz8FSDUSVOwTqL4rk9gtCvnI2Y6s
6e6DRLEg6zRSfBLnvNqAJOST4BwXyxZVMOLtZq8gcUazMOtkIUaJrHozUYKo3C2hWm6cgwtQu5/c
qV2Y6h1VINUMv4C8nfUuoBnyOALOHSzU6GWaQOOBLntmZue2XDLMe4rYpHWVwcbu8XK1uv4uTNXZ
zb1j/z+thkJS1xtj3Tu4W+bxYq22JWEgyZ1APoPaPhbSQ9YCSFC+rbYVE//xLC4OXTAhQR08ASTi
7bqr1AkJDr59YziiUP7zarIplt6cu8zUcTjKu8Gp1idxsCjXg/qB/d1yrzxO6pVuJcyQS6VCBWEh
GNpiBYYfoSiLz/0IYM6gg/rO9qbAwOJzJmVrgd0l3pdEGFXGbUP6EWAA2LwDwtC8jpAAAAAASUVO
RK5CYII=',
		    T => 'iVBORw0KGgoAAAANSUhEUgAAAEkAAABJCAYAAABxcwvcAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAAALRJREFUeNrs09ENQDAUQFHEXlhAYgJWMJnEBLqBUWxQFkCC/si5yftq
mzYnaR5jzM4KXXu++J9CNc311YYi022QIEGCBAkSJEiCBCll5c16k+DO4Zj+4dnxmPXj92xvkZYE
SPWLs2uiN/lukCBBggQJkiBBggQJEiRIggQJEiRIkCBBEiRIkCBBggRJkCBBggQJEiRICCBBggQJ
EiRIggQJEiRIkCAJEiRIkCBBggRJ1+0CDAAzsw5U48snWgAAAABJRU5ErkJggg==',
		    N => 'iVBORw0KGgoAAAANSUhEUgAAAEkAAABJCAYAAABxcwvcAAAACXBIWXMAAAsTAAALEwEAmpwYAAAK
T2lDQ1BQaG90b3Nob3AgSUNDIHByb2ZpbGUAAHjanVNnVFPpFj333vRCS4iAlEtvUhUIIFJCi4AU
kSYqIQkQSoghodkVUcERRUUEG8igiAOOjoCMFVEsDIoK2AfkIaKOg6OIisr74Xuja9a89+bN/rXX
Pues852zzwfACAyWSDNRNYAMqUIeEeCDx8TG4eQuQIEKJHAAEAizZCFz/SMBAPh+PDwrIsAHvgAB
eNMLCADATZvAMByH/w/qQplcAYCEAcB0kThLCIAUAEB6jkKmAEBGAYCdmCZTAKAEAGDLY2LjAFAt
AGAnf+bTAICd+Jl7AQBblCEVAaCRACATZYhEAGg7AKzPVopFAFgwABRmS8Q5ANgtADBJV2ZIALC3
AMDOEAuyAAgMADBRiIUpAAR7AGDIIyN4AISZABRG8lc88SuuEOcqAAB4mbI8uSQ5RYFbCC1xB1dX
Lh4ozkkXKxQ2YQJhmkAuwnmZGTKBNA/g88wAAKCRFRHgg/P9eM4Ors7ONo62Dl8t6r8G/yJiYuP+
5c+rcEAAAOF0ftH+LC+zGoA7BoBt/qIl7gRoXgugdfeLZrIPQLUAoOnaV/Nw+H48PEWhkLnZ2eXk
5NhKxEJbYcpXff5nwl/AV/1s+X48/Pf14L7iJIEyXYFHBPjgwsz0TKUcz5IJhGLc5o9H/LcL//wd
0yLESWK5WCoU41EScY5EmozzMqUiiUKSKcUl0v9k4t8s+wM+3zUAsGo+AXuRLahdYwP2SycQWHTA
4vcAAPK7b8HUKAgDgGiD4c93/+8//UegJQCAZkmScQAAXkQkLlTKsz/HCAAARKCBKrBBG/TBGCzA
BhzBBdzBC/xgNoRCJMTCQhBCCmSAHHJgKayCQiiGzbAdKmAv1EAdNMBRaIaTcA4uwlW4Dj1wD/ph
CJ7BKLyBCQRByAgTYSHaiAFiilgjjggXmYX4IcFIBBKLJCDJiBRRIkuRNUgxUopUIFVIHfI9cgI5
h1xGupE7yAAygvyGvEcxlIGyUT3UDLVDuag3GoRGogvQZHQxmo8WoJvQcrQaPYw2oefQq2gP2o8+
Q8cwwOgYBzPEbDAuxsNCsTgsCZNjy7EirAyrxhqwVqwDu4n1Y8+xdwQSgUXACTYEd0IgYR5BSFhM
WE7YSKggHCQ0EdoJNwkDhFHCJyKTqEu0JroR+cQYYjIxh1hILCPWEo8TLxB7iEPENyQSiUMyJ7mQ
AkmxpFTSEtJG0m5SI+ksqZs0SBojk8naZGuyBzmULCAryIXkneTD5DPkG+Qh8lsKnWJAcaT4U+Io
UspqShnlEOU05QZlmDJBVaOaUt2ooVQRNY9aQq2htlKvUYeoEzR1mjnNgxZJS6WtopXTGmgXaPdp
r+h0uhHdlR5Ol9BX0svpR+iX6AP0dwwNhhWDx4hnKBmbGAcYZxl3GK+YTKYZ04sZx1QwNzHrmOeZ
D5lvVVgqtip8FZHKCpVKlSaVGyovVKmqpqreqgtV81XLVI+pXlN9rkZVM1PjqQnUlqtVqp1Q61Mb
U2epO6iHqmeob1Q/pH5Z/YkGWcNMw09DpFGgsV/jvMYgC2MZs3gsIWsNq4Z1gTXEJrHN2Xx2KruY
/R27iz2qqaE5QzNKM1ezUvOUZj8H45hx+Jx0TgnnKKeX836K3hTvKeIpG6Y0TLkxZVxrqpaXllir
SKtRq0frvTau7aedpr1Fu1n7gQ5Bx0onXCdHZ4/OBZ3nU9lT3acKpxZNPTr1ri6qa6UbobtEd79u
p+6Ynr5egJ5Mb6feeb3n+hx9L/1U/W36p/VHDFgGswwkBtsMzhg8xTVxbzwdL8fb8VFDXcNAQ6Vh
lWGX4YSRudE8o9VGjUYPjGnGXOMk423GbcajJgYmISZLTepN7ppSTbmmKaY7TDtMx83MzaLN1pk1
mz0x1zLnm+eb15vft2BaeFostqi2uGVJsuRaplnutrxuhVo5WaVYVVpds0atna0l1rutu6cRp7lO
k06rntZnw7Dxtsm2qbcZsOXYBtuutm22fWFnYhdnt8Wuw+6TvZN9un2N/T0HDYfZDqsdWh1+c7Ry
FDpWOt6azpzuP33F9JbpL2dYzxDP2DPjthPLKcRpnVOb00dnF2e5c4PziIuJS4LLLpc+Lpsbxt3I
veRKdPVxXeF60vWdm7Obwu2o26/uNu5p7ofcn8w0nymeWTNz0MPIQ+BR5dE/C5+VMGvfrH5PQ0+B
Z7XnIy9jL5FXrdewt6V3qvdh7xc+9j5yn+M+4zw33jLeWV/MN8C3yLfLT8Nvnl+F30N/I/9k/3r/
0QCngCUBZwOJgUGBWwL7+Hp8Ib+OPzrbZfay2e1BjKC5QRVBj4KtguXBrSFoyOyQrSH355jOkc5p
DoVQfujW0Adh5mGLw34MJ4WHhVeGP45wiFga0TGXNXfR3ENz30T6RJZE3ptnMU85ry1KNSo+qi5q
PNo3ujS6P8YuZlnM1VidWElsSxw5LiquNm5svt/87fOH4p3iC+N7F5gvyF1weaHOwvSFpxapLhIs
OpZATIhOOJTwQRAqqBaMJfITdyWOCnnCHcJnIi/RNtGI2ENcKh5O8kgqTXqS7JG8NXkkxTOlLOW5
hCepkLxMDUzdmzqeFpp2IG0yPTq9MYOSkZBxQqohTZO2Z+pn5mZ2y6xlhbL+xW6Lty8elQfJa7OQ
rAVZLQq2QqboVFoo1yoHsmdlV2a/zYnKOZarnivN7cyzytuQN5zvn//tEsIS4ZK2pYZLVy0dWOa9
rGo5sjxxedsK4xUFK4ZWBqw8uIq2Km3VT6vtV5eufr0mek1rgV7ByoLBtQFr6wtVCuWFfevc1+1d
T1gvWd+1YfqGnRs+FYmKrhTbF5cVf9go3HjlG4dvyr+Z3JS0qavEuWTPZtJm6ebeLZ5bDpaql+aX
Dm4N2dq0Dd9WtO319kXbL5fNKNu7g7ZDuaO/PLi8ZafJzs07P1SkVPRU+lQ27tLdtWHX+G7R7ht7
vPY07NXbW7z3/T7JvttVAVVN1WbVZftJ+7P3P66Jqun4lvttXa1ObXHtxwPSA/0HIw6217nU1R3S
PVRSj9Yr60cOxx++/p3vdy0NNg1VjZzG4iNwRHnk6fcJ3/ceDTradox7rOEH0x92HWcdL2pCmvKa
RptTmvtbYlu6T8w+0dbq3nr8R9sfD5w0PFl5SvNUyWna6YLTk2fyz4ydlZ19fi753GDborZ752PO
32oPb++6EHTh0kX/i+c7vDvOXPK4dPKy2+UTV7hXmq86X23qdOo8/pPTT8e7nLuarrlca7nuer21
e2b36RueN87d9L158Rb/1tWeOT3dvfN6b/fF9/XfFt1+cif9zsu72Xcn7q28T7xf9EDtQdlD3YfV
P1v+3Njv3H9qwHeg89HcR/cGhYPP/pH1jw9DBY+Zj8uGDYbrnjg+OTniP3L96fynQ89kzyaeF/6i
/suuFxYvfvjV69fO0ZjRoZfyl5O/bXyl/erA6xmv28bCxh6+yXgzMV70VvvtwXfcdx3vo98PT+R8
IH8o/2j5sfVT0Kf7kxmTk/8EA5jz/GMzLdsAAAAgY0hSTQAAeiUAAICDAAD5/wAAgOkAAHUwAADq
YAAAOpgAABdvkl/FRgAAAh1JREFUeNrsmk1xwzAQRr8RgYRBwqBhkDJoGbQMagZ1GbgMVAYNA5dB
wsBm4CBwL9Wx0Uwk7593Z3z0SHmRn3fXi3me8d8FoAUw33kdQB/9PXu9xWCeZ4QFN9zBSCwJ6Qig
cUj5aAFsHdLt2Fh47ALBGi8AHh2ScYlTQXrQLPFAuJZaiVNC2gCIDikfTxolHhjWjA4pH7s/Pzmk
TDQA9g7JUCYeGNdWI/HAvH50SEYkHgTs4V26xIOQfUSHlI8jgGeHlI9OagEsCdIOQtspQdh+REo8
CPzjokNSKPGlIJ0qnKatdUgdgJ/CArhdw+NW+qZ6A888ASmkM4DPCifSvLhbANdCib9ahzRV+JHs
mThFCvCtXeJUeVLpaWKVOBWkAcCH1kycMuPuAIwF97PNE1BCqiHxlkPi1LVbX1iysHyK4ihwm8Lc
iXwojAPSUOE0dNYhJbdctEics5/UVAC9tQ6pB/BVKPFoHVINiZPME3BDmirUZdE6pPSmKimAF58n
kPIhoKlw/946pDPKupiLZuKSPim1FSR+sA6pRgG8sQ4JKO9iYg2QAAGNfw2QBpR3Mc1DSrnT6JCW
l7h5SKkAPjmk5QvgVUAaIGAeQDqklImPDkl47qQFUo+yLuYqILFKXBOkCUzTJZogpUz84pAESlwj
pDPKZzHNQ0q509Uh5SXeOKR8RBB1MTVDIpO4dkgDCLqY2iGl3Gl0SMwS/x0AsYSfWCRqIfIAAAAA
SUVORK5CYII='
		    };
my $MMCHART_B2 = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAAGCAYAAAACEPQxAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAAABdJREFUeNpiYGBg+M/w//9/BmwEQIABANxBD/HRDNRSAAAAAElFTkSu
QmCC';
my $FREQCHART_L = 'iVBORw0KGgoAAAANSUhEUgAAAC8AAABvCAIAAADzHQ6XAAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJ
bWFnZVJlYWR5ccllPAAABWNJREFUeNrkm91V7DoMhTOsaQBKgBKgBCjhUgKUAI/wBiUwJUAJUAKU
ACVACed+K/tcYxz5J5nYk3uOHmYFxokVWdrakj17v1rJ6+trdkw3y0x3d3dd172/v5vfPj8/d/8J
Iytq488U0+bw8PDq6oqLh4cHhn1+fta1jXQytWGB+ErLhB5cPz4+xp6z11WWr68vPvf39/WJnT4+
PmKD17W10dwoYX57c3Oji9vbW2xTXRvpgYVknoQ2fFZfKSkhC6ETFzE7tdDm+PiY6V9eXrh+enri
8/T0NDq6RkxxwZ/EczneSJNuIVgsbVa62rmsVqsfeHN5eanV3aU4W8n5cTqWNobx9ST0G2AbPZzD
c+HcsJk2ht8ACYTiSy8Y7J9eUmG5hQRYbMQUy+SMhDasnUyVyL1NV4rgRCcs1E6btBdDUPiqnTYX
FxeA5q8dyZxYXMJ5W2jDmrqcHHMsTO7GxFDD0OauF/cI1i4bR4yRPylTmixTMOHIvPlMafMjM1xf
X/tgAOqcnZ2JSsaEMcwkVoVab29vJhlVqOozwUQ7n0sHb8br8lqZmsMzPnZivuEYgWc6MKXJ2mev
AeAe9pKwTZrzOkHdo14Yr9ceYvEP24gQDUODmbK2cRY1bYOX8BDZj1kwtpmSQy8+7oXbPnvhTiwc
u9n3UKcut6uK80UVXfb1Qm0wTGDzrNP49ojFlJZGz9EKmEFu4w1Dr3rh/pI0aeJNwIt5lCtfhsZb
LhYvixev/XC9v78fwpdPM2rLtzbALgqJW+2YpRMLXO+cUayDKrWxLWwsFnZl0aWdbYAEcjiuE7Q2
wLdsJprZb9QDMyXRGKtYT8VytdkFql6H7/eCThR1akVJSp5lkqyt6nAWSyxO2YQFwq+zkF/Cixnj
orU0T4k/EFYiBjwCCOY/6dxZwot5DmN4jskwDW38cY6mcH+WVGT5TbZHbLB0k4mqvZtm6X6z02Tp
6gltNptVL6TCfJ5ynUsfjtXCTDhyCS/WyzCGV0cnKhOmcK+dwmJEzovBHRNN+02WFysDuofElt5g
okGe4s7semf9Rh7p8nGpNtOYaJYXi5lrjBh7Uc0wYx0e8GIfb0bU4QvixScnJ2ZQtGGiYZ4ys08a
bCr2i83dSnyiavNxnN80I4QGL06AaRte/O3FQStJfkOWAUIc06jtxSkm6naOF8FE23BQ3zZL3Z/a
Rubixd8xFcNiJ/gybjQkUxTwDrWHA3wydHR0RJ4i0eb3fMWtxmIx7MmdkSAqE+PPz89HYPGwI6ee
VAKLy89ICNZLc3jA0AqxONh4jnUY1XfWFGlt1sHSDtsU2/eLWUEVOiN6FGrN8wbq0GJ/uVuiqSPb
OIuatsEYSr0BBRvdoc0eJSrxm2EETeTFJYzO2SN77micbbY/KeVe2px4ijbcg7XH7s7Pv3c3rUcx
bw4Pz1GoR6GiGpDI7pZV4cWTexRVbDO5R1Gllz6tR1G3XzyhR1GrgpnWo6h74mWz2ciHUAVQb0aN
F8yLcVgdxtstL86wrWl5akj8Rp/3U0CN1SZ7dli4qmSpMWa2Cb2YzMBiDelIYg+GBaLSUDCy1gcH
B8MymTF+qOIfZtmQr3yzezCFvNhkZyleLJQbuwdTyItdoGB+LOebKorFk6lWmhf7NNc8gWIzinpn
h3EdmAmq+AEYPf/rVwVjgzzLi0VLsp2X3xHu++kE1MnyYh9s8jXDltrMy4vn6ZjMtZG+1y1JlqXN
KrBzAOSNO/vrYO7GRymi/eI/pwv5Z3rx36pNEXduUCst63cwJb+7a6RNYU95zqyZ3Wz7n/3uTgDY
tXHhEu7cYqX++t/dbY83vyOr2WnHEu78rwADABaBbeIZChwYAAAAAElFTkSuQmCC';

my $CSS_STYLE = '
html, body, div, span, p, img {
    margin: 0;
    padding: 0;
    border: 0;
    outline: 0;
    font-size: 100%;
    vertical-align: baseline;
    background: transparent;
}

html, body {
    font-family: Arial, Verdana;
    color: #40454b;
    font-size: 12px;
    text-align: center;
}

img {
    padding: 0px; margin: 0px; border: none;
}

.info-panel {
    margin-top: 10px;
    margin-bottom: 10px;
    width: 740px;
    text-align: left;
}

.info-header {
    padding-top: 20px;
    padding-top: 10px;
}

.info-header-title {
    color: #126499;
    text-decoration: none;
    font-family: sans-serif;
    font-weight: bold;
    font-size: 16px;
    vertical-align: baseline;
    margin-right: 20px;
    margin-bottom: 25px;
    margin-top: 15px;
}

.info-content {
    padding: 2px;
    font-family: "lucida grande",sans-serif,arial;
    margin-top: 15px;
    margin-bottom: 15px;
}

.info-table-type {
    min-width: 70px;
    padding: 4px;
    vertical-align: top;
}

.info-table-value {
    font-weight: bold;
    padding-top: 4px;
    padding-left: 10px;
    padding-right: 10px;
    vertical-align: top;
}

hr {
    background-color: #E0E0E0;
    border: medium none;
    color: #E0E0E0;
    height: 1px;
    outline: medium none;
}

.sequencetext {
    font-family: courier, "courier new";
    font-weight: normal;
}
';

my $VERSION = '0.6';
my $WHAT = 'graphs';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print "PRINSEQ-$WHAT $VERSION\n"; exit; },
            'i=s',
	    'o=s',
	    'png_all',
	    'html_all',
            'log:s',
            'web:s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

PRINSEQ - PReprocessing and INformation of SEQuence data

=head1 VERSION

PRINSEQ-graphs 0.6

=head1 SYNOPSIS

perl prinseq-graphs.pl [-h] [-help] [-version] [-man] [-verbose] [-i input_graph_data_file] [-png_all] [-html_all] [-log file]

=head1 DESCRIPTION

PRINSEQ will help you to preprocess your genomic or metagenomic sequence data in FASTA (and QUAL) or FASTQ format. The graphs version allows users of the lite version to generate graphs similar to the web version.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-i> <file>

Input file containing the graph data generated by the lite version.

=item B<***** OUTPUT OPTIONS *****>

=item B<-o> <string>

By default, the output files are created in the same directory as the input file with an additional "_prinseq_graphs_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option. The file extension will be added automatically.

=item B<-png_all>

Use this option to generate PNG files with the graphs.

=item B<-html_all>

Use this option to generate a HTML file with the graphs and tables.

=item B<-log> <file>

Log file to keep track of parameters, errors, etc. The log file name is optional. If no file name is given, the log file name will be "inputname.log". If the log file already exists, new content will be added to the file.

=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >> or use http://sourceforge.net/tracker/?group_id=315449 so that I can make PRINSEQ better.

=head1 COPYRIGHT

Copyright (C) 2011-2012  Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

my ($file1,$command,@dataread);

#Check if input file exists and check if file format is correct
if(exists $params{i}) {
    $command .= ' -i '.$params{i};
    $file1 = $params{i};
    if($params{i} eq 'stdin') {
        my $format = &checkInputFormat();
        unless($format eq 'gd') {
            &printError('input data for -i is in '.uc($format).' format not in graph data format');
        }
    } elsif(-e $params{i}) {
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'gd') {
            &printError('input file for -i is in '.uc($format).' format not in graph data format');
        }
    } else {
        &printError("could not find input file \"".$params{i}."\"");
    }
} else {
    &printError("you did not specify an input file containing the graph data");
}

#check output file name prefix
if(exists $params{o}) {
    $command .= ' -o '.$params{o};
}

#check for output format
unless(exists $params{png_all} || exists $params{html_all}) {
    &printError("No output format specified. Use -png_all and/or -html_all to generate graphs.");
}
if(exists $params{png_all}) {
    $command .= ' -png_all';
}
if(exists $params{html_all}) {
    $command .= ' -html_all';
}
if(exists $params{web}) {
    $command .= ' -web'.($params{web} ? ' '.$params{web} : '');
}

#add remaining to log command
if(exists $params{log}) {
    $command .= ' -log'.($params{log} ? ' '.$params{log} : '');

    unless($params{log}) {
	$params{log} = join("__",$file1||'nonamegiven').'.log';
    }
    $params{log} = cwd().'/'.$params{log} unless($params{log} =~ /^\//);
    &printLog("Executing PRINSEQ with command: \"perl prinseq-".$WHAT.".pl".$command."\"");
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

my $filename = $file1;
while($filename =~ /[\w\d]+\.[\w\d]+$/) {
    $filename =~ s/\.[\w\d]+$//;
    last if($filename =~ /\/[^\.]+$/);
}

if(exists $params{png_all}) {
    my $graphs = &generateGraphs($params{i},$params{o});
    if(exists $params{web} && $params{web} ne 'nozip') {
        #png files
        if(scalar(@$graphs)) {
            system("zip -j -r ".dirname($params{o})."/png_graphs.zip ".dirname($params{o}).' -i \*.png') == 0 or &printError("Cannot generate graphs ZIP file");
        }
    }
}
if(exists $params{html_all}) {
    &generateHtml($params{i},$params{o});
}

&printWeb("STATUS: done");

##
#################################################################################
### MISC FUNCTIONS
#################################################################################
##

sub printError {
    my $msg = shift;
    print STDERR "\nERROR: ".$msg.".\n\nTry \'perl prinseq-".$WHAT.".pl -h\' for more information.\nExit program.\n";
    &printLog("ERROR: ".$msg.". Exit program.\n");
    exit(0);
}

sub printWarning {
    my $msg = shift;
    print STDERR "WARNING: ".$msg.".\n";
    &printLog("WARNING: ".$msg.".\n");
}

sub printWeb {
    my $msg = shift;
    if(exists $params{web}) {
        print STDERR "\n".&getTime()."$msg\n";
    }
}

sub getTime {
    return sprintf("[%02d/%02d/%04d %02d:%02d:%02d] ",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
}

sub printLog {
    my $msg = shift;
    if(exists $params{log}) {
        my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
        open(FH, ">>", $params{log}) or die "ERROR: Can't open file ".$params{log}.": $! \n";
        flock(FH, LOCK_EX) or die "ERROR: Cannot lock file ".$params{log}.": $! \n";
        print FH "[prinseq-".$WHAT."-$VERSION] [$time] $msg\n";
        flock(FH, LOCK_UN) or die "ERROR: cannot unlock ".$params{log}.": $! \n";
        close(FH);
    }
}

sub addCommas {
    my $num = shift;
    return unless(defined $num);
    return $num if($num < 1000);
    $num = scalar reverse $num;
    $num =~ s/(\d{3})/$1\,/g;
    $num =~ s/\,$//;
    $num = scalar reverse $num;
    return $num;
}

sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual,$gd,$aa);
    $count = 3;
    $fasta = $fastq = $qual = $gd = $aa = 0;
    $format = 'unknown';

    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/))) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/))) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        } elsif(!$gd && /^\{\"numseqs\"\:/) {
	    $gd = 1;
	}
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    } elsif($gd == 1) {
	$format = 'gd';
    }

    return $format;
}

sub checkInputFormat {
    my ($format,$count,$id,$fasta,$fastq,$qual,$gd,$aa);
    $count = 3;
    $fasta = $fastq = $qual = $gd = $aa = 0;
    $format = 'unknown';

    while (<STDIN>) {
        push(@dataread,$_);
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/))) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/))) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        } elsif(!$gd && /^\{\"numseqs\"\:/) {
	    $gd = 1;
	}
    }

    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    } elsif($gd == 1) {
	$format = 'gd';
    }

    return $format;
}

sub readGdFile {
    my $file = shift;
    my $data;

    open(DATA,"<$file") or &printError("Could not open file $file: $!");
    while(<DATA>) {
	next if(/^\#/);
        chomp();
	if(length($_)) {
	    $data = from_json($_);
	}
    }
    close(DATA);

    return $data;
}

sub getFileName {
    my $ext = shift;
    my ($file,$fh);
    if(exists $params{o}) {
	$file = $params{o}.$ext;
	open(OUT,">$file") or &printError('cannot open output file');
	close(OUT);
    } else {
	$fh = File::Temp->new( TEMPLATE => $filename.'_prinseq_graphs_XXXX',
			       SUFFIX => $ext,
			       UNLINK => 0);
	$file = $fh->filename;
	$fh->close();
    }
    return $file;
}

sub generateGraphs {
    my ($in,$out) = @_;
    my ($file,$data,$surface,@graphs);
    $data = &readGdFile($in);

    #length plot
    if(exists $data->{counts}->{length}) {
	$file = &getFileName('_ld.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{length},1),$data->{stats}->{length},'Length Distribution','Read Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{counts2} && exists $data->{counts2}->{length}) {
	$file = &getFileName('_ld-2.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{length},1),$data->{stats2}->{length},'Length Distribution','Read Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #tail plot
    if(exists $data->{tail}) {
	$file = &getFileName('_td5.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{tail5},1),undef,'Poly-A/T Tail Distribution (> 4bp)','5\' Tail Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
	$file = &getFileName('_td3.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{tail3},1),undef,'Poly-A/T Tail Distribution (> 4bp)','3\' Tail Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{tail2}) {
	$file = &getFileName('_td5-2.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{tail5},1),undef,'Poly-A/T Tail Distribution (> 4bp)','5\' Tail Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
	$file = &getFileName('_td3-2.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{tail3},1),undef,'Poly-A/T Tail Distribution (> 4bp)','3\' Tail Length in bp','# Sequences',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Ns plot
    if(exists $data->{counts}->{ns}) {
	$file = &getFileName('_ns.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{ns},1),undef,'Percentage of N\'s (> 0%)','Percentage of N\'s per Read (1-100%)','# Sequences',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{counts2} && exists $data->{counts2}->{ns}) {
	$file = &getFileName('_ns-2.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{ns},1),undef,'Percentage of N\'s (> 0%)','Percentage of N\'s per Read (1-100%)','# Sequences',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #GC content plot
    if(exists $data->{counts}->{gc}) {
	$file = &getFileName('_gc.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{gc},0),$data->{stats}->{gc},'GC Content Distribution','GC Content (0-100%)','Number of Sequences',$file,1);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{counts2} && exists $data->{counts2}->{gc}) {
	$file = &getFileName('_gc-2.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{gc},0),$data->{stats2}->{gc},'GC Content Distribution','GC Content (0-100%)','Number of Sequences',$file,1);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Sequence complexity plot - dust
    if(exists $data->{compldust}) {
	$file = &getFileName('_cd.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{compldust},0),undef,'Sequence complexity distribution','Mean sequence complexity (DUST scores)','Number of sequences',$file,1);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Sequence complexity plot - entropy
    if(exists $data->{complentropy}) {
	$file = &getFileName('_ce.png');
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{complentropy},0),undef,'Sequence complexity distribution','Mean sequence complexity (Entropy values)','Number of sequences',$file,1);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Dinucleotide odd ratio PCA plot - microbial/viral
    #Odds ratio plot
    if(exists $data->{dinucodds}) {
	my @new = map {$data->{dinucodds}->{$_}} sort keys %{$data->{dinucodds}};
	$file = &getFileName('_pm.png');
	$surface = &createPCAPlot(&convertToPCAValues(\@new,'m'),'PCA','1st Principal Component Score','2nd Principal Component Score',$file);
	$surface->write_to_png($file);
        push(@graphs,$file);
	$file = &getFileName('_pv.png');
	$surface = &createPCAPlot(&convertToPCAValues(\@new,'v'),'PCA','1st Principal Component Score','2nd Principal Component Score',$file);
	$surface->write_to_png($file);
        push(@graphs,$file);
	$file = &getFileName('_or.png');
	$surface = &createOddsRatioPlot($data->{dinucodds},'Odds ratios','Dinucleotide','Odds ratio',$file);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Qual plot
    if(exists $data->{quals}) {
	$file = &getFileName('_qd.png');
	$surface = &createBoxPlot(&convertToBoxValues($data->{quals},4),'Base Quality Distribution','Read position in %','Quality score',$file);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{quals2}) {
	$file = &getFileName('_qd-2.png');
	$surface = &createBoxPlot(&convertToBoxValues($data->{quals2},4),'Base Quality Distribution','Read position in %','Quality score',$file);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Qualbin plot
    if(exists $data->{qualsbin}) {
	$file = &getFileName('_qd2.png');
	$surface = &createBoxPlot(&convertToBoxValues($data->{qualsbin},4),'Base Quality Distribution','Read position in bp','Quality score',$file,0,'bp',$data->{binval});
       $surface->write_to_png($file);
       push(@graphs,$file);
    }
    if(exists $data->{qualsbin2}) {
	$file = &getFileName('_qd2-2.png');
	$surface = &createBoxPlot(&convertToBoxValues($data->{qualsbin2},4),'Base Quality Distribution','Read position in bp','Quality score',$file,0,'bp',$data->{binval});
       $surface->write_to_png($file);
       push(@graphs,$file);
    }

    #Qualmean plot
    if(exists $data->{qualsmean}) {
	$file = &getFileName('_qd3.png');
	$surface = &createBarPlot(&convertToBarValues($data->{qualsmean},5,1),'Sequence Quality Distribution','Mean of quality scores per sequence','Number of sequences',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{qualsmean2}) {
	$file = &getFileName('_qd3-2.png');
	$surface = &createBarPlot(&convertToBarValues($data->{qualsmean2},5,1),'Sequence Quality Distribution','Mean of quality scores per sequence','Number of sequences',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    #Sequence duplicate plots
    if(exists $data->{dubscounts}) {
	$file = &getFileName('_df.png');
	$surface = &createStackBarPlot(&convertOdToStackBinMatrix($data->{dubscounts},5,1,100),'Sequence duplication level','Number of duplicates','Number of sequences',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{dubslength}) {
	$file = &getFileName('_dl.png');
	$surface = &createStackBarPlot(&convertOdToStackBinMatrix($data->{dubslength},5,1),'Sequence duplication level','Read Length in bp','Number of duplicates',$file,0,' bp');
	$surface->write_to_png($file);
        push(@graphs,$file);
    }
    if(exists $data->{dubscounts}) {
	my %dubsmax;
	my $count = 1;
	foreach my $n (sort {$b <=> $a} keys %{$data->{dubscounts}}) {
	    foreach my $s (keys %{$data->{dubscounts}->{$n}}) {
		foreach my $i (1..$data->{dubscounts}->{$n}->{$s}) {
		    $dubsmax{$count++}->{$s} = $n;
		    last unless($count <= 100);
		}
		last unless($count <= 100);
	    }
	    last unless($count <= 100);
	}
	$file = &getFileName('_dm.png');
	$surface = &createStackBarPlot(&convertOdToStackBinMatrix(\%dubsmax,5,1,100),'Sequence duplication level','Sequence','Number of duplicates',$file,0);
	$surface->write_to_png($file);
        push(@graphs,$file);
    }

    return \@graphs;
}

sub convertOdToBinMatrix {
    my ($data,$min,$max,$nonice) = @_;

    my ($num,$ymax,$xmax,$xmin,$step,%vals,$tmp,@matrix,$bin,$tmpbin);

    #make nice xmax value
    if(defined $max) {
	$xmax = $max;
    } else {
	$xmax = (sort {$b <=> $a} keys %$data)[0];
    }
    $bin = &getBinVal($xmax);
    $xmax = $bin*100;
    $xmin = (defined $min ? $min : 0);

    #get data to bin and find y axis max value
    $ymax = 0;
    $tmp = 0;
    $tmpbin = $bin;
    foreach my $i ($xmin..$xmax) {
	if(exists $data->{$i}) {
	    $tmp += $data->{$i};
	}
	if(--$tmpbin <= 0) {
	    $tmpbin = $bin;
	    $ymax = &max($ymax,$tmp);
	    push(@matrix,$tmp);
	    $tmp = 0;
	}
    }

    #make nice ymax value
    unless($nonice) {
	$ymax = sprintf("%d",($ymax/4)+1)*4 if($ymax % 4);
#        $step = ($ymax <= 10 ? 10 : ($ymax < 40 ? 40 : ($ymax < 100 ? 100 : ($ymax < 1000 ? 100 : 100))));
#        $ymax = sprintf("%d",($ymax/$step)+1)*$step if($ymax % $step);
    }

    return (\@matrix,$xmax,$ymax);
}

sub getBinVal {
    my $val = shift;
    my $step;
    if(!$val || $val <= 100) {
	return 1;
    } elsif($val < 10000) {
	return int($val/100)+($val % 100 ? 1 : 0);
    } elsif($val < 100000) {
	return 1000;
    } else {
	$step = 1000000;
	my $xmax = ($val % $step ? sprintf("%d",($val/$step+1))*$step : $val);
	return ($xmax/100);
    }
}

sub max {
    my ($a,$b) = @_;
    return ($a < $b ? $b : $a);
}

sub min {
    my ($a,$b) = @_;
    return ($a > $b ? $b : $a);
}

sub createAnnotBarPlot {
    my ($matrix,$xmax,$ymax,$annot,$title,$xlab,$ylab,$file,$zero,$add) = @_;

    my $bin = 1;
    if($xmax > 100) {
	$bin = $xmax / 100;
	$xmax = 100;
    }

    my @barcol = (127/255, 127/255, 255/255, 1); #b2b2ff
    my @meancol = (255/255, 127/255, 127/255, 1); #ffb2b2
    my @stdcol = (178/255, 178/255, 255/255, 0.8); #7f7fff
    my @std1col = (0, 0, 0, 0.04); #ff7f7f
    my @std2col = (0, 0, 0, 0.03); #ff7f7f
    my @linecol = (0, 0, 0, 0.4);
    my @helplinecol = (1, 1, 1, 0.9);
    my @background = (0.95, 0.95, 0.95, 1);
    my @tickcol = (0, 0, 0, 0.8);
    my @labelcol = (0, 0, 0, 1);

    #create new image
    my $size   = 6;
    my $offset = 20;
    my $left   = 40;
    my $bottom = 15;
    my $top = 20;
    my $height = 200;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+($xmax+$zero)*$size,$bottom+$top+$offset*2+$height); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+2*200+20);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face ('sans', 'normal', 'normal');

    $cr->save;

    #set up work space
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, ($xmax+$zero)*$size-1, $height);
    $cr->set_source_rgba(@background);
    $cr->fill;

    #draw ticks
    #x-axis
    $cr->set_source_rgba(@tickcol);
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%5) == 0 && $i > 1 && $i < $xmax) {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+3);
	} else {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+1);
	}
    }
    $cr->stroke;

    #y-axis
    $cr->move_to($left+$offset, $top+$offset);
    $cr->line_to($left+$offset-3, $top+$offset);
    $cr->move_to($left+$offset, $top+$offset+$height-1);
    $cr->line_to($left+$offset-3, $top+$offset+$height-1);
    $cr->stroke;

    #helplines
    $cr->set_source_rgba(@helplinecol);
    foreach my $j (1..3) {
        $cr->move_to($left+$offset, $top+$offset+$height*$j/4-($j ? 1 : 0));
        $cr->line_to($left+$offset+($xmax+$zero)*$size, $top+$offset+$height*$j/4-($j ? 1 : 0));
    }
    $cr->stroke;

    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(@tickcol);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%10) == 0 && $i > 1 && $i < $xmax) {
	    $extents = $cr->text_extents($i*$bin);
	    $cr->move_to($left+$offset+int($size/2+1)+$size*$i-($zero ? 0 : $size)-$extents->{width}/2-1-($i == 1 ? 1 : 0), $top+$offset+$height+$fontheight+2);
	    $cr->show_text($i*$bin);
	}
    }
    #y-axis
    $extents = $cr->text_extents(&addCommas($ymax));
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2);
    $cr->show_text(&addCommas($ymax));
    $extents = $cr->text_extents(0);
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2+$height);
    $cr->show_text(0);

    $cr->save;

    #labels
    $cr->set_font_size (14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};

    #axis labels
    $cr->set_source_rgba(@labelcol);
    $extents = $cr->text_extents($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->move_to($left+$offset+($xmax+$zero)*$size/2-$extents->{width}/2, $top+$offset+$height+$fontheight+15);
    $cr->show_text($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab.($bin>1 ? ' (per bin)' : ''));
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2),$offset+10);
    $cr->show_text($ylab.($bin>1 ? ' (per bin)' : ''));

    $cr->restore;

    #draw annotations
    if($annot) {
	$cr->set_antialias('none');
	my ($std1l,$std2l,$std1r,$std2r);
	#std boxes
	$std1l = int($annot->{mean})-int($annot->{std});
	$std2l = int($annot->{mean})-2*int($annot->{std});
	$std1r = int($annot->{mean})+int($annot->{std});
	$std2r = int($annot->{mean})+2*int($annot->{std});
        unless($std1l == $std1r) {
            if($std1l < 0) {
                $std1l = 0;
            } else {
                $std1l = int($std1l/$bin);
            }
            if($std2l < 0) {
                $std2l = 0;
            } else {
                $std2l = int($std2l/$bin);
            }
            if($std1r/$bin > 100) {
                $std1r = 100;
            } else {
                $std1r = int($std1r/$bin);
            }
            if($std2r/$bin > 100) {
                $std2r = 100;
            } else {
                $std2r = int($std2r/$bin);
            }
            $cr->rectangle($left+$offset+$std2l*$size+2, $top+$offset, ($std2r-$std2l)*$size, $height);
            $cr->set_source_rgba(@std2col);
            $cr->fill;
            $cr->rectangle($left+$offset+$std1l*$size+2, $top+$offset, ($std1r-$std1l)*$size, $height);
            $cr->set_source_rgba(@std1col);
            $cr->fill;
            #mean line
            $cr->set_source_rgba(@meancol);
            $cr->move_to($left+$offset+int(int($annot->{mean})/$bin)*$size+2, $top+$offset-5);
            $cr->line_to($left+$offset+int(int($annot->{mean})/$bin)*$size+2, $top+$offset+$height);
            $cr->stroke;
            #std lines
            $cr->set_source_rgba(@stdcol);
            if($std1l > 0) {
                $cr->move_to($left+$offset+$std1l*$size+2, $top+$offset-5);
                $cr->line_to($left+$offset+$std1l*$size+2, $top+$offset+$height);
            }
            if($std2l > 0) {
                $cr->move_to($left+$offset+$std2l*$size+2, $top+$offset-5);
                $cr->line_to($left+$offset+$std2l*$size+2, $top+$offset+$height);
            }
            if($std1r < 100) {
                $cr->move_to($left+$offset+$std1r*$size+2, $top+$offset-5);
                $cr->line_to($left+$offset+$std1r*$size+2, $top+$offset+$height);
            }
            if($std2r < 100) {
                $cr->move_to($left+$offset+$std2r*$size+2, $top+$offset-5);
                $cr->line_to($left+$offset+$std2r*$size+2, $top+$offset+$height);
            }
            $cr->stroke;
            #labels
            $cr->set_antialias('default');
            $cr->set_source_rgba(@tickcol);
            $extents = $cr->text_extents('M');
            $cr->move_to($left+$offset+int(int($annot->{mean})/$bin)*$size+2-$extents->{width}/2, $top+$offset-10);
            $cr->show_text('M');
            if($std1l > 0) {
                $extents = $cr->text_extents('1SD');
                $cr->move_to($left+$offset+$std1l*$size-$extents->{width}/2+2, $top+$offset-10);
                $cr->show_text('1SD');
            }
            if($std2l > 0) {
                $extents = $cr->text_extents('2SD');
                $cr->move_to($left+$offset+$std2l*$size-$extents->{width}/2+3, $top+$offset-10);
                $cr->show_text('2SD');
            }
            if($std1r < 100) {
                $extents = $cr->text_extents('1SD');
                $cr->move_to($left+$offset+$std1r*$size-$extents->{width}/2+2, $top+$offset-10);
                $cr->show_text('1SD');
            }
            if($std2r < 100) {
                $extents = $cr->text_extents('2SD');
                $cr->move_to($left+$offset+$std2r*$size-$extents->{width}/2+3, $top+$offset-10);
                $cr->show_text('2SD');
            }
        }
    }

    #draw boxes
    $cr->set_antialias('none');
    $cr->set_source_rgba(@barcol);
    foreach my $pos (0..$xmax-($zero ? 0 : 1)) {
        next unless($matrix->[$pos]);
        my $tmp = $matrix->[$pos] / $ymax;
	#unique
        if($tmp) {
            $cr->rectangle($left+$offset+$pos*$size, $top+$offset+$height, $size-1, -$tmp*$height);
            $cr->fill;
        }
    }

    #write image
    $cr->show_page;
    return $surface;
}

sub convertToPCAValues {
    my ($new,$type) = @_;

    my @data = ($type eq 'v' ? @$DINUCODDS_VIR : @$DINUCODDS_MIC);

    push(@data,$new);

    my $pca = Statistics::PCA->new;

    #suppress output from PCA module
    my $output = '';
    open(TOOUTPUT, '>', \$output) or &printError("Can't open TOOUTPUT: $!");
    select TOOUTPUT;

    $pca->load_data({format => 'table', data => \@data});
    $pca->pca();

    my @variances = $pca->results('proportion');
    my @list = $pca->results('transformed');

    #end suppress output from PCA module
    select STDOUT;
    close(TOOUTPUT);

    my ($xmin,$xmax,$ymin,$ymax);
    $xmax = $ymax = -100;
    $xmin = $ymin = 100;

    #get min/max values for PC1
    foreach my $v (@{$list[0]}) {
        $xmax = &max($xmax,$v);
        $xmin = &min($xmin,$v);
    }
    #get min/max values for PC2
    foreach my $v (@{$list[1]}) {
        $ymax = &max($ymax,$v);
        $ymin = &min($ymin,$v);
    }

    return ([$list[0],$list[1]],sprintf("%d",$variances[0]*100),sprintf("%d",$variances[1]*100),$xmin,$xmax,$ymin,$ymax,$type);
}

sub createPCAPlot {
    my ($data,$var1,$var2,$xmin,$xmax,$ymin,$ymax,$type,$title,$xlab,$ylab,$file) = @_;

    my @linecol = (0, 0, 0, 0.4);
    my @helplinecol1 = (1,1,1, 0.9);
    my @helplinecol2 = (1,1,1, 0.5);

    #create new image
    my $size   = 5;
    my $offset = 20;
    my $left   = 25;
    my $bottom = 15;
    my $top = ($type eq 'v' ? 35 : 20);
    my $height = 500;
    my $space = 10;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+$height+2*$space,$top+$bottom+$offset*2+$height+2*$space); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+$height+2*$space,$top+$bottom+$offset*2+$height+2*$space);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face('sans-serif', 'normal', 'normal');

    $cr->save;

    #set up work space
    my ($dx, $dy);
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, $height+2*$space, $height+2*$space);
    $cr->set_source_rgba(0.95, 0.95, 0.95, 1);
    $cr->fill;

    #get infos
    my $num = scalar(@{$data->[0]})-1;
    my $xrange = ($xmax-$xmin);
    my $yrange = ($ymax-$ymin);
    my $data_info = ($type eq 'v' ? $DATA_VIR : $DATA_MIC);

    #draw ticks
    #x-axis
    $cr->set_source_rgba(0, 0, 0, 0.8);
    $cr->move_to($left+$offset+$space, $top+$offset+$height+2*$space);
    $cr->line_to($left+$offset+$space, $top+$offset+$height+2*$space+3);
    $cr->move_to($left+$offset+$space+$height, $top+$offset+$height+2*$space);
    $cr->line_to($left+$offset+$space+$height, $top+$offset+$height+2*$space+3);
    $cr->move_to($left+$offset+$space+int(abs($xmin)/$xrange*$height), $top+$offset+$height+2*$space);
    $cr->line_to($left+$offset+$space+int(abs($xmin)/$xrange*$height), $top+$offset+$height+2*$space+3);
    $cr->stroke;
    #y-axis
    $cr->move_to($left+$offset, $top+$offset+$space);
    $cr->line_to($left+$offset-3, $top+$offset+$space);
    $cr->move_to($left+$offset, $top+$offset+$height+$space);
    $cr->line_to($left+$offset-3, $top+$offset+$height+$space);
    $cr->move_to($left+$offset, $top+$offset+$space+int(abs($ymax)/$yrange*$height));
    $cr->line_to($left+$offset-3, $top+$offset+$space+int(abs($ymax)/$yrange*$height));
    $cr->stroke;

    #helplines
    $cr->set_source_rgba(@helplinecol1);
    $cr->move_to($left+$offset+$space+int(abs($xmin)/$xrange*$height), $top+$offset);
    $cr->line_to($left+$offset+$space+int(abs($xmin)/$xrange*$height), $top+$offset+$height+2*$space);
    $cr->stroke;
    $cr->move_to($left+$offset, $top+$offset+$space+int(abs($ymax)/$yrange*$height));
    $cr->line_to($left+$offset+2*$space+$height, $top+$offset+$space+int(abs($ymax)/$yrange*$height));
    $cr->stroke;


    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(0, 0, 0, 0.8);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
    $extents = $cr->text_extents(sprintf("%.2f",$xmin));
    $cr->move_to($left+$offset+$space-$extents->{width}/2-1, $top+$offset+$height+2*$space+$fontheight+2);
    $cr->show_text(sprintf("%.2f",$xmin));
    $extents = $cr->text_extents(sprintf("%.2f",$xmax));
    $cr->move_to($left+$offset+$space+$height-$extents->{width}/2-1, $top+$offset+$height+2*$space+$fontheight+2);
    $cr->show_text(sprintf("%.2f",$xmax));
    $extents = $cr->text_extents(0);
    $cr->move_to($left+$offset+$space+int(abs($xmin)/$xrange*$height)-$extents->{width}/2, $top+$offset+$height+2*$space+$fontheight+2);
    $cr->show_text(0);
    #y-axis
    $extents = $cr->text_extents(sprintf("%.2f",$ymax));
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$space+$fontheight/2-2);
    $cr->show_text(sprintf("%.2f",$ymax));
    $extents = $cr->text_extents(sprintf("%.2f",$ymin));
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$height+$space+$fontheight/2-2);
    $cr->show_text(sprintf("%.2f",$ymin));
    $extents = $cr->text_extents(0);
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$space+int(abs($ymax)/$yrange*$height)+$fontheight/2-2);
    $cr->show_text(0);

    $cr->save;

    #labels
    $cr->set_font_size (14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};

    #add type
    $cr->set_source_rgba(0, 0, 0, 0.5);
    $extents = $cr->text_extents(uc($type));
    $cr->arc($offset/2+$extents->{width}/2, $offset-5, 10, 0, 2*$PI);
    $cr->fill;
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->move_to($offset/2-($type eq 'm' ? 1 : 0), $offset);
    $cr->show_text(uc($type));

    #axis labels
    $cr->set_source_rgba(0, 0, 0, 1);
    $extents = $cr->text_extents($xlab.'  ('.$var1.'%)');
    $cr->move_to($left+$offset+$height/2-$extents->{width}/2+$space, $top+$offset+$height+$fontheight+15+2*$space);
    $cr->show_text($xlab.'  ('.$var1.'%)');
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab.'  ('.$var2.'%)');
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2)+$space,$offset);
    $cr->show_text($ylab.'  ('.$var2.'%)');

    $cr->restore;

    #draw dots
    $cr->set_antialias('default');
    $cr->set_font_size (10);
    foreach my $i (0..$num) {
        $cr->set_source_rgba(@{$data_info->[$i]->[3]});
        $cr->arc(($left+$offset+$space+int(($data->[0]->[$i]+abs($xmin))/$xrange*$height)), ($space+$top+$offset+int(($data->[1]->[$i]+abs($ymin))/$yrange*$height)), $size, 0, 2*$PI);
        $cr->fill;
    }
    $cr->set_source_rgba(0, 0, 0, 1);
    foreach my $i (0..$num) {
        $extents = $cr->text_extents($data_info->[$i]->[1]);
	$cr->move_to(($left+$offset+$space+int(($data->[0]->[$i]+abs($xmin))/$xrange*$height))+$size+1, ($space+$top+$offset+int(($data->[1]->[$i]+abs($ymin))/$yrange*$height))+$size*2);
	$cr->show_text($data_info->[$i]->[1]);
    }

    #draw legend
    my %labels;
    foreach  my $i (0..$num) {
	$labels{$data_info->[$i]->[1]} = $data_info->[$i]->[2];
    }
    $cr->set_font_size(10);
    $fontheight = $font_extents->{height};
    $cr->set_source_rgba(0, 0, 0, 1);
    my $x = $left+$offset+$space;
    my $y = int($offset/2);
    foreach my $n (sort {$a <=> $b} keys %labels) {
	if($x+$cr->text_extents($n.' - '.$labels{$n})->{width}+15 >= $left+$offset+$space+$height) {
	    $x = $left+$offset+$space;
	    $y += $fontheight;
	}
	$cr->move_to($x,$y);
	$cr->show_text($n.' - '.$labels{$n});
	$x += $cr->text_extents($n.' - '.$labels{$n})->{width}+15;

    }

    #write image
    $cr->show_page;
    return $surface;
}

sub createOddsRatioPlot {
    my ($data,$title,$xlab,$ylab,$file) = @_;

    my @yvalues = (0.5,0.78,1.00,1.23,1.5);

    my @linecol = (0, 0, 0, 0.4);
    my @helplinecol1 = (1,1,1, 0.9);
    my @helplinecol2 = (1,1,1, 0.5);

    #create new image
    my $size   = 40;
    my $offset = 20;
    my $left   = 35;
    my $right = 90;
    my $bottom = 20;
    my $top = 0;
    my $height = 100;
    my $width = $size*10;
    my $space = 20;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+$width+$right,$top+$bottom+$offset*2+$height); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+$width+$right,$top+$bottom+$offset*2+$height);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face ('sans', 'normal', 'normal');

    $cr->save;

    #set up work space
    my ($dx, $dy);
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, $width, $height);
    $cr->set_source_rgba(0.95, 0.95, 0.95, 1);
    $cr->fill;

    #right side marks
    $cr->set_source_rgba(255/255, 127/255, 127/255, 0.6);
    $cr->rectangle($left+$offset+$width+8, $top+$offset, 3, 0.77/2*$height);
    $cr->fill;
    $cr->rectangle($left+$offset+$width+8, $top+$offset+$height-0.78/2*$height, 3, 0.78/2*$height);
    $cr->fill;

    #get infos
    my $num = scalar(keys %$data)-1;

    #draw ticks
    #x-axis
    $cr->set_source_rgba(0, 0, 0, 0.8);
    foreach my $i (0..$num) {
	$cr->move_to($left+$offset+$size/2+$i*$size, $top+$offset+$height);
	$cr->line_to($left+$offset+$size/2+$i*$size, $top+$offset+$height+3);
    }
    $cr->stroke;
    #y-axis
    foreach my $i (@yvalues) {
	$cr->move_to($left+$offset, $top+$offset+$height-$i/2*$height);
	$cr->line_to($left+$offset-3, $top+$offset+$height-$i/2*$height);
    }
    $cr->stroke;

    #helplines
    #x-axis
    $cr->set_source_rgba(@helplinecol1);
    foreach my $i (0..$num) {
	$cr->move_to($left+$offset+$size/2+$i*$size, $top+$offset);
	$cr->line_to($left+$offset+$size/2+$i*$size, $top+$offset+$height);
    }
    $cr->stroke;
    #yaxis
    foreach my $i (@yvalues) {
	$cr->set_source_rgba(0, 0, 0, ($i == 0.5 || $i == 1.00 || $i == 1.50 ? 0.1 : 0.3));
	$cr->move_to($left+$offset, $top+$offset+$height-$i/2*$height);
	$cr->line_to($left+$offset+$width, $top+$offset+$height-$i/2*$height);
	$cr->stroke;
    }

    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(0, 0, 0, 0.8);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
    my $xcur = 0;
    foreach my $dn (map {join("/",(m/../g ))} sort keys %$data) {
	$extents = $cr->text_extents($dn);
	$cr->move_to($left+$offset+$size/2-$extents->{width}/2-1+$size*$xcur++, $top+$offset+$height+$fontheight+2);
	$cr->show_text($dn);
    }
    #y-axis
    foreach my $i (@yvalues) {
	$cr->set_source_rgba(0, 0, 0, ($i == 0.5 || $i == 1.00 || $i == 1.50 ? 0.5 : 0.8));
	$extents = $cr->text_extents(sprintf("%.2f",$i));
	$cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$height-$i/2*$height+$fontheight/2-2);
	$cr->show_text(sprintf("%.2f",$i));
    }

    #label on right side
    $cr->set_source_rgba(0, 0, 0, 1);
    $extents = $cr->text_extents('Over-represented');
    $cr->move_to($left+$offset+$width+15, $top+$offset+$height-1.6/2*$height+$fontheight/2-2);
    $cr->show_text('Over-represented');
    $extents = $cr->text_extents('Under-represented');
    $cr->move_to($left+$offset+$width+15, $top+$offset+$height-0.4/2*$height+$fontheight/2-2);
    $cr->show_text('Under-represented');

    $cr->save;

    #labels
    $cr->set_font_size (14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};

    #axis labels
    #x-axis
    $cr->set_source_rgba(0, 0, 0, 1);
    $extents = $cr->text_extents($xlab);
    $cr->move_to($left+$offset+$width/2-$extents->{width}/2, $top+$offset+$height+$fontheight+15);
    $cr->show_text($xlab);
    #y-axis
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab);
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2),$offset);
    $cr->show_text($ylab);

    $cr->restore;

    #draw dots
    $cr->set_antialias('default');
    $xcur = 0;
    foreach my $dn (sort keys %$data) {
	if($data->{$dn} > 1.23 || $data->{$dn} < 0.78) {
	    $cr->set_source_rgba(255/255, 127/255, 127/255, 1);
	} else {
	    $cr->set_source_rgba(127/255, 127/255, 255/255, 1);
	}
        $cr->arc($left+$offset+$size/2+$size*$xcur++, $top+$offset+$height-$data->{$dn}/2*$height, 5, 0, 2*$PI);
        $cr->fill;
    }

    #write image
    $cr->show_page;
    return $surface;
}

sub convertToBoxValues {
    my ($data,$niceval) = @_;
    my ($xmax,$ymax,@matrix);
    $xmax = $ymax = 0;
    foreach my $i (sort {$a <=> $b} keys %$data) {
	$xmax++;
	push(@matrix,[$i,$data->{$i}->{min},$data->{$i}->{p25},$data->{$i}->{median},$data->{$i}->{p75},$data->{$i}->{max}]);
	$ymax = &max($ymax,$data->{$i}->{max});
    }

    if($niceval) {
        $ymax = sprintf("%d",($ymax/$niceval)+1)*$niceval if($ymax % $niceval);
    }

    return (\@matrix,$xmax,$ymax);
}

sub createBoxPlot {
    my ($matrix,$xmax,$ymax,$title,$xlab,$ylab,$file,$zero,$add,$bin) = @_;
    $bin = ($bin ? $bin : 1);
    $zero = 0 unless($zero);
    $add = '' unless(defined $add);
    if($xmax != 100) {
	$xmax = 100;
    }
    $ymax = 1 unless($ymax);
#    die Dumper $matrix;


    my @col0 = (178/255, 178/255, 255/255); #b2b2ff
    my @col1 = (255/255, 178/255, 178/255); #ffb2b2
    my @col3 = (127/255, 127/255, 255/255); #7f7fff
    my @col4 = (255/255, 127/255, 127/255); #ff7f7f
    my @linecol = (0, 0, 0, 0.4);
    my @linecol0 = (@col3, 1);
    my @linecol1 = (@col4, 1);
    my @boxcol = (@col3, 1);
    my @whiscol = (@col0, 0.9);
    my @medcol = (0,0,0, 0.5);
    my @helplinecol1 = (1,1,1, 0.9);
    my @helplinecol2 = (1,1,1, 0.5);

    #create new image
    my $size   = 6;
    my $offset = 20;
    my $left   = 25;
    my $bottom = 25;
    my $top = 5;
    my $height = 300;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+$height); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+2*200+20);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face ('sans', 'normal', 'normal');
#    $cr->set_font_size (30);

    $cr->save;

    #set up work space
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, ($xmax+$zero)*$size-1, $height);
    $cr->set_source_rgba(0.95, 0.95, 0.95, 1);
    $cr->fill;

    #draw legend
    $cr->set_font_size(10);
#    $font_extents = $cr->font_extents;
    my $x = $left+$offset+$size*50;
    foreach my $v ([\@whiscol,'Min/Max value'],[\@boxcol,'25th to 75th percentile'],[\@medcol,'Median']) {
	$cr->set_antialias('none');
	$cr->set_source_rgba(@{$v->[0]});
	$cr->rectangle($x, $top+5, 10, 10);
	$cr->fill;
	$x += 15;
	$cr->set_antialias('default');
	$cr->move_to($x,$top+5+9);
	$cr->set_source_rgba(0, 0, 0, 0.8);
	$cr->show_text($v->[1]);
	$x += $cr->text_extents($v->[1])->{width}+15;
    }

    $cr->set_antialias('none');

    #draw ticks
    #x-axis
    $cr->set_source_rgba(0, 0, 0, 0.8);
#    $cr->move_to($left+$offset+int($size/2+1), $top+$offset+$height);
#    $cr->line_to($left+$offset+int($size/2+1), $top+$offset+$height+3);
#    $cr->move_to($left+$offset+int($size/2+1), $top+$offset+$height+$space);
#    $cr->line_to($left+$offset+int($size/2+1), $top+$offset+$height+$space-3);
    foreach my $i (1..9) {
        $cr->move_to($left+$offset+int($size/2)+$size*10*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	$cr->line_to($left+$offset+int($size/2)+$size*10*$i-($zero ? 0 : $size)-1, $top+$offset+$height+3);
#	$cr->move_to($left+$offset+int($size/2)+$size*10*$i-($zero ? 0 : $size)-1, $top+$offset);
#	$cr->line_to($left+$offset+int($size/2)+$size*10*$i-($zero ? 0 : $size)-1, $top+$offset-3);
    }
    $cr->stroke;
    #y-axis
    foreach my $j (0..4) {
	$cr->move_to($left+$offset, $top+$offset+$height*$j/4-($j ? 1 : 0));
	$cr->line_to($left+$offset-3, $top+$offset+$height*$j/4-($j ? 1 : 0));
#	$cr->move_to($left+$offset+($xmax+$zero)*$size, $top+$offset+$height*$j/4-($j ? 1 : 0));
#	$cr->line_to($left+$offset+($xmax+$zero)*$size+3, $top+$offset+$height*$j/4-($j ? 1 : 0));
    }
    $cr->stroke;

    #helplines
    $cr->set_source_rgba(@helplinecol1);
    foreach my $j (1..3) {
        $cr->move_to($left+$offset, $top+$offset+$height*$j/4-($j ? 1 : 0));
        $cr->line_to($left+$offset+($xmax+$zero)*$size, $top+$offset+$height*$j/4-($j ? 1 : 0));
    }
    $cr->stroke;

    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(0, 0, 0, 0.8);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
#    $extents = $cr->text_extents(1);
#    $cr->move_to($left+$offset+int($size/2+1)-$extents->{width}, $top+$offset+$height+$fontheight+2);
#    $cr->show_text(1);
    foreach my $i (1..9) {
        $extents = $cr->text_extents($i*10*$bin);
	$cr->move_to($left+$offset+int($size/2+1)+$size*10*$i-($zero ? 0 : $size)-$extents->{width}/2-1-($i == 1 ? 1 : 0), $top+$offset+$height+$fontheight+2);
        $cr->show_text($i*10*$bin);
    }
    #y-axis
    foreach my $j (0..4) {
	$extents = $cr->text_extents(&addCommas($ymax*$j/4));
	$cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2+$height*(4-$j)/4);
	$cr->show_text(&addCommas($ymax*$j/4));
    }

    $cr->save;

    #axis labels
    $cr->set_source_rgba(0, 0, 0, 1);
    $cr->set_font_size (14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    $extents = $cr->text_extents($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->move_to($left+$offset+($xmax+$zero)*$size/2-$extents->{width}/2, $top+$offset+$height+$fontheight+15);
    $cr->show_text($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab);
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2),$offset);
    $cr->show_text($ylab);

    $cr->restore;

    #draw boxes
    my $factor = $height/$ymax;
    $cr->set_antialias('none');
    foreach my $v (@$matrix) {
	#wiskers
	$cr->set_source_rgba(@whiscol);
	if($v->[1] != $v->[2]) {
	    $cr->move_to($left+$offset+$size*$v->[0]+1, $top+$offset+$height-$v->[1]*$factor-1);
	    $cr->line_to($left+$offset+$size*$v->[0]+$size-2, $top+$offset+$height-$v->[1]*$factor-1);
	    $cr->stroke;
	}
	if($v->[4] != $v->[5]) {
	    $cr->move_to($left+$offset+$size*$v->[0]+1, $top+$offset+$height-$v->[5]*$factor);
	    $cr->line_to($left+$offset+$size*$v->[0]+$size-2, $top+$offset+$height-$v->[5]*$factor);
	    $cr->stroke;
	}
	$cr->save;
	$cr->set_dash(1,4,3);
	if($v->[1] != $v->[2]) {
	    $cr->move_to($left+$offset+$size*$v->[0]+int($size/2)-1, $top+$offset+$height-$v->[2]*$factor);
	    $cr->line_to($left+$offset+$size*$v->[0]+int($size/2)-1, $top+$offset+$height-$v->[1]*$factor);
	    $cr->stroke;
	}
	if($v->[4] != $v->[5]) {
	    $cr->move_to($left+$offset+$size*$v->[0]+int($size/2)-1, $top+$offset+$height-$v->[5]*$factor);
	    $cr->line_to($left+$offset+$size*$v->[0]+int($size/2)-1, $top+$offset+$height-$v->[4]*$factor-1);
	    $cr->stroke;
	}
	$cr->restore;
	#box
	if(($v->[2] != $v->[3]) || ($v->[4] != $v->[3])) {
	    $cr->set_source_rgba(@whiscol);
	    $cr->rectangle($left+$offset+$size*$v->[0], $top+$offset+$height-$v->[2]*$factor, $size-1, -($v->[4]-$v->[2])*$factor);
	    $cr->fill;
	    $cr->stroke;
	    $cr->set_source_rgba(@boxcol);
	    $cr->rectangle($left+$offset+$size*$v->[0], $top+$offset+$height-$v->[2]*$factor, $size-2, -($v->[4]-$v->[2])*$factor);
	    $cr->stroke;
	} else {
	    $cr->set_source_rgba(@boxcol);
	    $cr->move_to($left+$offset+$size*$v->[0], $top+$offset+$height-$v->[3]*$factor);
	    $cr->line_to($left+$offset+$size*$v->[0]+$size-1, $top+$offset+$height-$v->[3]*$factor);
	    $cr->stroke;
	}
	#median
	$cr->set_source_rgba(@medcol);
	$cr->move_to($left+$offset+$size*$v->[0]+1, $top+$offset+$height-$v->[3]*$factor);
	$cr->line_to($left+$offset+$size*$v->[0]+$size-2, $top+$offset+$height-$v->[3]*$factor);
	$cr->stroke;
    }

    #write image
    $cr->show_page;
    return $surface;
}

sub convertToBarValues {
    my ($data,$niceval,$start,$max) = @_;
    my ($xmax,$ymax,@matrix,$tmp);
    $xmax = $ymax = 0;

    #get xmax value
    if($max) {
	$xmax = $max;
    } else {
	foreach my $q (keys %$data) {
	    $xmax = &max($xmax,$q);
	}
    }
    if($niceval) {
        $xmax = sprintf("%d",($xmax/$niceval)+1)*$niceval if($xmax % $niceval);
    }

    #get matrix values
    foreach my $q ($start..$xmax) {
	$tmp = (exists $data->{$q} ? $data->{$q} : 0);
	$ymax = &max($ymax,$tmp);
	push(@matrix,$tmp);
    }

    $ymax = sprintf("%d",($ymax/4)+1)*4 if($ymax % 4);

    return (\@matrix,$xmax,$ymax);
}

sub createBarPlot {
    my ($matrix,$xmax,$ymax,$title,$xlab,$ylab,$file,$zero) = @_;

    my @col0 = (178/255, 178/255, 255/255); #b2b2ff
    my @col1 = (255/255, 178/255, 178/255); #ffb2b2
    my @col3 = (127/255, 127/255, 255/255); #7f7fff
    my @col4 = (255/255, 127/255, 127/255); #ff7f7f
    my @linecol = (0, 0, 0, 0.4);
    my @linecol0 = (@col3, 1);
    my @linecol1 = (@col4, 1);
    my @barcol0 = (@col3, 1);
    my @barcol1 = (@col4, 1);
    my @helplinecol1 = (1,1,1, 0.9);
    my @helplinecol2 = (1,1,1, 0.5);

    #create new image
    my $size   = ($xmax <= 50 ? 10 : ($xmax <= 100 ? 6 : 3));
    my $offset = 20;
    my $left   = 25;
    my $bottom = 15;
    my $top = 0;
    my $height = 200;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+$height); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+2*200+20);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face ('sans', 'normal', 'normal');
#    $cr->set_font_size (30);

    $cr->save;

    #set up work space
    my ($dx, $dy);
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, ($xmax+$zero)*$size-1, $height);
    $cr->set_source_rgba(0.95, 0.95, 0.95, 1);
    $cr->fill;

    #draw ticks
    #x-axis
    $cr->set_source_rgba(0, 0, 0, 0.8);
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%5) == 0 && $i > 1 && $i < $xmax) {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+3);
	} else {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+1);
	}
    }
    $cr->stroke;

    #y-axis
    $cr->move_to($left+$offset, $top+$offset);
    $cr->line_to($left+$offset-3, $top+$offset);
    $cr->move_to($left+$offset, $top+$offset+$height-1);
    $cr->line_to($left+$offset-3, $top+$offset+$height-1);
    $cr->stroke;

    #helplines
    $cr->set_source_rgba(@helplinecol1);
    foreach my $j (1..3) {
        $cr->move_to($left+$offset, $top+$offset+$height*$j/4-($j ? 1 : 0));
        $cr->line_to($left+$offset+($xmax+$zero)*$size, $top+$offset+$height*$j/4-($j ? 1 : 0));
    }
    $cr->stroke;

    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(0, 0, 0, 0.8);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%5) == 0 && $i > 1 && $i < $xmax) {
	    $extents = $cr->text_extents($i);
	    $cr->move_to($left+$offset+int($size/2+1)+$size*$i-($zero ? 0 : $size)-$extents->{width}/2-1-($i == 1 ? 1 : 0), $top+$offset+$height+$fontheight+2);
	    $cr->show_text($i);
	}
    }
    #y-axis
    $extents = $cr->text_extents(&addCommas($ymax));
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2);
    $cr->show_text(&addCommas($ymax));
    $extents = $cr->text_extents(0);
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2+$height);
    $cr->show_text(0);

    $cr->save;

    #labels
    $cr->set_font_size (14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};

    #axis labels
    $cr->set_source_rgba(0, 0, 0, 1);
    $extents = $cr->text_extents($xlab);
    $cr->move_to($left+$offset+($xmax+$zero)*$size/2-$extents->{width}/2, $top+$offset+$height+$fontheight+15);
    $cr->show_text($xlab);
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab);
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2),$offset);
    $cr->show_text($ylab);

    $cr->restore;

    #draw boxes
    $cr->set_antialias('none');
    foreach my $pos (0..$xmax-($zero ? 0 : 1)) {
        next unless($matrix->[$pos+($zero ? 0 : 1)]);
        my $tmp = $matrix->[$pos+($zero ? 0 : 1)] / $ymax;
	#unique
        if($tmp) {
            $cr->set_source_rgba(@barcol0);
            $cr->rectangle($left+$offset+($pos+($zero ? 0 : 1))*$size, $top+$offset+$height, $size-1, -$tmp*$height);
            $cr->fill;
        }
    }

    #write image
    $cr->show_page;
    return $surface;
}

sub convertOdToStackBinMatrix {
    my ($data,$stacks,$min,$max,$nonice) = @_;

    my ($num,$ymax,$xmax,$xmin,$step,%vals,%sums,$sum,@matrix,$bin,$tmpbin);

    #make nice xmax value
    if(defined $max) {
	$xmax = $max;
    } else {
	$xmax = (sort {$b <=> $a} keys %$data)[0];
    }
    $bin = &getBinVal($xmax);
    $xmax = $bin*100;
    $xmin = (defined $min ? $min : 0);

    #get data to bin and find y axis max value
    $ymax = 0;
    foreach my $s (0..$stacks-1) {
	$sums{$s} = 0;
    }
    $sum = 0;
    $tmpbin = $bin;
    foreach my $i ($xmin..$xmax) {
	foreach my $s (0..$stacks-1) {
	    next unless(exists $data->{$i}->{$s});
	    $sums{$s} += $data->{$i}->{$s};
	    $sum += $data->{$i}->{$s};
	}
	if(--$tmpbin <= 0) {
	    $tmpbin = $bin;
	    $ymax = &max($ymax,$sum);
	    $sum = 0;
	    foreach my $s (0..$stacks-1) {
		push(@{$matrix[$s]},$sums{$s});
		$sums{$s} = 0;
	    }
	}
    }

    #make nice ymax value
    unless($nonice) {
	$ymax = sprintf("%d",($ymax/4)+1)*4 if($ymax % 4);
#        $step = ($ymax <= 10 ? 10 : ($ymax < 40 ? 40 : ($ymax < 100 ? 100 : ($ymax < 1000 ? 100 : 100))));
#        $ymax = sprintf("%d",($ymax/$step)+1)*$step if($ymax % $step);
    }

    return (\@matrix,$xmax,$ymax,$stacks);
}

sub createStackBarPlot {
    my ($matrix,$xmax,$ymax,$stacks,$title,$xlab,$ylab,$file,$zero,$add) = @_;

    my $bin = 1;
    if($xmax > 100) {
	$bin = $xmax / 100;
	$xmax = 100;
    }

    my @legend = ('Exact dupl.','5\' dupl.','3\' dupl.','Rev. compl. exact dupl.','Rev. compl. 5\'/3\' dupl.');
    my @cols = ([69/255, 114/255, 167/255, 1],
		[137/255, 1165/255, 78/255, 1],
		[170/255, 70/255, 67/255, 1],
		[147/255, 169/255, 207/255, 1],
		[51/255, 102/255, 102/255, 1]);
    my @barcol = (127/255, 127/255, 255/255, 1); #b2b2ff
    my @meancol = (255/255, 127/255, 127/255, 1); #ffb2b2
    my @stdcol = (178/255, 178/255, 255/255, 0.8); #7f7fff
    my @std1col = (0, 0, 0, 0.02); #ff7f7f
    my @std2col = (0, 0, 0, 0.02); #ff7f7f
    my @linecol = (0, 0, 0, 0.4);
    my @helplinecol = (1, 1, 1, 0.9);
    my @background = (0.95, 0.95, 0.95, 1);
    my @tickcol = (0, 0, 0, 0.8);
    my @labelcol = (0, 0, 0, 1);

    #create new image
    my $size   = 6;
    my $offset = 20;
    my $left   = 40;
    my $bottom = 15;
    my $top = 20;
    my $height = 200;
    my $surface = Cairo::ImageSurface->create('argb32', $left+$offset*2+($xmax+$zero)*$size,$bottom+$top+$offset*2+$height); #format, width, height
    my $cr = Cairo::Context->create($surface);

    my ($font_extents,$extents,$fontheight,$fontdescent);

    #background
    $cr->rectangle(0, 0, $left+$offset*2+($xmax+$zero)*$size,$bottom+$offset*2+2*200+20);
    $cr->set_source_rgba(1, 1, 1, 1);
    $cr->fill;

    #fonts
    $cr->select_font_face ('sans', 'normal', 'normal');

    $cr->save;

    #set up work space
    $cr->set_antialias('none');
    $cr->set_line_width(1);

    #background for plot
    $cr->rectangle($left+$offset, $top+$offset, ($xmax+$zero)*$size-1, $height);
    $cr->set_source_rgba(@background);
    $cr->fill;

    #draw legend
    $cr->set_font_size(10);
#    $font_extents = $cr->font_extents;
    my $x = $left+$offset+$size*100-5;
    foreach my $i (reverse (0..scalar(@legend)-1)) {
	$cr->set_antialias('default');
	$x -= $cr->text_extents($legend[$i])->{width};
	$cr->move_to($x,$top+5+9);
	$cr->set_source_rgba(@tickcol);
	$cr->show_text($legend[$i]);
	$x -= 15;
	$cr->set_antialias('none');
	$cr->set_source_rgba(@{$cols[$i]});
	$cr->rectangle($x, $top+5, 10, 10);
	$cr->fill;
	$x -= 15;
    }

    #draw ticks
    $cr->set_antialias('none');
    #x-axis
    $cr->set_source_rgba(@tickcol);
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%5) == 0 && $i > 1 && $i < $xmax) {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+3);
	} else {
	    $cr->move_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height);
	    $cr->line_to($left+$offset+int($size/2)+$size*$i-($zero ? 0 : $size)-1, $top+$offset+$height+1);
	}
    }
    $cr->stroke;

    #y-axis
    $cr->move_to($left+$offset, $top+$offset);
    $cr->line_to($left+$offset-3, $top+$offset);
    $cr->move_to($left+$offset, $top+$offset+$height-1);
    $cr->line_to($left+$offset-3, $top+$offset+$height-1);
    $cr->stroke;

    #helplines
    $cr->set_source_rgba(@helplinecol);
    foreach my $j (1..3) {
        $cr->move_to($left+$offset, $top+$offset+$height*$j/4-($j ? 1 : 0));
        $cr->line_to($left+$offset+($xmax+$zero)*$size, $top+$offset+$height*$j/4-($j ? 1 : 0));
    }
    $cr->stroke;

    $cr->set_antialias('default');

    #tick labels
    $cr->set_source_rgba(@tickcol);
    $cr->set_font_size(10);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};
    #x-axis
    foreach my $i (($zero ? 0 : 1)..$xmax) {
	if(($i%10) == 0 && $i > 1 && $i < $xmax) {
	    $extents = $cr->text_extents($i*$bin);
	    $cr->move_to($left+$offset+int($size/2+1)+$size*$i-($zero ? 0 : $size)-$extents->{width}/2-1-($i == 1 ? 1 : 0), $top+$offset+$height+$fontheight+2);
	    $cr->show_text($i*$bin);
	}
    }
    #y-axis
    $extents = $cr->text_extents(&addCommas($ymax));
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2);
    $cr->show_text(&addCommas($ymax));
    $extents = $cr->text_extents(0);
    $cr->move_to($left+$offset-5-$extents->{width}, $top+$offset+$fontheight/2-2+$height);
    $cr->show_text(0);

    $cr->save;

    #labels
    $cr->set_font_size(14);
    $font_extents = $cr->font_extents;
    $fontheight = $font_extents->{height};

    #axis labels
    $cr->set_source_rgba(@labelcol);
    $extents = $cr->text_extents($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->move_to($left+$offset+($xmax+$zero)*$size/2-$extents->{width}/2, $top+$offset+$height+$fontheight+15);
    $cr->show_text($xlab.($bin>1 ? ' (Bin size: '.$bin.($add ? $add.')' : '') : ''));
    $cr->rotate($PI * 3 / 2);
    $extents = $cr->text_extents($ylab.($bin>1 ? ' (per bin)' : ''));
    $cr->move_to(-($top+$offset+$height/2+$extents->{width}/2+($bin>1 ? 12 : 0)),$offset+10);
    $cr->show_text($ylab.($bin>1 ? ' (per bin)' : ''));

    $cr->restore;

    #draw boxes
    $cr->set_antialias('none');
    foreach my $pos (0..$xmax-($zero ? 0 : 1)) {
	my $tmp = 0;
	foreach my $s (0..$stacks-1) {
	    next unless($matrix->[$s]->[$pos]);
	    my $cur = $matrix->[$s]->[$pos] / $ymax;
	    $cr->set_source_rgba(@{$cols[$s]});
	    if($cur) {
		$cr->rectangle($left+$offset+$pos*$size, $top+$offset+$height-$tmp*$height, $size-1, -$cur*$height);
		$cr->fill;
	    }
	    $tmp += $cur;
	}
    }

    #write image
    $cr->show_page;
    return $surface;
}

sub header {
    return '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
<title>PRINSEQ-'.$WHAT.' Report</title>
<style type="text/css">
<!--/* <![CDATA[ */
'.$CSS_STYLE.'
/* ]]> */-->
</style>
</head>
<body>
<center>';
}

sub footer {
    return '</center></body></html>';
}

sub generateHtml {
    my ($in,$out) = @_;
    my ($file,$data,$surface,$html,$png);
    $data = &readGdFile($in);
    my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));

    $html .= &header();
    $html .= '<h2>PRINSEQ-'.$WHAT.' v'.$VERSION.' HTML Report&nbsp;&nbsp;&nbsp;</h2>[Generated: '.$time.']<br /><br />';
    $html .= '<div class="info-panel">';

    #input info
    if(exists $data->{numseqs}) {
        $html .= '<div class="info-header"><span class="info-header-title">Input Information</span></div>';
        $html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Input file(s):</td><td class="info-table-value">'.($data->{filename1} ? &convertIntToString($data->{filename1}) : '-').($data->{filename2} ? ' and '.&convertIntToString($data->{filename2}) : '').'</td></tr><tr><td class="info-table-type">Input format(s):</td><td class="info-table-value">'.($data->{format1} ? uc($data->{format1}) : '-').($data->{format2} ? ' and '.uc($data->{format2}) : '').'</td></tr>';
        if(exists $data->{pairedend} && $data->{pairedend}) {
            my $singletons1 = ($data->{numseqs}||0)-($data->{pairs}||0);
            my $singletons2 = ($data->{numseqs2}||0)-($data->{pairs}||0);
            $html .= '<tr><td class="info-table-type"># Sequences (file 1):</td><td class="info-table-value">'.&addCommas($data->{numseqs}||'-').'</td></tr><tr><td class="info-table-type">Total bases (file 1):</td><td class="info-table-value">'.&addCommas($data->{numbases}||'-').'</td></tr><tr><td class="info-table-type"># Sequences (file 2):</td><td class="info-table-value">'.&addCommas($data->{numseqs2}||'-').'</td></tr><tr><td class="info-table-type">Total bases (file 2):</td><td class="info-table-value">'.&addCommas($data->{numbases2}||'-').'</td></tr><tr><td class="info-table-type"># Pairs:</td><td class="info-table-value">'.&addCommas($data->{pairs}||'-').($data->{pairs} ? '&nbsp;&nbsp;('.sprintf("%.2f",(100*(2*$data->{pairs})/(($data->{numseqs}||0)+($data->{numseqs2}||0)))).'% of sequences)' : '').'</td></tr></tr><tr><td class="info-table-type"># Singletons (file 1):</td><td class="info-table-value">'.&addCommas($singletons1).($singletons1 ? '&nbsp;&nbsp;('.sprintf("%.2f",(100*$singletons1/$data->{numseqs})).'%)' : '').'</td></tr><tr><td class="info-table-type"># Singletons (file 2):</td><td class="info-table-value">'.&addCommas($singletons2).($singletons2 ? '&nbsp;&nbsp;('.sprintf("%.2f",(100*$singletons2/$data->{numseqs2})).'%)' : '').'</td></tr>';
        } else {
            $html .= '<tr><td class="info-table-type"># Sequences:</td><td class="info-table-value">'.&addCommas($data->{numseqs}||'-').'</td></tr><tr><td class="info-table-type">Total bases:</td><td class="info-table-value">'.&addCommas($data->{numbases}||'-').'</td></tr>';
        }
        $html .= '</tbody></table></div><hr>';
    }

    #length plot
    if(exists $data->{counts}->{length} && keys %{$data->{counts}->{length}}) {
        $html .= '<div class="info-header"><span class="info-header-title">Length Distribution</span></div>';
        if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<div class="info-content"><b>File 1</b><br /><table border="0" cellpadding="0" cellspacing="0"> <tbody><tr><td class="info-table-type">Mean sequence length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{mean} ? sprintf("%.2f",$data->{stats}->{length}->{mean}) : '-').' &plusmn; '.(exists $data->{stats}->{length}->{std} ? sprintf("%.2f",$data->{stats}->{length}->{std}) : '-').' bp</td></tr><tr><td class="info-table-type">Minimum length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{min} ? &addCommas($data->{stats}->{length}->{min}) : '-').' bp</td></tr> <tr><td class="info-table-type">Maximum length:</td><td class="info-table-value">'.(exists $data->{stats}->{length}->{max} ? &addCommas($data->{stats}->{length}->{max}) : '-').' bp</td></tr> <tr><td class="info-table-type">Length range:</td><td class="info-table-value">'.(exists $data->{stats}->{length}->{range} ? &addCommas($data->{stats}->{length}->{range}) : '-').' bp</td></tr> <tr><td class="info-table-type">Mode length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{mode} ? &addCommas($data->{stats}->{length}->{mode}) : '-').' bp with '.(exists $data->{stats}->{length}->{modeval} ? &addCommas($data->{stats}->{length}->{modeval}) : '-').' sequences</td></tr></tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{length},1),$data->{stats}->{length},'Length Distribution','Read Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
            $html .= '<div class="info-content"><br /><b>File 2</b><br /><table border="0" cellpadding="0" cellspacing="0"> <tbody><tr><td class="info-table-type">Mean sequence length:</td> <td class="info-table-value">'.(exists $data->{stats2}->{length}->{mean} ? sprintf("%.2f",$data->{stats2}->{length}->{mean}) : '-').' &plusmn; '.(exists $data->{stats2}->{length}->{std} ? sprintf("%.2f",$data->{stats2}->{length}->{std}) : '-').' bp</td></tr><tr><td class="info-table-type">Minimum length:</td> <td class="info-table-value">'.(exists $data->{stats2}->{length}->{min} ? &addCommas($data->{stats2}->{length}->{min}) : '-').' bp</td></tr> <tr><td class="info-table-type">Maximum length:</td><td class="info-table-value">'.(exists $data->{stats2}->{length}->{max} ? &addCommas($data->{stats2}->{length}->{max}) : '-').' bp</td></tr> <tr><td class="info-table-type">Length range:</td><td class="info-table-value">'.(exists $data->{stats2}->{length}->{range} ? &addCommas($data->{stats2}->{length}->{range}) : '-').' bp</td></tr> <tr><td class="info-table-type">Mode length:</td> <td class="info-table-value">'.(exists $data->{stats2}->{length}->{mode} ? &addCommas($data->{stats2}->{length}->{mode}) : '-').' bp with '.(exists $data->{stats2}->{length}->{modeval} ? &addCommas($data->{stats2}->{length}->{modeval}) : '-').' sequences</td></tr></tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{length},1),$data->{stats2}->{length},'Length Distribution','Read Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
        } else {
            $html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"> <tbody><tr><td class="info-table-type">Mean sequence length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{mean} ? sprintf("%.2f",$data->{stats}->{length}->{mean}) : '-').' &plusmn; '.(exists $data->{stats}->{length}->{std} ? sprintf("%.2f",$data->{stats}->{length}->{std}) : '-').' bp</td></tr><tr><td class="info-table-type">Minimum length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{min} ? &addCommas($data->{stats}->{length}->{min}) : '-').' bp</td></tr> <tr><td class="info-table-type">Maximum length:</td><td class="info-table-value">'.(exists $data->{stats}->{length}->{max} ? &addCommas($data->{stats}->{length}->{max}) : '-').' bp</td></tr> <tr><td class="info-table-type">Length range:</td><td class="info-table-value">'.(exists $data->{stats}->{length}->{range} ? &addCommas($data->{stats}->{length}->{range}) : '-').' bp</td></tr> <tr><td class="info-table-type">Mode length:</td> <td class="info-table-value">'.(exists $data->{stats}->{length}->{mode} ? &addCommas($data->{stats}->{length}->{mode}) : '-').' bp with '.(exists $data->{stats}->{length}->{modeval} ? &addCommas($data->{stats}->{length}->{modeval}) : '-').' sequences</td></tr></tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{length},1),$data->{stats}->{length},'Length Distribution','Read Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
        }
        $html .= '<hr>';
    }

    #GC content
    if(exists $data->{counts}->{gc} && keys %{$data->{counts}->{gc}}) {
	$html .= '<div class="info-header"><span class="info-header-title">GC Content Distribution</span></div>';
        if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<div class="info-content"><b>File 1</b><br /><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Mean GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{mean} ? sprintf("%.2f",$data->{stats}->{gc}->{mean}) : '-').' &plusmn; '.(exists $data->{stats}->{gc}->{std} ? sprintf("%.2f",$data->{stats}->{gc}->{std}) : '-').' %</td></tr> <tr><td class="info-table-type">Minimum GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{min} ? $data->{stats}->{gc}->{min} : '-').' %</td></tr> <tr><td class="info-table-type">Maximum GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{max} ? $data->{stats}->{gc}->{max} : '-').' %</td></tr> <tr><td class="info-table-type">GC content range:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{range} ? $data->{stats}->{gc}->{range} : '-').' %</td></tr> <tr><td class="info-table-type">Mode GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{mode} ? $data->{stats}->{gc}->{mode} : '-').' % with '.(exists $data->{stats}->{gc}->{modeval} ? &addCommas($data->{stats}->{gc}->{modeval}) : '-').' sequences</td></tr> </tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{gc},0),$data->{stats}->{gc},'GC Content Distribution','GC Content (0-100%)','Number of Sequences','',1);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
            $html .= '<div class="info-content"><br /><b>File 2</b><br /><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Mean GC content:</td> <td class="info-table-value">'.(exists $data->{stats2}->{gc}->{mean} ? sprintf("%.2f",$data->{stats2}->{gc}->{mean}) : '-').' &plusmn; '.(exists $data->{stats2}->{gc}->{std} ? sprintf("%.2f",$data->{stats2}->{gc}->{std}) : '-').' %</td></tr> <tr><td class="info-table-type">Minimum GC content:</td> <td class="info-table-value">'.(exists $data->{stats2}->{gc}->{min} ? $data->{stats2}->{gc}->{min} : '-').' %</td></tr> <tr><td class="info-table-type">Maximum GC content:</td> <td class="info-table-value">'.(exists $data->{stats2}->{gc}->{max} ? $data->{stats2}->{gc}->{max} : '-').' %</td></tr> <tr><td class="info-table-type">GC content range:</td> <td class="info-table-value">'.(exists $data->{stats2}->{gc}->{range} ? $data->{stats2}->{gc}->{range} : '-').' %</td></tr> <tr><td class="info-table-type">Mode GC content:</td> <td class="info-table-value">'.(exists $data->{stats2}->{gc}->{mode} ? $data->{stats2}->{gc}->{mode} : '-').' % with '.(exists $data->{stats2}->{gc}->{modeval} ? &addCommas($data->{stats2}->{gc}->{modeval}) : '-').' sequences</td></tr> </tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{gc},0),$data->{stats2}->{gc},'GC Content Distribution','GC Content (0-100%)','Number of Sequences','',1);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
        } else {
            $html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Mean GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{mean} ? sprintf("%.2f",$data->{stats}->{gc}->{mean}) : '-').' &plusmn; '.(exists $data->{stats}->{gc}->{std} ? sprintf("%.2f",$data->{stats}->{gc}->{std}) : '-').' %</td></tr> <tr><td class="info-table-type">Minimum GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{min} ? $data->{stats}->{gc}->{min} : '-').' %</td></tr> <tr><td class="info-table-type">Maximum GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{max} ? $data->{stats}->{gc}->{max} : '-').' %</td></tr> <tr><td class="info-table-type">GC content range:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{range} ? $data->{stats}->{gc}->{range} : '-').' %</td></tr> <tr><td class="info-table-type">Mode GC content:</td> <td class="info-table-value">'.(exists $data->{stats}->{gc}->{mode} ? $data->{stats}->{gc}->{mode} : '-').' % with '.(exists $data->{stats}->{gc}->{modeval} ? &addCommas($data->{stats}->{gc}->{modeval}) : '-').' sequences</td></tr> </tbody></table><br>';
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{gc},0),$data->{stats}->{gc},'GC Content Distribution','GC Content (0-100%)','Number of Sequences','',1);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
            $html .= '</div>';
        }
        $html .= '<hr>';
    }

    #Base quality
    if(exists $data->{quals} || exists $data->{qualsmean} || exists $data->{qualsbin}) {
	$html .= '<div class="info-header"><span class="info-header-title">Base Quality Distribution</span></div><div class="info-content">';
        if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<b>File 1</b><br />';
        }
    }
    if(exists $data->{quals} && keys %{$data->{quals}}) {
	$surface = &createBoxPlot(&convertToBoxValues($data->{quals},4),'Base Quality Distribution','Read position in %','Quality score','');
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= &insert_image($png);
    }
    if(exists $data->{qualsbin} && keys %{$data->{qualsbin}}) {
	$surface = &createBoxPlot(&convertToBoxValues($data->{qualsbin},4),'Base Quality Distribution','Read position in bp','Quality score','',0,'bp',$data->{binval});
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
        $html .= '<br /><br />' if(exists $data->{quals});
	$html .= &insert_image($png);
    }
    if(exists $data->{qualsmean} && keys %{$data->{qualsmean}}) {
	$surface = &createBarPlot(&convertToBarValues($data->{qualsmean},5,1),'Sequence Quality Distribution','Mean of quality scores per sequence','Number of sequences','',0);
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= '<br /><br />' if(exists $data->{qualsbin});
	$html .= &insert_image($png);
    }
    if(exists $data->{pairedend} && $data->{pairedend}) {
        if(exists $data->{quals} || exists $data->{qualsmean} || exists $data->{qualsbin}) {
            $html .= '<br /><br /><br /><b>File 2</b><br />';
        }
        if(exists $data->{quals2} && keys %{$data->{quals2}}) {
            $surface = &createBoxPlot(&convertToBoxValues($data->{quals2},4),'Base Quality Distribution','Read position in %','Quality score','');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= &insert_image($png);
        }
        if(exists $data->{qualsbin2} && keys %{$data->{qualsbin2}}) {
            $surface = &createBoxPlot(&convertToBoxValues($data->{qualsbin2},4),'Base Quality Distribution','Read position in bp','Quality score','',0,'bp',$data->{binval});
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br /><br />' if(exists $data->{quals2});
            $html .= &insert_image($png);
        }
        if(exists $data->{qualsmean2} && keys %{$data->{qualsmean2}}) {
            $surface = &createBarPlot(&convertToBarValues($data->{qualsmean2},5,1),'Sequence Quality Distribution','Mean of quality scores per sequence','Number of sequences','',0);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br /><br />' if(exists $data->{qualsbin2});
            $html .= &insert_image($png);
        }
    }
    if(exists $data->{quals} || exists $data->{qualsmean} || exists $data->{qualsbin}) {
	$html .= '</div><hr>';
    }

    #Ns
    if((exists $data->{counts}->{ns} && keys %{$data->{counts}->{ns}}) || (exists $data->{counts2} && exists $data->{counts2}->{ns} && keys %{$data->{counts2}->{ns}})) {
        $html .= '<div class="info-header"><span class="info-header-title">Occurence of N</span></div><div class="info-content">';
        if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<b>File 1</b><br />';
        }
    }
    if(exists $data->{counts}->{ns} && keys %{$data->{counts}->{ns}}) {
	my $nscount = 0;
	foreach my $n (values %{$data->{counts}->{ns}}) {
	    $nscount += $n;
	}	
	$html .= '<table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Sequences with N:</td> <td class="info-table-value">'.($nscount ? &addCommas($nscount).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$nscount).' %)' : 0).'</td></tr><tr><td class="info-table-type">Max percentage of Ns per sequence:</td> <td class="info-table-value">'.(exists $data->{stats}->{ns}->{max} ? $data->{stats}->{ns}->{max} : 0).' %</td></tr></tbody></table>';
	if($nscount) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{ns},1),undef,'Percentage of N\'s (> 0%)','Percentage of N\'s per Read (1-100%)','# Sequences','',0);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
        }
    }
    if(exists $data->{pairedend} && $data->{pairedend} && exists $data->{counts2}->{ns} && keys %{$data->{counts2}->{ns}}) {
        $html .= '<br /><br /><br /><b>File 2</b><br />';
        my $nscount = 0;
        foreach my $n (values %{$data->{counts2}->{ns}}) {
            $nscount += $n;
        }	
        $html .= '<table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">Sequences with N:</td> <td class="info-table-value">'.($nscount ? &addCommas($nscount).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs2}*$nscount).' %)' : 0).'</td></tr><tr><td class="info-table-type">Max percentage of Ns per sequence:</td> <td class="info-table-value">'.(exists $data->{stats2}->{ns}->{max} ? $data->{stats2}->{ns}->{max} : 0).' %</td></tr></tbody></table>';
        if($nscount) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{ns},1),undef,'Percentage of N\'s (> 0%)','Percentage of N\'s per Read (1-100%)','# Sequences','',0);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
        }
    }
    if((exists $data->{counts}->{ns} && keys %{$data->{counts}->{ns}}) || (exists $data->{counts2} && exists $data->{counts2}->{ns} && keys %{$data->{counts2}->{ns}})) {
        $html .= '</div><hr>';
    }

    #tails
    if(exists $data->{tail} || exists $data->{tail2}) {
        $html .= '<div class="info-header"><span class="info-header-title">Poly-A/T Tails</span></div><div class="info-content">';
    }
    if(exists $data->{tail}) {
	my $tail5count = 0;
	foreach my $n (values %{$data->{counts}->{tail5}}) {
	    $tail5count += $n;
	}
	my $tail3count = 0;
	foreach my $n (values %{$data->{counts}->{tail3}}) {
	    $tail3count += $n;
	}
	if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<b>File 1</b><br />';
        }
	$html .= '<table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type"></td><td class="info-table-value">5\'-end</td> <td class="info-table-value">3\'-end</td></tr> <tr><td class="info-table-type">Sequences with tail:</td><td class="info-table-value">'.($tail5count ? &addCommas($tail5count).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$tail5count).' %)' : 0).'</td> <td class="info-table-value">'.($tail3count ? &addCommas($tail3count).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$tail3count).' %)' : 0).'</td></tr> <tr><td class="info-table-type">Maximum tail length:</td> <td class="info-table-value">'.(exists $data->{stats}->{tail5}->{max} ? $data->{stats}->{tail5}->{max} : 0).'</td> <td class="info-table-value">'.(exists $data->{stats}->{tail3}->{max} ? $data->{stats}->{tail3}->{max} : 0).'</td></tr></tbody></table>';
        if($tail5count) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{tail5},1),undef,'Poly-A/T Tail Distribution (> 4bp)','5\' Tail Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
            if($tail3count) {
                $html .= '<br />';
            }
        }
        if($tail3count) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts}->{tail3},1),undef,'Poly-A/T Tail Distribution (> 4bp)','3\' Tail Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
        }
    }
    if(exists $data->{pairedend} && $data->{pairedend} && exists $data->{tail2}) {
        my $tail5count = 0;
	foreach my $n (values %{$data->{counts2}->{tail5}}) {
	    $tail5count += $n;
	}
	my $tail3count = 0;
	foreach my $n (values %{$data->{counts2}->{tail3}}) {
	    $tail3count += $n;
	}
        $html .= '<br /><br /><br /><b>File 2</b><br />';
	$html .= '<table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type"></td><td class="info-table-value">5\'-end</td> <td class="info-table-value">3\'-end</td></tr> <tr><td class="info-table-type">Sequences with tail:</td><td class="info-table-value">'.($tail5count ? &addCommas($tail5count).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs2}*$tail5count).' %)' : 0).'</td> <td class="info-table-value">'.($tail3count ? &addCommas($tail3count).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs2}*$tail3count).' %)' : 0).'</td></tr> <tr><td class="info-table-type">Maximum tail length:</td> <td class="info-table-value">'.(exists $data->{stats2}->{tail5}->{max} ? $data->{stats2}->{tail5}->{max} : 0).'</td> <td class="info-table-value">'.(exists $data->{stats2}->{tail3}->{max} ? $data->{stats2}->{tail3}->{max} : 0).'</td></tr></tbody></table>';
        if($tail5count) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{tail5},1),undef,'Poly-A/T Tail Distribution (> 4bp)','5\' Tail Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
            if($tail3count) {
                $html .= '<br />';
            }
        }
        if($tail3count) {
            $surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{counts2}->{tail3},1),undef,'Poly-A/T Tail Distribution (> 4bp)','3\' Tail Length in bp','# Sequences','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
        }
    }
    if(exists $data->{tail} || exists $data->{tail2}) {
        $html .= '</div><hr>';
    }
    

    #tag sequence check
    if(exists $data->{freqs} || exists $data->{freqs2}) {
        $html .= '<div class="info-header"><span class="info-header-title">Tag Sequence Check</span></div><div class="info-content">';
    }
    if(exists $data->{freqs}) {
	my $tagmidseq;
	if(exists $data->{tagmidseq}) {
	    $tagmidseq = $data->{tagmidseq};
	    $tagmidseq =~ s/\,/\<br \/\>/g;
	}
	if(exists $data->{pairedend} && $data->{pairedend}) {
            $html .= '<b>File 1</b><br />';
        }
	$html .= '<table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type"></td><td class="info-table-value">5\'-end</td><td class="info-table-value">3\'-end</td></tr><tr><td class="info-table-type">Probability of tag sequence:</td><td class="info-table-value">'.(exists $data->{tagprob}->{5} ? $data->{tagprob}->{5}.' %' : '-').'</td><td class="info-table-value">'.(exists $data->{tagprob}->{3} ? $data->{tagprob}->{3}.' %' : '-').'</td></tr><tr><td class="info-table-type">GSMIDs or RLMIDs:</td><td class="info-table-value">'.(exists $data->{tagmidnum} ? ($data->{tagmidnum} == 0 ? 'none' : ($tagmidseq ? $tagmidseq : $data->{tagmidnum})) : '-').'</td><td class="info-table-value">&nbsp;</td></tr></tbody></table><br>';
	$html .= '<table border="0" cellspacing="0" cellpadding="0"><tr><td>'.&insert_image($FREQCHART_L,undef,undef,1).'</td>';
	foreach my $pos (sort {$a <=> $b} keys %{$data->{freqs}->{5}}) {
	    $html .= '<td valign="bottom" align="center">';
	    foreach my $base (qw(A C G T N)) {
		if($data->{freqs}->{5}->{$pos}->{$base}) {
		    $html .= &insert_image($BASE64_BASES->{$base},$data->{freqs}->{5}->{$pos}->{$base},14,1).'<br />';
		    #'<img height="'.$data->{freqs}->{5}->{pos}->{$base}.'px" border="0" src="'..'" alt="'.$base.'" width="14px" /><br />';
		}
	    }
	    $html .= &insert_image($MMCHART_B2,6,16,1).'</td>';
	}
	$html .= '<td align="center" valign="middle">&nbsp;...&nbsp;</td>';
	foreach my $pos (sort {$a <=> $b} keys %{$data->{freqs}->{3}}) {
	    $html .= '<td valign="bottom" align="center">';
	    foreach my $base (qw(A C G T N)) {
		if($data->{freqs}->{3}->{$pos}->{$base}) {
		    $html .= &insert_image($BASE64_BASES->{$base},$data->{freqs}->{3}->{$pos}->{$base},14,1).'<br />';
		}
	    }
	    $html .= &insert_image($MMCHART_B2,6,16,1).'</td>';
	}
	$html .= '</tr>';
	$html .= '<tr><td>&nbsp;</td>';
	foreach my $num (1,0,0,0,5,0,0,0,0,10,0,0,0,0,15,0,0,0,0,20,0,20,0,0,0,0,15,0,0,0,0,10,0,0,0,0,5,0,0,0,1) {
	    $html .= '<td valign="top" align="center" style="font-size: 10px;margin: 0;">'.($num ? $num : '').'&nbsp;</td>';
	}
	$html .= '</tr><tr><td align="left" valign="middle">&nbsp;</td><td align="center" valign="middle" colspan="41" class="pinfo"><b>Position from Sequence Ends</b></td></tr>';
	$html .= '</table>';
    }
    if(exists $data->{pairedend} && $data->{pairedend} && exists $data->{freqs2}) {
        $html .= '<br /><br /><br /><b>File 2</b><br />';
        $html .= '<table border="0" cellpadding="0" cellspacing="0"> <tbody><tr><td class="info-table-type"></td><td class="info-table-value">5\'-end</td><td class="info-table-value">3\'-end</td></tr><tr><td class="info-table-type">Probability of tag sequence:</td><td class="info-table-value">'.(exists $data->{tagprob2}->{5} ? $data->{tagprob2}->{5}.' %' : '-').'</td><td class="info-table-value">'.(exists $data->{tagprob2}->{3} ? $data->{tagprob2}->{3}.' %' : '-').'</td></tr></tbody></table><br>';
	$html .= '<table border="0" cellspacing="0" cellpadding="0"><tr><td>'.&insert_image($FREQCHART_L,undef,undef,1).'</td>';
	foreach my $pos (sort {$a <=> $b} keys %{$data->{freqs2}->{5}}) {
	    $html .= '<td valign="bottom" align="center">';
	    foreach my $base (qw(A C G T N)) {
		if($data->{freqs2}->{5}->{$pos}->{$base}) {
		    $html .= &insert_image($BASE64_BASES->{$base},$data->{freqs2}->{5}->{$pos}->{$base},14,1).'<br />';
		    #'<img height="'.$data->{freqs}->{5}->{pos}->{$base}.'px" border="0" src="'..'" alt="'.$base.'" width="14px" /><br />';
		}
	    }
	    $html .= &insert_image($MMCHART_B2,6,16,1).'</td>';
	}
	$html .= '<td align="center" valign="middle">&nbsp;...&nbsp;</td>';
	foreach my $pos (sort {$a <=> $b} keys %{$data->{freqs2}->{3}}) {
	    $html .= '<td valign="bottom" align="center">';
	    foreach my $base (qw(A C G T N)) {
		if($data->{freqs2}->{3}->{$pos}->{$base}) {
		    $html .= &insert_image($BASE64_BASES->{$base},$data->{freqs2}->{3}->{$pos}->{$base},14,1).'<br />';
		}
	    }
	    $html .= &insert_image($MMCHART_B2,6,16,1).'</td>';
	}
	$html .= '</tr>';
	$html .= '<tr><td>&nbsp;</td>';
	foreach my $num (1,0,0,0,5,0,0,0,0,10,0,0,0,0,15,0,0,0,0,20,0,20,0,0,0,0,15,0,0,0,0,10,0,0,0,0,5,0,0,0,1) {
	    $html .= '<td valign="top" align="center" style="font-size: 10px;margin: 0;">'.($num ? $num : '').'&nbsp;</td>';
	}
	$html .= '</tr><tr><td align="left" valign="middle">&nbsp;</td><td align="center" valign="middle" colspan="41" class="pinfo"><b>Position from Sequence Ends</b></td></tr>';
	$html .= '</table>';
    }
    if(exists $data->{freqs} || exists $data->{freqs2}) {
        $html .= '</div><hr>';
    }

    #Sequence duplicates
    if(exists $data->{dubslength} || exists $data->{dubscounts}) {
	$html .= '<div class="info-header"><span class="info-header-title">Sequence Duplication</span></div>';
    }
    my %dubs;
    if(exists $data->{dubscounts} && keys %{$data->{dubscounts}}) {
	my $exactonly = $data->{exactonly}||0;
	foreach my $n (keys %{$data->{dubscounts}}) {
	    foreach my $s (keys %{$data->{dubscounts}->{$n}}) {
		$dubs{$s}->{count} += $data->{dubscounts}->{$n}->{$s} * $n;
		$dubs{$s}->{max} = $n unless(exists $dubs{$s}->{max} && $dubs{$s}->{max} > $n);
		$dubs{all} += $data->{dubscounts}->{$n}->{$s} * $n;
	    }
	}
	$html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type"></td><td class="info-table-value"># Sequences</td> <td class="info-table-value">Max duplicates</td></tr><tr><td class="info-table-type">Exact duplicates:</td><td class="info-table-value">'.(exists $dubs{0}->{count} ? &addCommas($dubs{0}->{count}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{0}->{count}).' %)' : 0).'</td><td class="info-table-value">'.($dubs{0}->{max}||0).'</td></tr><tr><td class="info-table-type">Exact duplicates with reverse complements:</td><td class="info-table-value">'.(exists $dubs{3}->{count} ? &addCommas($dubs{3}->{count}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{3}->{count}).' %)' : 0).'</td> <td class="info-table-value">'.($dubs{3}->{max}||0).'</td></tr>';
        unless($exactonly) {
            $html .= '<tr><td class="info-table-type">5\' duplicates</td><td class="info-table-value">'.(exists $dubs{1}->{count} ? &addCommas($dubs{1}->{count}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{1}->{count}).' %)' : 0).'</td> <td class="info-table-value">'.($dubs{1}->{max}||0).'</td></tr><tr><td class="info-table-type">3\' duplicates</td><td class="info-table-value">'.(exists $dubs{2}->{count} ? &addCommas($dubs{2}->{count}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{2}->{count}).' %)' : 0).'</td> <td class="info-table-value">'.($dubs{2}->{max}||0).'</td></tr><tr><td class="info-table-type">5\'/3\' duplicates with reverse complements</td><td class="info-table-value">'.(exists $dubs{4}->{count} ? &addCommas($dubs{4}->{count}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{4}->{count}).' %)' : 0).'</td> <td class="info-table-value">'.($dubs{4}->{max}||0).'</td></tr>';
        }
        $html .= '<tr><td class="info-table-type">Total:</td><td class="info-table-value">'.(exists $dubs{all} ? &addCommas($dubs{all}).' &nbsp;('.sprintf("%.2f",100/$data->{numseqs}*$dubs{all}).' %)' : 0).'</td><td class="info-table-value">-</td></tr></tbody></table>';
    }
    if(exists $dubs{all} && $dubs{all}) {
        if(exists $data->{dubslength} && keys %{$data->{dubslength}}) {
            $surface = &createStackBarPlot(&convertOdToStackBinMatrix($data->{dubslength},5,1),'Sequence duplication level','Read Length in bp','Number of duplicates','',0,' bp');
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br />'.&insert_image($png);
        }
        if(exists $data->{dubscounts} && keys %{$data->{dubscounts}}) {
            $surface = &createStackBarPlot(&convertOdToStackBinMatrix($data->{dubscounts},5,1,100),'Sequence duplication level','Number of duplicates','Number of sequences','',0);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br /><br />' if(exists $data->{dubslength});
            $html .= &insert_image($png);
            my %dubsmax;
            my $count = 1;
            foreach my $n (sort {$b <=> $a} keys %{$data->{dubscounts}}) {
                foreach my $s (keys %{$data->{dubscounts}->{$n}}) {
                    foreach my $i (1..$data->{dubscounts}->{$n}->{$s}) {
                        $dubsmax{$count++}->{$s} = $n;
                        last unless($count <= 100);
                    }
                    last unless($count <= 100);
                }
                last unless($count <= 100);
            }
            $surface = &createStackBarPlot(&convertOdToStackBinMatrix(\%dubsmax,5,1,100),'Sequence duplication level','Sequence','Number of duplicates','',0);
            $png = '';
            $surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
            $html .= '<br /><br />' if(exists $data->{dubslength});
            $html .= &insert_image($png);
        }
    }
    if(exists $data->{dubslength} || exists $data->{dubscounts}) {
	$html .= '</div><hr>';
    }

    #Sequence complexity
    if(exists $data->{compldust} || exists $data->{complentropy}) {
	$html .= '<div class="info-header"><span class="info-header-title">Sequence Complexity</span></div>';
	if(exists $data->{complvals}) {
	    my $complseq;
	    foreach my $d (keys %{$data->{complvals}}) {
		foreach my $m ('minseq','maxseq') {
		    $complseq = $data->{complvals}->{$d}->{$m};
		    $complseq = substr($complseq,0,797).'...' if(length($complseq) > 800);
		    $complseq =~ s/(.{60})/$1\<br \/\>/g;
		    $data->{complvals}->{$d}->{$m} = $complseq;
		}
	    }
	}
	$html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type"></td><td class="info-table-value">Value</td><td class="info-table-value">Sequence</td></tr><tr><td class="info-table-type">Minimum DUST score:</td><td class="info-table-value">'.(exists $data->{complvals}->{dust}->{minval} ? $data->{complvals}->{dust}->{minval} : '-').'</td><td class="info-table-value sequencetext">'.(exists $data->{complvals}->{dust}->{minseq} ? $data->{complvals}->{dust}->{minseq} : '').'</td></tr><tr><td class="info-table-type">Maximum DUST score:</td><td class="info-table-value">'.(exists $data->{complvals}->{dust}->{maxval} ? $data->{complvals}->{dust}->{maxval} : '').'</td><td class="info-table-value sequencetext">'.(exists $data->{complvals}->{dust}->{maxseq} ? $data->{complvals}->{dust}->{maxseq} : '').'</td></tr><tr><td class="info-table-type">Minimum Entropy value:</td><td class="info-table-value">'.(exists $data->{complvals}->{entropy}->{minval} ? $data->{complvals}->{entropy}->{minval} : '').'</td><td class="info-table-value sequencetext">'.(exists $data->{complvals}->{entropy}->{minseq} ? $data->{complvals}->{entropy}->{minseq} : '').'</td></tr><tr><td class="info-table-type">Maximum Entropy value:</td><td class="info-table-value">'.(exists $data->{complvals}->{entropy}->{maxval} ? $data->{complvals}->{entropy}->{maxval} : '').'</td><td class="info-table-value sequencetext">'.(exists $data->{complvals}->{entropy}->{maxseq} ? $data->{complvals}->{entropy}->{maxseq} : '').'</td></tr> </tbody></table><br />';
    }
    if(exists $data->{compldust}) {
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{compldust},0),undef,'Sequence complexity distribution','Mean sequence complexity (DUST scores)','Number of sequences','',1);
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= &insert_image($png);
    }
    if(exists $data->{complentropy}) {
	$surface = &createAnnotBarPlot(&convertOdToBinMatrix($data->{complentropy},0),undef,'Sequence complexity distribution','Mean sequence complexity (Entropy values)','Number of sequences','',1);
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= '<br /><br />' if(exists $data->{compldust});
	$html .= &insert_image($png);
    }
    if(exists $data->{compldust} || exists $data->{complentropy}) {
	$html .= '</div><hr>';
    }

    #Dinucleotide odd ratio PCA - microbial/viral
    if(exists $data->{dinucodds} && keys %{$data->{dinucodds}}) {
	$html .= '<div class="info-header"><span class="info-header-title">Dinucleotide Odds Ratios</span></div>';
	$html .= '<div class="info-content"><table border="0" cellpadding="0" cellspacing="0"><tbody><tr><td class="info-table-type">&nbsp;</td>';
	foreach my $d (map {join("/",(m/../g ))} sort keys %{$data->{dinucodds}}) {
	    $html .= '<td class="info-table-value">'.$d.'</td>';
	}
	$html .= '</tr><tr><td class="info-table-type">Odds ratio</td>';
	foreach my $d (map {sprintf("%.4f",$data->{dinucodds}->{$_})} sort keys %{$data->{dinucodds}}) {
	    $html .= '<td class="info-table-value">'.$d.'</td>';
	}
	$html .= '</tr></tbody></table><br />';
	my @new = map {$data->{dinucodds}->{$_}} sort keys %{$data->{dinucodds}};
        $surface = &createOddsRatioPlot($data->{dinucodds},'Odds ratios','Dinucleotide','Odds ratio','');
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= &insert_image($png);
	$surface = &createPCAPlot(&convertToPCAValues(\@new,'m'),'PCA','1st Principal Component Score','2nd Principal Component Score','');
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= '<br /><br />';
	$html .= &insert_image($png);
	$surface = &createPCAPlot(&convertToPCAValues(\@new,'v'),'PCA','1st Principal Component Score','2nd Principal Component Score','');
	$png = '';
	$surface->write_to_png_stream(sub { my ($closure, $data) = @_; $png .= $data; });
	$html .= '<br /><br />';
	$html .= &insert_image($png);
	$html .= '</div>';
    }

    $html .= '</div>';
    $html .= &footer();

    #write html to file
    $file = &getFileName('.html');
    open(FH, ">$file") or &printError("Can't open file ".$file.": $!");
    print FH $html;
    close(FH);
    &printLog("Done with HTML data");
}

sub insert_image {
  my ($data, $height, $width, $noencode) = @_;
  my $content .= '<img border="0" '.($height ? 'height="'.$height.'"' : '').($width ? 'width="'.$width.'"' : '').' src="'.($noencode ? "data:image/png;base64,".$data : &inline_image($data)).'" />'."\n";
  return $content;
}

sub inline_image {
    return "data:image/png;base64,".MIME::Base64::encode_base64($_[0]);
}

sub convertIntToString {
    my $int = shift;
    $int =~ s/(.{2})/chr(hex($1))/eg;
    return $int;
}
