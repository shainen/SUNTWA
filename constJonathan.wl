(* ::Package:: *)

(*tscale=10;*)


(*tminExp=-2;
tmaxExp=4;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;
split=100;
(*splitTimes=Partition[times,steps/split];*)
splitTimes=Split[times,!Or@@Table[#1<=m tmax/split<=#2,{m,split-1}]&];*)


tmax=20;
steps=1000;
dt=tmax/steps;
times=N[Range[0,tmax-dt,dt]];
split=1;
splitTimes=Partition[times,steps/split];


(*tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;*)


runs=1;


(* ::Subsubsection:: *)
(*vars*)


jcoup=1/Sqrt[3];
delta=1/Sqrt[3];


length=6;


sites=length;


suLocalDim=2;


maxGroupSize=2;


(*dis=10;*)


(*localPot=RandomReal[{-dis,dis},sites];*)


(*localPot={-0.495460963974,0.229981111878,-0.707592970791,-0.669264172184,-0.309590954355,1.09176755111,-0.747575234176,0.0972739232244,-1.43377636215,0.273191494177,-1.83010670053,0.65801897364,0.512345988527,0.122983644592,-1.44248842988,0.780950867261,-0.0876437995615,-1.55356320629,-0.458311511781,-0.244775319171,-0.255088527936,1.29166945103,1.41098713916,1.12901903104,1.27575216702,-0.998962634945,-1.41634540577,1.29083781132,-0.171546775935,-1.13372878017,-1.17806459389,0.128384470004,0.0449326143957,0.608721596483,1.31593606743,-0.314752808226,-1.52345102415,0.847472424367,-0.312480647777,-0.471583309149};*)


localPot=delta*{1.23860, 1.55957,-0.02625, 0.28480,-.78682,2.59600};


(*diffs=Table[1/Abs[(localPot[[n]]-localPot[[n+1]])],{n,length-1}];
bsites={};
While[diffs!=Table[0,{length-1}],
max=Position[diffs,Max[diffs]][[1,1]];
AppendTo[bsites,max];
diffs[[max]]=0;
If[max!=length-1,diffs[[max+1]]=0];
If[max!=1,diffs[[max-1]]=0];
]*)


(*ii=1;
clustSites={};
Table[
AppendTo[clustSites,Range[ii,ii+jj-1]];
ii+=jj;
,{jj,clustSizes}];*)


(*ii=1;
clustSites={};
While[ii<=sites,
If[MemberQ[bsites,ii],
AppendTo[clustSites,Range[ii,ii+1]];
ii+=2;,
AppendTo[clustSites,Range[ii,ii]];
ii+=1;
]
]*)


(*clustSizes={1,2,1};*)


(*diffs=Table[1/Abs[(localPot[[n]]-localPot[[n+1]])],{n,length-1}];
clustSites={};
While[diffs!=Table[0,{length-1}],
max=Position[diffs,Max[diffs]][[1,1]];
g1=Select[clustSites,MemberQ[#,max]&];
g2=Select[clustSites,MemberQ[#,max+1]&];
newGroup=Sort[DeleteDuplicates[Flatten[{max,max+1,g1,g2}]]];
If[Length[newGroup]<=maxGroupSize,
clustSites=Complement[clustSites,g1];
clustSites=Complement[clustSites,g2];
AppendTo[clustSites,newGroup];
];
diffs[[max]]=0;
]
If[!MemberQ[Flatten[clustSites],#],AppendTo[clustSites,{#}]]&/@Range[sites];
clustSites=SortBy[clustSites,First];*)


clustSites=Partition[Range[sites],maxGroupSize];


clustSizes=Length/@clustSites;


suClustDims=suLocalDim^clustSizes;


s2c[ss_]:=Position[clustSites,ss][[1,1]]


s2p[ss_]:=Position[clustSites,ss][[1,2]]


s2l[ss_]:=clustSizes[[s2c[ss]]]


numClust=Length[clustSites];


(*containedSite[clust_,num_]:=clustSize(clust-1)+num*)


(*containedSites[clust_]:=Table[containedSite[clust,n],{n,clustSize}]*)


(*ordPairs=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDim}],{ii,suClustDim}],1];*)


realPairs[clustNum_]:=Flatten[Table[Table[{ii,jj},{jj,ii,suClustDims[[clustNum]]}],{ii,suClustDims[[clustNum]]}],1];


imPairs[clustNum_]:=Flatten[Table[Table[{ii,jj},{jj,ii+1,suClustDims[[clustNum]]}],{ii,suClustDims[[clustNum]]}],1];


coToLiR[cc_,co_]:=Position[realPairs[cc],co]


coToLiI[cc_,co_]:=Position[imPairs[cc],co]


(*bonds=Table[{n,Mod[n+1,length,1]},{n,length-1}];*)


(*siteToClustNum[site_]:=Quotient[site,clustSize,1]+1*)


(*siteToClustPos[site_]:=Mod[site,clustSize,1]*)


(*addl[num_]:=Mod[num,length]*)


(*nfc[coord_]:=FromDigits[addl[coord],length]+1*)


(*bondsHoriz=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{0,1}]},{xx,0,length-1},{yy,0,length-2}],1];*)


(*bondsVert=Flatten[Table[{nfc[{xx,yy}],nfc[{xx,yy}+{1,0}]},{xx,0,length-2},{yy,0,length-1}],1];*)


bonds=Table[{ii,Mod[ii+1,sites,1]},{ii,sites}];


extBonds=Complement[bonds,Flatten[Tuples[#,2]&/@clustSites,1]];


intBonds=Complement[bonds,extBonds];


(*clustFromSite[ss_]:=Quotient[ss-1,clustSize]+1*)


(*clustBonds=Map[clustFromSite,extBonds,{2}];*)


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Matrices*)


localMat=PauliMatrix[3];


hopMats={PauliMatrix[1],PauliMatrix[2],Sqrt[0.1]*PauliMatrix[3]};


(* ::Subsubsection:: *)
(*Ham*)


clustOp[op_,ss_]:=KroneckerProduct[IdentityMatrix[suLocalDim^(s2p[ss]-1)],op,IdentityMatrix[suLocalDim^(s2l[ss]-s2p[ss])]]


(*siteOp[op_,ss_]:=clustOp[op,siteToClustPos[ss]]*)


(*disConst=4*dis;*)


(*crosCoup[s1_,s2_]:=Abs[s1-s2]^-\[Alpha]+(length-Abs[s1-s2])^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=Min[Abs[s1-s2],(length-Abs[s1-s2])]^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=If[Abs[s1-s2]\[Equal]1,1,0]*)


crosCoup[s1_,s2_]:=jcoup


(*localPot=Flatten[disConst*corrands+harmV];*)


(*localPot=Array[h,length];*)


(*potHam=Table[Sum[localPot[[containedSite[cc,ss]]]
clustOp[numOp,cc,ss]
,{ss,clustSize}],{cc,numClust}];*)


potHam=Table[Sum[localPot[[ss]]
clustOp[localMat,ss]
,{ss,cc}],{cc,clustSites}];


(*potHam=0;*)


(*selfCoup=Table[
Sum[
If[MemberQ[intBonds,{containedSite[cc,s1],containedSite[cc,s2]}],
crosCoup[s1,s2](clustOp[sX,s1].clustOp[sX,s2]+clustOp[sY,s1].clustOp[sY,s2]),
0
]
,{s1,clustSize},{s2,clustSize}],{cc,numClust}];*)


selfCoup=Table[
Sum[
s1=pp[[1]];
s2=pp[[2]];
If[MemberQ[intBonds,pp],
crosCoup[s1,s2]Sum[clustOp[mats,s1].clustOp[mats,s2],{mats,hopMats}],
0
]
,{pp,Tuples[cc,2]}],{cc,clustSites}];


localHam=potHam+selfCoup;


(*crosHam1=PauliMatrix[1];
crosHam2=PauliMatrix[2];*)


(* ::Input:: *)
(*(*crosHam={clustOp[PauliMatrix[#],clustSize]&/@Range[3],clustOp[PauliMatrix[#],1]&/@Range[3]};*)*)


crosHamFunc[site_]:=Table[clustOp[mats,site],{mats,hopMats}]


(*initKetSite=Riffle[{0,1}&/@Range[sites/2],{1,0}&/@Range[sites/2]];*)


(*initKet=Table[
ketList=List/@initKetSite[[cc]];
If[Length[ketList]>1,
Flatten[KroneckerProduct@@ketList,1],
Flatten[ketList,2]
]
,{cc,clustSites}];*)


(*initKet={{-2.810657460437625543*10^(-01)+7.015253844296336083*10^(-02)I,  7.852878689403355872*0.01+-6.192005576846262294*0.1I}, {-3.648629763895638867*0.1+-5.324742131074827745*0.1I,  2.572602126900596850*0.01+2.307222487417846801*0.1I}, {5.902328944158318214*0.1+-3.067415813462650021*0.1I,  -4.490762813675505116*0.1+-7.282147742146949376*0.01I}, {1.994789645845160830*0.1+-1.310334457212311288*0.1I,  -4.919171265126216497*0.1+3.279853370808094026*0.1I}}\[Transpose];*)


initKet={{-2.426353554217702824*0.1+-2.648005089919371335*0.1I,7.951051132552851286*0.01+-7.373185215600709663*0.01I,-3.887888245264342069*0.1+-1.396259490061620978*0.1I},{3.107066623255042362*0.01+3.961251073785703736*0.1I,4.271567656597460849*0.1+1.591008888680371658*0.1I,-6.694773248319145775*0.1+-2.732814324808480189*0.01I},{1.722040133043985233*0.1+-3.294586310365372195*0.1I,-4.980668084377977145*0.1+2.964585210922828784*0.1I,7.475278316445083115*0.01+-3.222444608459837934*0.1I},{-7.478410150428786984*0.1+1.251594670736896342*0.1I,-3.853562971838829787*0.1+5.440661918406550779*0.1I,3.439578536556726074*0.1+-3.907232424387242498*0.1I}}\[Transpose];


(* ::Subsubsection:: *)
(*Observable*)


(*siteObs=PauliMatrix[3]*)


(*basicObs=Table[cR[ss][#,#],{ss,numClust}]&/@Range[suClustDim];*)


allObs=Table[{cR[cc]@@@realPairs[cc],cI[cc]@@@imPairs[cc]},{cc,numClust}];


(*multBy=Table[Diagonal[clustOp[numOp,n]],{n,clustSize}];*)


(*multByCross=Table[Diagonal[clustOp[PauliMatrix[3],n[[1]]].clustOp[PauliMatrix[3],n[[2]]]],{n,Flatten[Table[Table[{ii, jj}, {jj, ii+1, clustSize}], {ii, clustSize}], 1]}];*)


locMatToVars[cc_,mat_, clVars_] := Total[(If[Length[coToLiR[cc,#1]] == 1, 2 Re[#2] clVars[[1, coToLiR[cc,#1][[1, 1]]]], 0]/(1 + KroneckerDelta[#1[[1]], #1[[2]]]) + If[Length[coToLiI[cc,#1]] == 1, -2 Im[#2] clVars[[2, coToLiI[cc,#1][[1, 1]]]], 0]) & @@@ Drop[ArrayRules[mat], -1]]


matsToVars[mat_, clVars_]:=
Flatten[
Table[
Table[
locMatToVars[cc,clustOp[mat,ss],clVars[[cc]]]
,{ss,clustSites[[cc]]}]
,{cc,numClust}]
]


clustMatsToVars[mat_, clVars_]:=
Flatten[
Table[
locMatToVars[cc,mat,clVars[[cc]]]
,{cc,numClust}]
]


(*observables=allObs;
obsfun=Function[{values},
(*avNum=Flatten[(#.values&/@multBy)\[Transpose],1];*)
sigZ=matsToVars[PauliMatrix[3],values];
spinSum=Table[(-1)^n,{n,sites}].sigZ;
{sigZ,spinSum}
];*)


orderedP=PauliMatrix/@{3,1,2,0};


basisMats=Flatten[Outer[KroneckerProduct,orderedP,orderedP,1],1];


observables=allObs;
obsfun=Function[{values},
jonBasis=Table[clustMatsToVars[basisMats[[n]],values],{n,16}];
jonBasis
];
