(* ::Package:: *)

(*tscale=10;*)


(*tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;
split=1;
(*splitTimes=Partition[times,steps/split];*)
splitTimes=Split[times,!Or@@Table[#1<=m tmax/split<=#2,{m,split-1}]&];*)


tmax=20;
steps=2000;
dt=tmax/steps;
times=Range[0+dt,tmax,dt];
split=1;
splitTimes=Partition[times,steps/split];


(*tminExp=-2;
tmaxExp=2;
tmax=10.^tmaxExp;
steps=1000;
tExps=Range[tminExp,tmaxExp,(tmaxExp-tminExp)/(steps-1)];
times=10.^#&/@tExps;*)


runs=100;


(* ::Subsubsection:: *)
(*vars*)


length=20;


sites=length;


suLocalDim=2;


maxGroupSize=1;


dis=6;


localPot=RandomReal[{-dis,dis},sites];


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


diffs=Table[1/Abs[(localPot[[n]]-localPot[[n+1]])],{n,length-1}];
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
clustSites=SortBy[clustSites,First];


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


bonds=Table[{ii,ii+1},{ii,sites-1}];


extBonds=Complement[bonds,Flatten[Tuples[#,2]&/@clustSites,1]];


intBonds=Complement[bonds,extBonds];


(*clustFromSite[ss_]:=Quotient[ss-1,clustSize]+1*)


(*clustBonds=Map[clustFromSite,extBonds,{2}];*)


(*bonds={{1,2}};*)


(* ::Subsubsection:: *)
(*Matrices*)


localMat=PauliMatrix[3];


hopMats={PauliMatrix[1],PauliMatrix[2],PauliMatrix[3]};


(* ::Subsubsection:: *)
(*Ham*)


clustOp[op_,ss_]:=KroneckerProduct[IdentityMatrix[suLocalDim^(s2p[ss]-1)],op,IdentityMatrix[suLocalDim^(s2l[ss]-s2p[ss])]]


(*siteOp[op_,ss_]:=clustOp[op,siteToClustPos[ss]]*)


(*disConst=4*dis;*)


(*crosCoup[s1_,s2_]:=Abs[s1-s2]^-\[Alpha]+(length-Abs[s1-s2])^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=Min[Abs[s1-s2],(length-Abs[s1-s2])]^-\[Alpha]*)


(*crosCoup[s1_,s2_]:=If[Abs[s1-s2]\[Equal]1,1,0]*)


crosCoup[s1_,s2_]:=1


(*localPot=Flatten[disConst*corrands+harmV];*)


(*localPot=Array[h,length];*)


(*potHam=Table[Sum[localPot[[containedSite[cc,ss]]]
clustOp[numOp,cc,ss]
,{ss,clustSize}],{cc,numClust}];*)


potHam=Table[Sum[localPot[[ss]]
clustOp[localMat,ss]
,{ss,cc}],{cc,clustSites}];


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


initKetSite=Riffle[{0,1}&/@Range[sites/2],{1,0}&/@Range[sites/2]];


initKet=Table[
ketList=List/@initKetSite[[cc]];
If[Length[ketList]>1,
Flatten[KroneckerProduct@@ketList,1],
Flatten[ketList,2]
]
,{cc,clustSites}];


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


observables=allObs;
obsfun=Function[{values},
(*avNum=Flatten[(#.values&/@multBy)\[Transpose],1];*)
sigZ=matsToVars[PauliMatrix[3],values];
spinSum=Table[(-1)^n,{n,sites}].sigZ;
{sigZ,spinSum}
];
