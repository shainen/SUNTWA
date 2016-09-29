(* ::Package:: *)

(* ::Subsection:: *)
(*import*)


(*import[rname_,runs_,ignore_,comb_]:=(
list=Partition[Complement[Range[0,runs-1],ignore],comb];
eachOne={};
eachVarImb={};
Do[
tempAll=0;
tempVar=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/spinchain.dat"];
data={TWASpDiscSU4,TWASpGauSU4,TWASpDelta};
AddTo[tempAll,data];
AddTo[tempVar,data[[All,3]]-data[[All,2]]^2];
,{kk,rr}];
AppendTo[eachOne,tempAll/comb];
AppendTo[eachVarImb,tempVar/comb];
,{rr,list}];
all=Total[eachOne]/Length[list];
avImb=all[[All,2]];
varImb=Total[eachVarImb]/Length[list];
)*)


(*import[rname_,runs_,ignore_,comb_]:=(
list=Partition[Complement[Range[0,runs-1],ignore],comb];
eachOne={};
Do[
tempAll=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/dataTWA.dat"];
data=allData;
AddTo[tempAll,data];
,{kk,rr}];
AppendTo[eachOne,tempAll/comb];
,{rr,list}];
mean=Total[eachOne]/Length[list];
stEr=Sqrt[Total[eachOne^2]/Length[list]-Total[eachOne/Length[list]]^2];
)*)


import[rname_,runs_,names_,ignore_]:=(
fullList={};
Table[If[FileExistsQ[ParentDirectory[]<>"/"<>rname<>".qsub.o"<>ToString[qname]<>"-"<>ToString[ll]],AppendTo[fullList,ll]],{ll,0,runs-1},{qname,names}];
fullList=Complement[fullList,ignore];
comb=Quotient[Length[fullList],10];
list=Partition[fullList,comb];
eachOne={};
Do[
tempAll=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/dataTWA.dat"];
data=allData;
AddTo[tempAll,data];
,{kk,rr}];
AppendTo[eachOne,tempAll/comb];
,{rr,list}];
mean=Total[eachOne]/Length[list];
stEr=Sqrt[Total[eachOne^2]/Length[list]-Total[eachOne/Length[list]]^2];
)


(*importE[rname_,bunch_,ignore_]:=(
list=Complement[Range[0,bunch-1],ignore];
tempAll=0;
tempSq=0;
numPoints=0;
Do[
Get["/data/shainen/"<>rname<>"/r"<>ToString[rr]<>"/dataTWA.dat"];
AddTo[tempAll,runs*allData];
AddTo[tempSq,runs*squares];
AddTo[numPoints,runs];
,{rr,list}];
mean=tempAll/numPoints;
stEr=Sqrt[tempSq/numPoints-mean^2]/Sqrt[numPoints];
)*)


dir=StringSplit[ParentDirectory[],"/"][[5]];


import[dir,100,{323525},{}];


(*importE[dir,100,{}];*)


(*Save["/data/shainen/"<>dir<>"compiled.dat",{avg,eachOne}];*)


Save["/data/shainen/"<>dir<>"_compiled.dat",{mean,stEr}];
