(* ::Package:: *)

(* ::Subsection:: *)
(*import*)


done[rname_,runs_]:=(
fullList={};
Table[If[FileExistsQ["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/dataTWA.dat"],AppendTo[fullList,kk]],{kk,0,runs-1}];
fullList
)


dir=StringSplit[ParentDirectory[],"/"][[5]];


done[dir,100];
