(* ::Package:: *)

(* ::Subsection:: *)
(*import*)


done[rname_,runs_]:=(
fullList={};
Table[If[FileExistsQ["/data/shainen/"<>rname<>"/r"<>ToString[kk]<>"/dataTWA.dat"],AppendTo[fullList,kk]],{kk,0,runs-1}];
notDone=Complement[Range[1,runs],fullList];
Save["doneThings.dat",notDone];
)


dir=StringSplit[ParentDirectory[],"/"][[5]];


done[dir,101]
