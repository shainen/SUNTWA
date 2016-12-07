(* ::Package:: *)

matR[clustNum_,i_,j_]:=(SparseArray[{{i,j}->1},{suClustDims[[clustNum]],suClustDims[[clustNum]]}]+SparseArray[{{j,i}->1},{suClustDims[[clustNum]],suClustDims[[clustNum]]}])/2


matI[clustNum_,i_,j_]:=(SparseArray[{{i,j}->-I},{suClustDims[[clustNum]],suClustDims[[clustNum]]}]+SparseArray[{{j,i}->I},{suClustDims[[clustNum]],suClustDims[[clustNum]]}])/2


randFromMat[mat_,init_]:=RandomChoice[Abs[init.#]^2&/@(Eigensystem[N[mat]][[2]])->(Eigensystem[N[mat]][[1]])]


meanFromMat[mat_,init_]:=init\[Conjugate].mat.init


(* ::Subsubsection:: *)
(*Real*)


discRandOR[s_,i_,j_]:=randFromMat[matR[s,i,j],initKet[[s]]]


discRandOI[s_,i_,j_]:=randFromMat[matI[s,i,j],initKet[[s]]]


discInitsOR := Table[{cR[ss][#1, #2][0] == discRandOR[ss,#1,#2] & @@@ realPairs[ss],cI[ss][#1, #2][0] == discRandOI[ss,#1,#2] & @@@ imPairs[ss]}, {ss, numClust}]


meanInitsOR := Table[{cR[ss][#1, #2][0] == meanFromMat[matR[ss,#1,#2],initKet[[ss]]] & @@@ realPairs[ss],cI[ss][#1, #2][0] == meanFromMat[matI[ss,#1,#2],initKet[[ss]]] & @@@ imPairs[ss]}, {ss, numClust}]


discInitsMid[st_,obs_] := Table[{cR[ss][#1, #2][st] == obs[[ss,1,coToLiR[ss,{#1,#2}][[1,1]]]] & @@@ realPairs[ss],cI[ss][#1, #2][st] ==obs[[ss,2,coToLiI[ss,{#1,#2}][[1,1]]]] & @@@ imPairs[ss]}, {ss, numClust}]
