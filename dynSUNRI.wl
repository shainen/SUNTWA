(* ::Package:: *)

vR[ss_][xx_,yy_]:=Which[yy<xx,cR[ss][yy,xx][t],True,cR[ss][xx,yy][t]]


vI[ss_][xx_,yy_]:=Which[xx==yy,0,yy<xx,-cI[ss][yy,xx][t],True,cI[ss][xx,yy][t]]


(*dotO=Function[{ss,a,b,sudim},
SparseArray[{{i_,a}->-1.0*I vO[ss][i,b]},{sudim,sudim}]+
SparseArray[{{b,i_}->1.0*I vO[ss][a,i]},{sudim,sudim}]
{eqsR,eqsI}
];*)


(*dotO=Function[{ss,a,b,sudim},
SparseArray[{{i_,a}->-1.0*I vO[ss][i,b]},{sudim,sudim}]+
SparseArray[{{b,i_}->1.0*I vO[ss][a,i]},{sudim,sudim}]
{eqsR,eqsI}
];*)


dotR=Function[{ss,a,b,sudim},
eqsR=SparseArray[{{b,i_}->-1/2 vI[ss][a,i]},{sudim,sudim}]+
SparseArray[{{a,i_}->-1/2 vI[ss][b,i]},{sudim,sudim}]+
SparseArray[{{i_,b}->-1/2 vI[ss][a,i]},{sudim,sudim}]+
SparseArray[{{i_,a}->-1/2 vI[ss][b,i]},{sudim,sudim}];
eqsI=SparseArray[{{b,i_}->1/2 vR[ss][a,i]},{sudim,sudim}]+
SparseArray[{{i_,b}->-1/2 vR[ss][a,i]},{sudim,sudim}]+
SparseArray[{{a,i_}->1/2 vR[ss][b,i]},{sudim,sudim}]+
SparseArray[{{i_,a}->-1/2 vR[ss][b,i]},{sudim,sudim}];
{eqsR,eqsI}
];


dotI=Function[{ss,a,b,sudim},
eqsR=SparseArray[{{b,i_}->1/2 vR[ss][a,i]},{sudim,sudim}]+
SparseArray[{{i_,b}->1/2 vR[ss][a,i]},{sudim,sudim}]+
SparseArray[{{a,i_}->-1/2 vR[ss][b,i]},{sudim,sudim}]+
SparseArray[{{i_,a}->-1/2 vR[ss][b,i]},{sudim,sudim}];
eqsI=SparseArray[{{b,i_}->1/2 vI[ss][a,i]},{sudim,sudim}]+
SparseArray[{{i_,a}->1/2 vI[ss][b,i]},{sudim,sudim}]+
SparseArray[{{i_,b}->-1/2 vI[ss][a,i]},{sudim,sudim}]+
SparseArray[{{a,i_}->-1/2 vI[ss][b,i]},{sudim,sudim}];
{eqsR,eqsI}
];


dotForHam[ss_, ham_] :=Total[(Re[#2] dotR[ss, #1[[1]], #1[[2]], Length[ham]]+Im[#2] dotI[ss, #1[[1]], #1[[2]], Length[ham]]) & @@@ Drop[ArrayRules[ham], -1]]


varForHam[ss_,ham_]:=Total[(Re[#2] vR[ss][ #1[[1]], #1[[2]]]+Im[#2] vI[ss][#1[[1]], #1[[2]]])&@@@Drop[ArrayRules[ham],-1]]


crossForHam[a_,b_,cHam1_,cHam2_]:=dotForHam[a,cHam1]varForHam[b,cHam2]


makeDSolveStart = Function[{localHam, crosHamFunc, observables},
   eqsRI = Table[
     matLocal = dotForHam[cc, localHam[[cc]]];
     matCross=
     Sum[
       Sum[
         otherSite=(Complement[pairs, {inSites}][[1]]);
         crosCoup[inSites,otherSite]Sum[
                                      crossForHam[cc,s2c[otherSite],crosHamFunc[inSites][[n]], crosHamFunc[otherSite][[n]]]
                                    ,{n,Length[crosHamFunc[inSites]]}]
       ,{pairs,Cases[extBonds,{_,inSites}|{inSites,_}]}]
     ,{inSites, clustSites[[cc]]}];
     matEq = matLocal + matCross;
     {cR[cc][#1, #2]'[t] == matEq[[1,#1, #2]] & @@@ realPairs[cc] , cI[cc][#1, #2]'[t] == matEq[[2,#1, #2]] & @@@ imPairs[cc]}
     , {cc, numClust}];
   initsRI = Table[{cR[cc][#1, #2][0] == 0 & @@@ realPairs[cc] , cI[cc][#1, #2][0] == 0 & @@@ imPairs[cc]}, {cc, numClust}];
                        First@NDSolve`ProcessEquations[Flatten[{eqsRI, initsRI}], observables, t, Method -> {"EquationSimplification" -> "Solve"}]
                        ];
