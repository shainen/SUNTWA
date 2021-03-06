(* ::Package:: *)

(* ::Section:: *)
(*All*)


(* ::Subsubsection:: *)
(*setup*)


(*SetDirectory[NotebookDirectory[]]*)


SetDirectory[Directory[]<>"/SUNTWA"];


SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> 1];
SetSystemOptions["MKLThreads" -> 1];


<<randomSeed.wl


<<dynSUNRI.wl


<<ndsolve.wl


<<inits.wl


(*<<gausInits.wl*)


<<constRandHeis.wl


(*<<constJonathan.wl*)


(*<<constRepSU3.wl*)


(*<<constHeisAllToAllSpin1.wl*)


(* ::Subsection:: *)
(*run TWA*)


Dynamic[rr]


(*Timing[start=makeDSolveStart[localHam,crosHamFunc,observables];]*)


(*Timing[eachTWA={};
Table[
AppendTo[eachTWA,singleRun[start,Flatten[gausInitsOR],obsfun]];
,{rr,runs}];
fullTWA=Total[eachTWA]/runs;
varTWA=Variance[eachTWA];]*)


(*Timing[
fullTWA=0;
squares=0;
Table[
newObs=singleRun[start,Flatten[discInitsOR],obsfun];
AddTo[fullTWA,newObs/runs];
AddTo[squares,newObs^2/runs];
,{rr,runs}];]*)


Timing[
fullTWA=0;
squares=0;
firstTime=First@splitTimes;
nextTimes=Drop[splitTimes,1];
Table[
<<constRandHeis.wl;
start=makeDSolveStart[localHam,crosHamFunc,observables];
stuff=singleRunShort[start,discInitsOR,firstTime];
(*stuff=singleRunShort[start,meanInitsOR,firstTime];*)
lastTime=Last@firstTime;
Table[
stuff=Join[stuff,singleRunShort[start,Flatten[discInitsMid[lastTime,Last@stuff]],trange]];
lastTime=Last@trange;
,{trange,nextTimes}];
newObs=Chop[obsfun/@stuff];
AddTo[fullTWA,newObs/runs];
AddTo[squares,newObs^2/runs];
,{rr,runs}];]


(*Timing[fullTWA=0;
Table[
<<constRandHeisLR.wl;
start=makeDSolveStart[localHam,crosHamFunc,observables];
AddTo[fullTWA,singleRun[start,Flatten[discInitsOR],obsfun]/runs];
,{rr,runs}];]*)


mmu=MaxMemoryUsed[]/10.^6;


SetDirectory[ParentDirectory[]];


allData=fullTWA[[All,2]];
spins=fullTWA[[All,1]];
(*stError=Sqrt[squares-fullTWA^2]/Sqrt[runs];*)


Save["dataTWA.dat",{mmu,allData,spins,squares,runs,localPot,clustSites}];


(*Save["jons6JinfFIXED.dat",fullTWA];*)
