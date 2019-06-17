#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


RESULTSPATHLOCAL=RESULTSPATH<>"2level/";
FILEPREFIX="2level-";
Print["Saving output to "<>RESULTSPATHLOCAL<>FILEPREFIX<>"..."]


graphlist=Flatten@Table[All2LevelGraphs[i,Range[2,i]],{i,5,32}]//Map[Function[gd,Module[{
kSys=VertexCount@gd["graph"],
graphs
},
graphs=EnvironmentPlot[#,kSys]&/@AllHairGraphs[gd["graph"],Subsets[#,{1}]&];
Thread[{kSys,graphs}]
]]]//Flatten[#,1]&//SortBy[First];
Print[ToString[Length@graphlist]<>" graphs precomputed"];


(* multi-resolution samples, concentration point where we "expect" a crossing roughly *)
samples=<|
"depolarizing"->MultiResolutionSamples["depolarizing",0.1904],
"BB84"->MultiResolutionSamples["BB84",0.1106],
"2Pauli"->MultiResolutionSamples["2Pauli",0.1135]
|>;
evaluatedSamples=Join@@KeyValueMap[Function[{channel,pts},
CHANNELS[channel]/.p->#&/@pts
],samples];
Print[ToString[Length@evaluatedSamples]<>" sample points per graph"];


graphlist//Map@Apply[Function[{kSys,G}, Module[{
filename=RESULTSPATHLOCAL<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[kSys],
prettyG=EnvironmentPlot[G,kSys],
pNumeric=Numeric@@evaluatedSamples,
pa
},
If[FileExistsQ[filename],Echo["skipping "<>filename];Return[0,Module]];
Echo["running "<>filename];

pa=Catch@PauliActionCC[G,kSys,pNumeric(* numeric poly solver only works with the symmetric solver *)];
Echo["done PA for "<>filename];
If[StringQ[pa],Echo[pa];Return[]]; (* display errors *)

With[{
pdep=CIThreshold[pa,pNumeric,;;Length@samples["depolarizing"]],
pbb84=CIThreshold[pa,pNumeric,Length@samples["depolarizing"]+1;;Length@samples["depolarizing"]+Length@samples["BB84"]],
p2pauli=CIThreshold[pa,pNumeric,-Length@samples["2Pauli"];;]
},
 <|
"graph"->G,
"kSys"->kSys,
"pauli action"->pa,
"samples"->samples,
"CI"->CIMult[pa],
"threshold dep"->pdep,
"threshold bb84"->pbb84,
"threshold 2pauli"->p2pauli
|>//Put[#,filename]&;

Echo["done "<>filename<>", pdep="<>ToString@pdep<>", pbb84="<>ToString@pbb84<>", p2pauli="<>ToString@p2pauli];
]]]];
