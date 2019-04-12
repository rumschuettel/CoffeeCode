#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


RESULTSPATHLOCAL=RESULTSPATH<>"357codes/";
FILEPREFIX="357code-";
Print["Saving output to "<>RESULTSPATHLOCAL<>FILEPREFIX<>"..."]


graphlist=Flatten@Table[
ConcatGraph[StarGraph[i],StarGraph[j],Range@k,1],
{i,{3,5,7}},{j,{3,5,7}},{k,1,i}
];
With[{
G=HairGraph[#,{1}],
kSys=VertexCount@#
},
Module[{
filename=RESULTSPATHLOCAL<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[kSys],
prettyG=EnvironmentPlot[G,kSys],
pa
},
If[FileExistsQ[filename],Echo["skipping "<>filename];Return[Module]];
Echo["running "<>filename];

pa=Catch@PauliActionCC[G,kSys,p,False(*symmetric solver*)];
Echo["done PA for "<>filename];
If[StringQ[pa],Echo[pa];Return[]]; (* handle errors *)

With[{
pdep=p/.CIThreshold[pa,p]
},
 <|
"graph"->G,
"kSys"->kSys,
"pauli action"->Compress[pa],
"threshold dep"->pdep
|>//Put[#,filename]&;

Echo["done "<>filename<>"\np="<>ToString@pdep];
]

]]&/@graphlist;
