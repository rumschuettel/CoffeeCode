#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


RESULTSPATHLOCAL=RESULTSPATH<>"polys/";
FILEPREFIX="polys-";
Print["Saving output to "<>RESULTSPATHLOCAL<>FILEPREFIX<>"..."]
inputGraph=StarGraph[10];(*Import[InputString[],"Graph6"];*)
kSys=VertexCount@inputGraph;
Print["read from STDIN: "<>ToString@inputGraph];


graphlist=AllHairGraphs[inputGraph]
Print[ToString[Length@graphlist]<>" graphs precomputed"];


graphlist//Map@Apply[Function[G, Module[{
filename=RESULTSPATHLOCAL<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[kSys],
prettyG=EnvironmentPlot[G,kSys],
pExact=Exact[q0,q1,q2,q3],
pa
},
If[FileExistsQ[filename],Echo["skipping "<>filename];Return[0,Module]];
Echo["running "<>filename];

pa=Catch@PauliActionCC[G,kSys,pExact,False (*use full solver*)];
Echo["done PA for "<>filename];
If[StringQ[pa],Echo[pa];Return[]]; (* display errors *)

 <|
"graph"->G,
"kSys"->kSys,
"pauli action"->pa,
"CI"->CIMult[pa]
|>//Put[#,filename]&;

Echo["done "<>filename];
]]];
