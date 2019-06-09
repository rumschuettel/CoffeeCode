#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


RESULTSPATHLOCAL=RESULTSPATH<>"polys/";
FILEPREFIX="polys-";
Print["Saving output to "<>RESULTSPATHLOCAL<>FILEPREFIX<>"..."]
inputGraph=ImportString[InputString[],"Graph6"];
kSys=VertexCount@inputGraph;
Print["read from STDIN: "<>ToString@inputGraph];


graphlist=AllHairGraphs[inputGraph]
Print[ToString[Length@graphlist]<>" graphs precomputed"];


graphlist//Map[Function[G, Module[{
mat=ExportAdjacencyMatrix[G]
},
Print[CCRELEASEPATHFULLSOLVERARBCH<>"CoffeeCode."<>ToString@kSys<>"."<>ToString[VertexCount@G-kSys]<>" "<>mat];
]]];
