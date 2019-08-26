#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


With[{
kSysMax=25
},

graphsToTest=ParallelTable[
With[{
graphs=All2LevelGraphs[kSys]
},
Thread[{VertexCount@#["graph"],#["legcount"],AllHairGraphs[#["graph"],{{1}}&]}]&/@graphs
],
{kSys,3,kSysMax}
]//Flatten[#,2]&//SortBy[First]
];

Print["Testing "<>ToString@Length@graphsToTest<>" graphs"]


Clear[collectData2Level]
collectData2Level[kSys_Integer,legcount_Integer,G_Graph]:=Module[{
FILEPREFIX="2level.root.",
filename,lockfilename,
jsonString
},
filename=RESULTSPATH<>"2level/"<>FILEPREFIX<>GraphHash[G,kSys]<>"-legcount-"<>ToString[legcount]<>"-kSys"<>ToString[kSys]<>"-kTot"<>ToString[kSys+1]<>".json";
lockfilename=filename<>".lock";

If[FileExistsQ[filename]\[Or]FileExistsQ[lockfilename],Print["Skipping "<>filename];Return[]];

Export[lockfilename,$MachineName,"Text"];
Pause[.5];

Print["Computing "<>ToString[kSys]<>"-graph code..."];
jsonString=If[VertexCount@G<15,
RunCC[G,kSys,"Multi","Full"],
RunCC[G,kSys,"Multi"]
];

Export[filename,jsonString,"RawJSON","Compact"->True];
DeleteFile[lockfilename];
Print["Done "<>filename];
];


collectData2Level@@@graphsToTest
