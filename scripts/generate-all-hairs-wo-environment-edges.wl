(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


(* ::Input::Initialization:: *)
Clear[AllGraphs,AllColorings]
AllGraphs[kSys_,kEnv_]:=With[{
cmd="~/nauty27rc2/geng -qc "<>ToString[kSys]
},
Print[cmd];
Import["!"<>cmd,"List"]
]

AllColorings[graphStr_String,kEnv_]:=With[{
cmd="~/nauty27rc2/vcolg -qT -m"<>ToString[kEnv+1]
},
Print[cmd];
RunProcess["cmd",All,graphStr]["StandardOutput"]
]


(* ::Input::Initialization:: *)
Table[Module[{
graphs=AllGraphs[kSys,kEnv]
},

Print["running kSys="<>ToString[kSys]<>", kEnv="<>ToString[kEnv]];

graphs//Map[Function[graph,With[{
allColorings=ImportString[
AllColorings[graph,kEnv],
"Table"
]
},

allColorings//Map[Function[out,Module[{
nv=out[[1]],
ne=out[[2]],
coloring=out[[3;;3+kSys-1]],
edges=out[[3+kSys;;]],
envEdges
},

(* edges to env *)
envEdges=MapIndexed[Function[{envTarget,idx},
If[envTarget==0,Nothing,
(First@idx-1)\[UndirectedEdge](envTarget+kSys)
]
],coloring];


ExportAdjacencyMatrix@CanonicalizeByOrbits[Graph[
Range[0,kSys+kEnv-1],
Join[
UndirectedEdge@@@Partition[edges,2],
envEdges
]],
kSys
]
]]]//Union
]]]

],
{kSys,2,4},{kEnv,kSys,kSys}
]



