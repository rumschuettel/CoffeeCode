(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


LaunchKernels[16];


(* ::Input::Initialization:: *)
Clear[AllGraphs,AllColorings]
AllGraphs[kSys_,kEnv_]:=With[{
cmd="~/nauty27rc2/geng -qc "<>ToString[kSys]
},
Print[cmd];
Import["!"<>cmd,"List"]
]

AllColorings[graphStr_String,kEnv_]:=With[{
cmd="~/nauty27rc2/vcolg -Tq -m"<>ToString[kEnv+1]
},
Print[cmd];
Import["!echo "<>graphStr<>" | "<>cmd, "Table"]
]


(* ::Input::Initialization:: *)
Table[Module[{
graphs=AllGraphs[kSys,kEnv]
},

Print["found "<>ToString@Length@graphs<>" graphs for kSys="<>ToString[kSys]];

graphs//ParallelMap[Function[graphStr,Module[{
allColorings,
pruned,
graph=ImportString[graphStr,"Graph6"],
fn
},

fn=RESULTSPATH<>"/all-env/all-env-"<>GraphHash[graph,kSys]<>"."<>ToString@kSys<>".adjm";
If[FileExistsQ[fn],
Print["skipping "<>graphStr],
allColorings=AllColorings[graphStr,kEnv];
Print["found "<>ToString@Length@allColorings<>" colorings for "<>graphStr];


pruned=allColorings//Map[Function[out,Module[{
nv=out[[1]],
ne=out[[2]],
coloring=out[[3;;3+kSys-1]],
edges=out[[3+kSys;;]],
envEdges,
g
},

If[Max[coloring]==0,Nothing,(* no environment vertex *)

(* edges to env *)
envEdges=MapIndexed[Function[{envTarget,idx},
If[envTarget==0,Nothing,
(First@idx-1)\[UndirectedEdge](envTarget+kSys)
]
],coloring];


g=Graph[
Join[
UndirectedEdge@@@Partition[edges,2],
envEdges
]];

(* destroy any legacy information from graph, and canonicalize *)
Graph@EdgeList@CanonicalizeByOrbits[AdjacencyGraph@AdjacencyMatrix@g,kSys]

]
]],#]&//Union;

Print["pruned to "<>ToString@Length@pruned<>" colorings for "<>graphStr];

Print["exporting "<>fn];
Export[fn, {kSys, VertexCount@#-kSys, StringTrim@ExportString[#,"Graph6"],ExportAdjacencyMatrix[#]}&/@pruned,"Table"];

]; (* file existed *)
]], #, Method->"FinestGrained"]&;

],
{kSys,2,8},{kEnv,kSys,kSys}
];

Print["done"];



