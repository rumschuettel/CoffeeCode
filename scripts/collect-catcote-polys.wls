#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]
Quiet@Needs["CCInterface`",CurrentDir<>"CCInterface.m"]


Clear[concatRepCode]
concatRepCode[n1_Integer,n2_Integer,keepEnvVertices_:False]:=Module[{roots,connectors,edgeSet},
(* n1: inner code, n2: outer code *)
(* resulting state will have n1(n2+1) qubits: n1*n2 system qubits and n1 environment qubits *)
roots=Table[n1*(i-1)+1,{i,2,n2+1}];
connectors=Range[n1];
edgeSet=Join[Flatten[Table[UndirectedEdge[#,#+i]&/@roots,{i,n1-1}]],Flatten[Table[UndirectedEdge[roots[[i]],#]&/@connectors,{i,Length[roots]}]]];
If[keepEnvVertices,
Return[Graph[edgeSet]];,
Return[Subgraph[Graph[edgeSet],Range[n1*n2+1]]];
];
];


Clear[collectDataRep]
collectDataRep[pair_,writeToFile_:True]:=Module[{G,kSys,FILEPREFIX="concatcodes.",filename,n1=pair[[1]],n2=pair[[2]],jsonString},
G=concatRepCode[n1,n2];
kSys=n1*n2;
filename=RESULTSPATH<>"/concat_codes/"<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[n1]<>"-in-"<>ToString[n2];
Print["Computing "<>ToString[n1]<>"-in-"<>ToString[n2]<>" code..."];
jsonString=If[VertexCount@G<15,
RunCC[G,kSys,"Multi","Full"],
RunCC[G,kSys,"Multi"]
];
If[writeToFile,
Export[filename,jsonString,"RawJSON","Compact"->True]
];
];


Clear[collectDataGraph]
collectDataGraph[G_Graph,kSys_Integer,writeToFile_:True]:=Module[{FILEPREFIX="graphcodes.",filename,jsonString},
filename=RESULTSPATH<>"/graph_codes/"<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[kSys];
Print["Computing "<>ToString[kSys]<>"-graph code..."];
jsonString=If[VertexCount@G<15,
RunCC[G,kSys,"Multi","Full"],
RunCC[G,kSys,"Multi"]
];
If[writeToFile,
Export[filename,jsonString,"RawJSON","Compact"->True]
];
];


singleReps=Table[{1,i},{i,1,60}];
concatPairs=Flatten[Table[{i,j},{i,2,7},{j,2,7}],1];
Join[singleReps,concatPairs]//Map[collectDataRep];


G13=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,1\[UndirectedEdge]12,1\[UndirectedEdge]13,1\[UndirectedEdge]14}];
G15=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,11\[UndirectedEdge]12,11\[UndirectedEdge]13,1\[UndirectedEdge]14,1\[UndirectedEdge]15,1\[UndirectedEdge]16}];
G16=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,11\[UndirectedEdge]12,11\[UndirectedEdge]13,1\[UndirectedEdge]14,1\[UndirectedEdge]15,1\[UndirectedEdge]16,1\[UndirectedEdge]17}];
G17=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,11\[UndirectedEdge]12,11\[UndirectedEdge]13,1\[UndirectedEdge]14,14\[UndirectedEdge]15,14\[UndirectedEdge]16,1\[UndirectedEdge]17,1\[UndirectedEdge]18}];
G18=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,11\[UndirectedEdge]12,11\[UndirectedEdge]13,1\[UndirectedEdge]14,14\[UndirectedEdge]15,14\[UndirectedEdge]16,1\[UndirectedEdge]17,1\[UndirectedEdge]18,1\[UndirectedEdge]19}];
G18b=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,2\[UndirectedEdge]5,2\[UndirectedEdge]6,1\[UndirectedEdge]7,7\[UndirectedEdge]8,7\[UndirectedEdge]9,1\[UndirectedEdge]10,10\[UndirectedEdge]11,10\[UndirectedEdge]12,1\[UndirectedEdge]13,13\[UndirectedEdge]14,13\[UndirectedEdge]15,1\[UndirectedEdge]16,1\[UndirectedEdge]17,1\[UndirectedEdge]18,1\[UndirectedEdge]19}];
G19=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,1\[UndirectedEdge]5,5\[UndirectedEdge]6,5\[UndirectedEdge]7,1\[UndirectedEdge]8,8\[UndirectedEdge]9,8\[UndirectedEdge]10,1\[UndirectedEdge]11,11\[UndirectedEdge]12,11\[UndirectedEdge]13,1\[UndirectedEdge]14,14\[UndirectedEdge]15,14\[UndirectedEdge]16,1\[UndirectedEdge]17,1\[UndirectedEdge]18,1\[UndirectedEdge]19,1\[UndirectedEdge]20}];
G19b=Graph[{1\[UndirectedEdge]2,2\[UndirectedEdge]3,2\[UndirectedEdge]4,2\[UndirectedEdge]5,2\[UndirectedEdge]6,1\[UndirectedEdge]7,7\[UndirectedEdge]8,7\[UndirectedEdge]9,1\[UndirectedEdge]10,10\[UndirectedEdge]11,10\[UndirectedEdge]12,1\[UndirectedEdge]13,13\[UndirectedEdge]14,13\[UndirectedEdge]15,1\[UndirectedEdge]16,1\[UndirectedEdge]17,1\[UndirectedEdge]18,1\[UndirectedEdge]19,1\[UndirectedEdge]20}];
G21=HairGraph[RootedGraphProduct[StarGraph[7],StarGraph[3],Range@7,1],{1}];
{G13,G15,G16,G17,G18,G18b,G19,G19b,G21}//Map[Function[graph,collectDataGraph[graph,VertexCount[graph]-1];]];
