#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]


(* ::Input::Initialization:: *)
(* fix paths in this file *)
Get["CCInterfacePaths.m",Path->CurrentDir];


(* ::Input::Initialization:: *)
Assert[On];(* enable sanity checks instead of failing silently *)


(* ::Input::Initialization:: *)
(* avoid race conditions when building in parallel *)
Clear[ParallelBuildDirectory]
ParallelBuildDirectory[]:=Module[{
targetFolder=PARALLELBUILDPATH<>"release-"<>ToString@$SessionID<>"."<>ToString@$KernelID<>"/"
},
If[!DirectoryQ[targetFolder],
Echo["setting up parallel build directory: "<>targetFolder];
RunProcess[{PARALLELBUILDSETUPSCRIPT,targetFolder}] (* bug with commands like "source" *)
];
targetFolder
]


(* ::Input::Initialization:: *)
Clear[BuglessRunProcess]
BuglessRunProcess[exec_String,what_,input_String:"",opts:OptionsPattern[]]:=Module[{
process=StartProcess[exec,opts],
out
},
If[input!="",WriteLine[process,input]];
While[ProcessStatus[process,"Running"],Pause[1]];
out=ReadString[process];
<|
"ExitCode"->ProcessInformation[process]["ExitCode"],
"StandardError"->"",
"StandardOutput"->out
|>
]

Unprotect[StringForm];
StringForm[a_,b_Parallel`Kernels`kernel]:=ToString@a



(* ::Input::Initialization:: *)
Needs["IGraphM`"]
ParallelNeeds["IGraphM`"];


(* ::Input::Initialization:: *)
Clear[GroupForGraph]
GroupForGraph[graph_Graph,kSys_Integer,extraStabilizedVertices_List:{}]:=GroupForGraph[graph,kSys,extraStabilizedVertices]=With[{
kTot=VertexCount@graph,
(* {1,3,11} for a graph of size 14 would be mapped to {1,0,2,0,0,0,0,0,0,0,3,0,0,0} *)
extraColors=Normal@SparseArray[Thread[extraStabilizedVertices->Range@Length@extraStabilizedVertices],VertexCount@graph]
},
PermutationGroup[
PermutationCycles/@IGBlissAutomorphismGroup[{graph,"VertexColors"->ConstantArray[0,kSys]~Join~ConstantArray[kTot^2,kTot-kSys]+extraColors}]
]
]


(* ::Input::Initialization:: *)
Clear[Subgroup]
Subgroup[group_PermutationGroup,kSys_Integer]:=Subgroup[group,kSys]=group/.{p_Integer/;p>kSys->Nothing}/.Cycles[{}]->Nothing


(* ::Input::Initialization:: *)
Clear[EnvironmentPlot]
EnvironmentPlot[g_Graph,systemSize_]:=Graph[g,VertexStyle->(
Thread[(Range[VertexCount@g-systemSize]+systemSize)->Red]
),VertexSize->.2,VertexLabels->None]


(* ::Input::Initialization:: *)
Clear[IsomorphicHairChoiceQ,HairColoring]
HairColoring[hairs_]:=HairColoring[Sort@hairs]=Length/@GroupBy[hairs,Identity]
IsomorphicHairChoiceQ[g_Graph,hairs1_,hairs2_]:=IsomorphicHairChoiceQ[g,hairs1,hairs2]=IsomorphicHairChoiceQ[g,hairs2,hairs1]=IGBlissIsomorphicQ[
{g,"VertexColors"->HairColoring@hairs1},
{g,"VertexColors"->HairColoring@hairs2}
]
IsomorphicHairChoiceQ[g_Graph]:=IsomorphicHairChoiceQ[g,#1,#2]&


(* ::Input::Initialization:: *)
Clear[HairGraph,AllHairGraphs]
HairGraph[g_Graph,hairs_]:=HairGraph[g,hairs]=With[{
new=VertexCount@g+Range@Length@hairs
},
EdgeList[g]~Join~Thread[hairs\[UndirectedEdge]new]//Graph
]
AllHairGraphs[g_Graph,generator_:Hold[Subsets[#][[2;;]]&]]:=AllHairGraphs[g,generator]=With[{
(* first check for graph isomorphism using a vertex coloring *)
uniqueHairEndpoints=DeleteDuplicates[
ReleaseHold[generator][VertexList@g],
IsomorphicHairChoiceQ[g]
]
},
EnvironmentPlot[HairGraph[g,#],VertexCount@g]&/@uniqueHairEndpoints
]


(* ::Input::Initialization:: *)
Clear[SGSTransversal]
(* the head element of the group stabilizer chain already delivers a transversal of a strong generating set *)
SGSTransversal[group_PermutationGroup]:=With[{
chain=GroupStabilizerChain[group]
},
Partition[chain,2,1]/.{Rule[sA_,gA_],Rule[sB_,gB_]}:>(
Complement[sB,sA]->PermutationGroup@Complement[GroupGenerators@gA,GroupGenerators@gB]
)/.{Rule[s_,PermutationGroup[{}]]:>Nothing}
]


(* ::Input::Initialization:: *)
Clear[ProductOfSymmetricGroupsQ]
ProductOfSymmetricGroupsQ[group_PermutationGroup]:=With[{
orbits=GroupOrbits@group
},
(GroupOrder/@(
GroupStabilizer[group,Complement[Flatten@orbits,#]]&/@orbits
)
)==(Factorial/@Length/@orbits)
]



(* ::Input::Initialization:: *)
Clear[CanonicalizeByOrbits]
(* sort graph vertices such that environment all the way to the right, and orbits of size > 1 come up first *)
CanonicalizeByOrbits[graph_Graph,kSys_Integer,extraStabilizedVertices_List:{}]:=With[{
orbits=GroupOrbits[GroupForGraph[graph,kSys,extraStabilizedVertices],Range@VertexCount@graph],
kTot=VertexCount@graph
},

With[{
coloring=Thread/@MapIndexed[Function[{orbit,i},
(* environment last, length 1 next, nontrivial orbits first; stable sort so add i *)
If[Min@orbit>kSys,
orbit->10kTot+First@i,
If[Length@orbit==1,
orbit->5kTot+First@i,
orbit->First@i
]]
],orbits]//Flatten//SortBy[First]
},
Graph[
IGBlissCanonicalGraph[{graph,"VertexColors"->Last/@coloring}],
VertexLabels->"Index"
]
]
]


(* ::Input::Initialization:: *)
Clear[OrbitsHaveMultipleHairsQ]
OrbitsHaveMultipleHairsQ[group_PermutationGroup,g_Graph,kSys_Integer]:=Module[{
environmentVertices=VertexList[g][[kSys+1;;]],
hairNeighbourIndices
},
hairNeighbourIndices=VertexIndex[g,#]&/@VertexList@Subgraph[NeighborhoodGraph[g,environmentVertices],Complement[VertexList@g,environmentVertices]];

Intersection[#,hairNeighbourIndices]&/@GroupOrbits[group]//AnyTrue[#,Length@##>1&]&
]


(* ::Input::Initialization:: *)
Clear[ExportAdjacencyMatrix]
ExportAdjacencyMatrix[g_Graph]:=With[{
Avec=g//AdjacencyMatrix//Normal//Flatten
},
StringJoin[ToString/@Avec]
]


(* ::Input::Initialization:: *)
Clear[SGSTransversalMMForm]
SGSTransversalMMForm[list_List,kSys_Integer]:=list/.{
Rule[{idx_Integer},PermutationGroup[cycles_List]]:>Module[{
permutationIndices=PermutationList[#,kSys]&/@cycles
},

(* the pivot point is the first index for which the SGSGenerator acts nontrivially on (so one point past the points that are stabilized).
If the pivot point for the SGSGenerator is 2, that means we want to check indices 0 and 1 in C++; in MM this would correspond to indices {1, 2}, which is what \[LeftDoubleBracket];;2\[RightDoubleBracket] truncates to.
Therefore we validly don't subtract 1 from the pivot point, but do subtract 1 from all other indices *)
StringRiffle[Flatten@{
"\tSGSGenerator<"<>ToString@(#1)<>", Group<",
("\t\tPermutation<"<>StringRiffle[ToString/@(##-1),","]<>">")&/@#2,
"\t>>"
},"\n"]&@@{idx,permutationIndices}

]
}//Flatten//("using sgs = SGSTransversal<\n"<>StringRiffle[#,",\n"]<>"\n>;")&

(* for the trivial transversal we simply pass an orbit, sort the graph canonically etc *)
Clear[TrivialSGSTransversalMMForm]
TrivialSGSTransversalMMForm[list_List,kSys_Integer]:=With[{
orbits=list/.{n_Integer}:>Nothing//SortBy[First]
},
With[{
lengths=Accumulate[Length/@orbits]
},
With[{
orbitBounds={0}~Join~lengths//Partition[#,2,1]&
},

StringRiffle[(
"\tTrivialSGSOrbit<"<>ToString[#1]<>", "<>ToString[#2]<>">"
)&@@@orbitBounds,
",\n"]
]//("using sgs = TrivialSGSTransversal<"<>ToString[kSys]<>",\n"<>#<>"\n>;")&
]
]


(* ::Input::Initialization:: *)
Clear[ExportSymmetricCCInstance]
ExportSymmetricCCInstance[graph_Graph,kSys_Integer,extraStabilizedVertices_List:{},name_:"graphstate_instance"]:=With[{
group=Subgroup[GroupForGraph[graph,kSys,extraStabilizedVertices],kSys],
kTot=VertexCount@graph
},

If[!OrbitsHaveMultipleHairsQ[group,graph,kSys]\[And]ProductOfSymmetricGroupsQ@group,
(* write Orbits as expected for TrivialSGS Solver *)
With[{
orbits=GroupOrbits@group,
graphCanonicalized=CanonicalizeByOrbits[graph,kSys,extraStabilizedVertices]
},

If[Length@orbits==0\[Or]Max[Length/@orbits]==1,Echo["WARN: symmetry group is 0. Not a valid symmetric problem instance."];Throw["invalid instance"]
];
(*Echo["Product symmetry group with single hairs detected. Compiling for trivial canonical image provider."];*)

"struct "<>name<>" {\n"<>
TrivialSGSTransversalMMForm[orbits,kSys]<>"\n"<>
"constexpr static size_t k_sys = "<>ToString[kSys]<>", k_env = "<>ToString[kTot-kSys]<>";\n"<>
"constexpr static AdjacencyMatrixT<"<>ToString[kTot]<>"> adjacency_matrix{"<>ToString[Normal@AdjacencyMatrix@graphCanonicalized]<>"};\n"<>
"};"
]
,
(* write SGSGenerators as expected for Nauty *)
With[{
sgs=SGSTransversal@group
},

If[Length@sgs==0,Echo["WARN: symmetry group is 0. Not a valid symmetric problem instance."];
Throw["invalid instance"]
];
Echo["Nested symmetry group or multiple hairs per orbit detected. Compiling for nauty's canonical image provider."];

"struct "<>name<>" {\n"<>
SGSTransversalMMForm[sgs,kSys]<>"\n"<>
"constexpr static size_t k_sys = "<>ToString[kSys]<>", k_env = "<>ToString[kTot-kSys]<>";\n"<>
"constexpr static AdjacencyMatrixT<"<>ToString[kTot]<>"> adjacency_matrix{"<>ToString[Normal@AdjacencyMatrix@graph]<>"};\n"<>
"};"
]
]
]


(* ::Input::Initialization:: *)
Clear[MakeCC]
MakeCC[kSys_Integer,kEnv_Integer]:=With[{
out=RunProcess[
{
"make",
"-B",
"K_SYS="<>ToString@kSys,
"K_ENV="<>ToString@kEnv
},
ProcessDirectory->ParallelBuildDirectory[]
]
},
If[StringTrim@out["StandardError"]!="",Echo[out["StandardError"]]];
If[out["ExitCode"]!=0(*\[Or]out["StandardError"]\[NotEqual]""*),Echo[out];Throw["make error"]];
]
MakeCC[graph_Graph,kSys_Integer,extraStabilizedVertices_List:{}]:=With[{
instanceFileContent=ExportSymmetricCCInstance[graph,kSys,extraStabilizedVertices]
},

Export[ParallelBuildDirectory[]<>CCCUSTOMINSTANCEFILE,instanceFileContent,"Text"];
MakeCC[kSys,VertexCount@graph-kSys]
]


(* ::Input::Initialization:: *)
Clear[RunCC]
RunCC[graph_Graph,kSys_Integer,extraStabilizedVertices:_List:{}]:=With[{},
If[!MakeCC[graph,kSys,extraStabilizedVertices],Return[False]];

With[{
ccResult=BuglessRunProcess["CoffeeCode",All,ProcessDirectory->ParallelBuildDirectory[]]
},
If[ccResult["ExitCode"]!=0\[Or]ccResult["StandardError"]!="",Echo[ccResult];Throw["run error"]];

ImportString[ccResult["StandardOutput"],"RawJSON"]
]
]


(* ::Input::Initialization:: *)
Clear[CCPath]
CCPath[kSys_Integer,kEnv_Integer]:=CCPath[kSys,kEnv]=CCRELEASEPATHFULLSOLVER<>"CoffeeCode."<>ToString[kSys]<>"."<>ToString[kEnv]

(* RunCC overload *)
RunCC[graph_Graph,kSys_Integer,"Full",extraStabilizedVertices_List:{}(*ignored*)]:=Module[{
adjacencyMatrix=ExportAdjacencyMatrix@graph,
kEnv=VertexCount@graph-kSys
},
With[{
executable=CCPath[kSys,kEnv]
},

Assert[FileExistsQ[executable]];

With[{
ccResult=BuglessRunProcess[executable,All,adjacencyMatrix]
},
If[ccResult["ExitCode"]!=0\[Or]ccResult["StandardError"]!="",Echo[ccResult];Throw["run error"]];

ImportString[ccResult["StandardOutput"],"RawJSON"]
]
]
]


(* ::Input::Initialization:: *)
Clear[CoeffArrayToPoly,MultArrayToPoly]
CoeffArrayToPoly[{q1_,__},{coeffs__Integer}]:=Sum[{coeffs}[[i]] q1^(i-1),{i,1,Length@{coeffs}}]
CoeffArrayToPoly[{q1_,q2_,q3_},{coeffs__List}]:=(#1 q1^#2[[1]] q2^#2[[2]] q3^#2[[3]])&@@@{coeffs}//Total
MultArrayToPoly[{q0_,qs__}]:=Function[{mult,list},{
mult,
q0 CoeffArrayToPoly[{qs},list]
}]


(* ::Input::Initialization:: *)
Clear[PauliActionCC]
PauliActionCC[graph_Graph,kSys_Integer,p_,symmetricSolver_:True,extraStabilizedVertices:_List:{}]:=Module[{
kTot=VertexCount@graph,
kEnv=VertexCount@graph-kSys,

pp=If[ListQ@p,p,{1-p,p/3,p/3,p/3}],
\[Lambda],\[Lambda]a,\[Lambda]m,\[Lambda]ma
},

If[kEnv>kSys\[Or]kSys<1\[Or]kTot<2,Throw["invalid params"]];

(* get lambda and lambda_a *)
{{\[Lambda]m,\[Lambda]},{\[Lambda]ma,\[Lambda]a}}=With[{
result=If[symmetricSolver,
RunCC[graph,kSys,extraStabilizedVertices],
RunCC[graph,kSys,"Full"]
],
toPolyFn=MultArrayToPoly[{pp[[1]]^kSys,pp[[2]]/pp[[1]],pp[[3]]/pp[[1]],pp[[4]]/pp[[1]]}]
},


{
toPolyFn@@@result["lambda"]//Transpose,
toPolyFn@@@result["lambda_a"]//Transpose
}
];

(* final expression in the variables given *)
Thread/@Thread[{\[Lambda],\[Lambda]a/2^(kTot-kSys)}->{\[Lambda]m,\[Lambda]ma}]
]


(* ::Input::Initialization:: *)
(* presets for important channels *)
CHANNELGENERIC={q0,q1,q2,q3};
CHANNELDEPOLARIZING={1-3p,p,p,p};
CHANNELBB84={1-2p-p^2,p-p^2,p^2,p-p^2};
CHANNEL2PAULI={1-2p,p,0,p};


(* ::Input::Initialization:: *)
Clear[ShannonEntropyMult,CIMult]
ShannonEntropyMult[l_List]:=l/.Rule[0,_]:>0/.Rule[term_,mult_]:>-mult*term Log2[term]//Total
CIMult[\[Lambda]_,\[Lambda]a_]:=(ShannonEntropyMult[\[Lambda]a]-ShannonEntropyMult[\[Lambda]])/Log2@Total@(Last/@List@@\[Lambda]a)


(* ::Input::Initialization:: *)
Clear[CIThreshold]
CIThreshold[pa_List,p_,compile_:False]:=Module[{
targetFunction
},
If[!compile,
targetFunction=Function[p,CIMult@@pa] (* this scope conflict is intended *)
];

FindRoot[targetFunction[p],
{p,.04,.16},
Method->"Secant",
AccuracyGoal->9,
MaxIterations->1000,
WorkingPrecision->$MachinePrecision
]
]

CIThreshold[g_Graph,kSys_Integer,p_,symmetricSolver_:True,extraStabilizedVertices_List:{}]:=With[{
pa=PauliActionCC[g,kSys,p,symmetricSolver,extraStabilizedVertices]
},
CIThreshold[pa,p]
]


(* ::Input::Initialization:: *)
Clear[ConcatGraph]
ConcatGraph[g1_Graph,g2_Graph,v1_Integer,v2_Integer]:=With[{
vcount1=VertexCount@g1,
vertices2=VertexList@g2,
edges2=EdgeList@g2
},
With[{
vmap=If[#==v2,v1,If[#<v2,#+vcount1,#+vcount1-1]]&
},
GraphUnion[g1,Graph[vmap/@vertices2,Map[vmap,edges2,{2}]]]
]
]

ConcatGraph[g1_Graph,g2_Graph,v1_List,v2_Integer]:=
Fold[f,{g1,g2,v2},v1]//.f[{gg1_,gg2_,vv2_},vv1_]:>{ConcatGraph[gg1,gg2,vv1,vv2],gg2,vv2}//First


(* ::Input::Initialization:: *)
Clear[All2LevelGraphs]
All2LevelGraphs[n_,minChildSize_:2]:=With[{
partitions=Cases[
IntegerPartitions[n-2,{2,n},Range[minChildSize,n]],
{x_,xs__}/;x>=Length[{xs}]
]
},
With[{
head=StarGraph[First@#+1,VertexLabels->"Index"],
legs=StarGraph[##+1,VertexLabels->"Index"]&/@Rest@#
},
<|
"legcount"->Length@legs,
"graph"->Fold[f,head,{Range@Length@legs+1,legs}\[Transpose]]//.f[g1_Graph,{i_Integer,g2_Graph}]:>ConcatGraph[g1,g2,{i},1]
|>
]&/@partitions
]


(* ::Input::Initialization:: *)
Clear[GraphHash]
GraphHash[g_Graph,kSys_Integer]:=StringTrim@ExportString[ExportString[CanonicalizeByOrbits[g,kSys],"Graph6"],"Base64"]//StringReplace[#,{"/"->"-"}]&


(* ::Input::Initialization:: *)
With[{
kSysMax=10
},

graphsToTest=ParallelTable[
With[{
graphs=All2LevelGraphs[kSys]
},
Thread[{VertexCount@#["graph"],#["legcount"],AllHairGraphs[#["graph"],Subsets[##,{1}]&]}]&/@graphs
],
{kSys,8,kSysMax}
]//Flatten[#,2]&//SortBy[First]
];

Print["Testing "<>ToString@Length@graphsToTest<>" graphs"]


(* ::Input::Initialization:: *)
count=0;
SetSharedVariable[count];
ParallelMap[Apply@Function[{kSys,legcount,graph},With[{
localcounter=count++//(Echo["starting "<>ToString@#];#)&,
pa=Catch@PauliActionCC[
graph,
kSys,
p,
kSys>14  (* symmetric solver *),
Range@legcount+1 (* extra stabilized vertices *)
],
hash=GraphHash[graph,kSys]
},

Echo["done PA for "<>ToString@localcounter];
If[StringQ[pa],Echo[pa];Return[]]; (* handle errors *)

With[{
pdep=p/.CIThreshold[pa,p]
},
 <|
"graph"->graph,
"kSys"->kSys,
"pauli action"->pa,
"threshold dep"->pdep
|>//Put[#,RESULTSPATH<>"graph-"<>hash<>"."<>ToString@kSys]&;

Echo["done "<>ToString@localcounter<>": "<>hash<>", p="<>ToString@pdep];
]

]],graphsToTest,Method->"FinestGrained"];

