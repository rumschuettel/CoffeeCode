#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*CoffeeCode Mathematica*)


(* ::Text:: *)
(*CoffeeCode Mathematica Interface, \[Copyright] 2019 Johannes Bausch*)
(*For joint Graph State project with Felix Leditzky*)


(* ::Text:: *)
(*This file can be exported as a module, and then imported in a script.*)


(* ::Title:: *)
(*CoffeeCode Interface*)


(* ::Section:: *)
(*Setup [Build paths to be fixed here; contains more instructions]*)


(* ::Text:: *)
(*Set the CoffeeCode release path below (where you execute make), and the file you configured as cc-instance-custom.h during ccmake. You can use the file open dialogue helper to obtain said paths.*)
(**)
(*Note that depending on whether the build directory is set up with SYMMETRIC_SOLVER or not will NOT determine what variant of CoffeeCode is used to execute the instance; you have to do this manually. Also note that the output will be wrong for general channels if OPTIMIZE_FOR_DEPOLARIZING is set as build flag.*)
(**)
(*Finally, you'll have to edit scripts/make-build-dir.sh to activate a valid conda environment for building with a modern gcc. You can remove all the source... blabla lines there if the compilation environment is installed globally anyhow.*)


(* ::Input:: *)
(*SystemDialogInput["Directory"]*)
(*SystemDialogInput["FileOpen"]*)


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


(* ::Input:: *)
(*ParallelBuildDirectory[]*)


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



(* ::Section:: *)
(*Basic Graph Functionality [install IGraphM once below, then forget]*)


(* ::Subsection:: *)
(*IGraphM*)


(* ::Text::Initialization:: *)
(*Uncomment the last line in the following cell, then execute to install IGraphM. This only has to be done once per machine.*)


(* ::Input:: *)
(*updateIGraphM[]:=Module[{json,download,target,msg},Check[json=Import["https://api.github.com/repos/szhorvat/IGraphM/releases","JSON"];*)
(*download=Lookup[First@Lookup[First[json],"assets"],"browser_download_url"];*)
(*msg="Downloading IGraph/M "<>Lookup[First[json],"tag_name"]<>" ...";*)
(*target=FileNameJoin[{CreateDirectory[],"IGraphM.paclet"}];*)
(*If[$Notebooks,PrintTemporary@Labeled[ProgressIndicator[Appearance->"Necklace"],msg,Right],Print[msg]];*)
(*URLSave[download,target],Return[$Failed]];*)
(*If[FileExistsQ[target],PacletManager`PacletInstall[target],$Failed]]*)
(*(*updateIGraphM[]*) (* uncomment this once and run to install IGraphM locally, then comment again *)*)
(**)


(* ::Input::Initialization:: *)
Needs["IGraphM`"]
ParallelNeeds["IGraphM`"];


(* ::Subsection::Closed:: *)
(*Group for Graph State*)


(* ::Text:: *)
(*GroupForGraph gives the full automorphism group of the underlying graph; we color the environment vertices in a different color to make sure the permutation group will be a product across this partitioning.*)


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


(* ::Input:: *)
(*(* for iterating over tuples, we'll only need the system vertices subgroup, not the full group *)*)


(* ::Input::Initialization:: *)
Clear[Subgroup]
Subgroup[group_PermutationGroup,kSys_Integer]:=Subgroup[group,kSys]=group/.{p_Integer/;p>kSys->Nothing}/.Cycles[{}]->Nothing


(* ::Subsection::Closed:: *)
(*Plot Environment Vertices in Red*)


(* ::Input::Initialization:: *)
Clear[EnvironmentPlot]
EnvironmentPlot[g_Graph,systemSize_]:=Graph[g,VertexStyle->(
Thread[(Range[VertexCount@g-systemSize]+systemSize)->Red]
),VertexSize->.2,VertexLabels->None]


(* ::Subsection::Closed:: *)
(*Create all Possible Choices of Hairs to the Environment*)


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


(* ::Input:: *)
(*AllHairGraphs[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}, {Null, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}, {1, 11}, {1, 12}, {1, 13}, {1, 14}, {1, 15}, {1, 16}, {2, 17}, {2, 18}, {2, 19}}}, {FormatType -> TraditionalForm, ImageSize -> {296., Automatic}}],Subsets[#,{1}]&]*)


(* ::Subsection::Closed:: *)
(*Strong Generating Set Transversal*)


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


(* ::Input:: *)
(*SGSTransversal[GroupForGraph[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, {Null, {{1, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}, {2, 7}, {3, 8}, {3, 9}, {3, 10}, {4, 11}, {4, 12}, {4, 13}, {1, 14}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {14 -> RGBColor[1, 0, 0]}}],13]]*)


(* ::Subsection::Closed:: *)
(*Detecting Products of Symmetric Groups*)


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



(* ::Input:: *)
(*GroupForGraph[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}, {Null, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}, {1, 11}, {1, 12}, {12, 13}, {12, 14}, {12, 15}, {12, 16}, {12, 17}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {17 -> RGBColor[1, 0, 0]}}],16]*)
(*ProductOfSymmetricGroupsQ@%*)
(*GroupOrbits@%%*)
(**)


(* ::Input:: *)
(*GroupForGraph[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, {Null, {{1, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}, {2, 7}, {3, 8}, {3, 9}, {3, 10}, {4, 11}, {4, 12}, {4, 13}, {1, 14}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {14 -> RGBColor[1, 0, 0]}}],13]*)
(*ProductOfSymmetricGroupsQ@%*)


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


(* ::Input:: *)
(*CanonicalizeByOrbits[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}, {Null, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}, {1, 11}, {1, 12}, {12, 13}, {12, 14}, {12, 15}, {12, 16}, {12, 17}}}, {FormatType -> TraditionalForm, ImageSize -> {373., Automatic}, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {17 -> RGBColor[1, 0, 0]}}],16]*)


(* ::Subsection::Closed:: *)
(*Detecting Multiple Hairs per Orbit*)


(* ::Input::Initialization:: *)
Clear[OrbitsHaveMultipleHairsQ]
OrbitsHaveMultipleHairsQ[group_PermutationGroup,g_Graph,kSys_Integer]:=Module[{
environmentVertices=VertexList[g][[kSys+1;;]],
hairNeighbourIndices
},
hairNeighbourIndices=VertexIndex[g,#]&/@VertexList@Subgraph[NeighborhoodGraph[g,environmentVertices],Complement[VertexList@g,environmentVertices]];

Intersection[#,hairNeighbourIndices]&/@GroupOrbits[group]//AnyTrue[#,Length@##>1&]&
]


(* ::Input:: *)
(*Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {1, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {3 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}]*)
(*Subgroup[GroupForGraph[%,2],2]*)
(*OrbitsHaveMultipleHairsQ[%,%%,2]*)


(* ::Section:: *)
(*CoffeeCode Link [See individual sections for tests]*)


(* ::Subsection:: *)
(*Exporting MM symbols*)


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


(* ::Input:: *)
(*Subgroup[GroupForGraph[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {1, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {3 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2],2]*)
(*TrivialSGSTransversalMMForm[GroupOrbits@%,2]*)


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


(* ::Input:: *)
(*(* can enforce to use the non-nauty symmetric solver by flagging certain vertices as not symmetric *)ExportSymmetricCCInstance[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,{2}]*)


(* ::Subsection:: *)
(*Build Interface*)


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


(* ::Input:: *)
(*MakeCC[9,1]*)


(* ::Input:: *)
(*MakeCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7]*)


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


(* ::Input:: *)
(*RunCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,{2}]*)


(* ::Subsection:: *)
(*Full Solver Interface*)


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


(* ::Input:: *)
(*RunCC[Graph[{1, 2, 3, 4, 5, 6, 7}, {Null, {{1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 6}, {4, 7}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {7 -> RGBColor[1, 0, 0]}}],6,"Full"]*)


(* ::Subsection:: *)
(*Importing JSON format*)


(* ::Input::Initialization:: *)
Clear[CoeffArrayToPoly,MultArrayToPoly]
CoeffArrayToPoly[{q1_,__},{coeffs__Integer}]:=Sum[{coeffs}[[i]] q1^(i-1),{i,1,Length@{coeffs}}]
CoeffArrayToPoly[{q1_,q2_,q3_},{coeffs__List}]:=(#1 q1^#2[[1]] q2^#2[[2]] q3^#2[[3]])&@@@{coeffs}//Total
MultArrayToPoly[{q0_,qs__}]:=Function[{mult,list},{
mult,
q0 CoeffArrayToPoly[{qs},list]
}]


(* ::Input:: *)
(*MultArrayToPoly[{q0,q1,q2,q3}]@@@{*)
(*{8000,{{1,{1,1,1}},{2,{2,1,1}},{11,{0,0,0}},{3,{15,2,2}}}},*)
(*{1300,{{1,{1,1,1}},{2,{2,1,1}},{11,{0,0,0}},{3,{15,2,2}}}}*)
(*}*)


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


(* ::Input:: *)
(*PauliActionCC[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {2, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {4 -> RGBColor[1, 0, 0]}}],3,p]*)


(* ::Input:: *)
(*1//AbsoluteTiming*)


(* ::Chapter:: *)
(*User Interface*)


(* ::Section:: *)
(*Entropy and CI*)


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

FindRoot[targetFunction[x],
{x,.04,.16},
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


(* ::Input:: *)
(*ConcatGraph[StarGraph[15],StarGraph[15],{2},1]*)
(*CIThreshold[%,VertexCount@%-1,p]*)


(* ::Section:: *)
(*Interface Tutorial*)


(* ::Subsection:: *)
(*Manual Export*)


(* ::Subsubsection:: *)
(*Create graphs by hand*)


(* ::Text:: *)
(*To create graphs by hand, look up the following functions:*)


(* ::Input:: *)
(*AdjacencyGraph*)
(*SparseArray*)
(*PathGraph*)


(* ::Text:: *)
(*Or directly use MM's Graph object, where edges are entered with ESC ue ESC*)


(* ::Input:: *)
(*Graph[{1\[UndirectedEdge]2,3\[UndirectedEdge]4,3\[UndirectedEdge]5,1\[UndirectedEdge]3}]*)


(* ::Subsubsection:: *)
(*Create a random graph, add a hair, and plot.*)


(* ::Input:: *)
(*CompleteGraph[10];*)
(*HairGraph[%,{10,9}];*)
(*graph=EnvironmentPlot[%,10]*)


(* ::Subsubsection:: *)
(*Export instance*)


(* ::Text:: *)
(*Just adjacency matrix for Full Solver*)


(* ::Input:: *)
(*ExportAdjacencyMatrix[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {2, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {4 -> RGBColor[1, 0, 0]}}]]*)


(* ::Text:: *)
(*SGS for Symmetric Solver*)


(* ::Input:: *)
(*ExportSymmetricCCInstance[Graph[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}, {Null, {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}, {1, 11}, {1, 12}, {12, 13}, {12, 14}, {12, 15}, {12, 16}, {12, 17}}}, {FormatType -> TraditionalForm, ImageSize -> {373., Automatic}, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {17 -> RGBColor[1, 0, 0]}}],16]*)


(* ::Subsection:: *)
(*Automatic Export, Build and Import*)


(* ::Subsubsection:: *)
(*Get CI for a Graph*)


(* ::Text:: *)
(*Remove the semicolon to see full analytic expression. If you give "p", that corresponds to giving CHANNELDEPOLARIZING where all the terms are equal; alternatively you can specify an array of four individual Pauli coefficients, e.g. {1-3p^2,p^2,p^2,p^2}, or use a predefined block (see section "Entropy and CI" above).*)


(* ::Input:: *)
(*CIMult@@PauliActionCC[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {2, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {4 -> RGBColor[1, 0, 0]}}],3,p];*)


(* ::Text:: *)
(*For plotting, use e.g.*)


(* ::Input:: *)
(*CIMult@@PauliActionCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,p];*)
(*Plot[%,{p,0,1}]*)


(* ::Text:: *)
(*The same graph is small enough to be faster with the full solver:*)


(* ::Input:: *)
(*CIMult@@PauliActionCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,p,False];*)
(*Plot[%,{p,0,1}]*)


(* ::Text:: *)
(*If the symmetry group is 0, we cannot solve with the symmetric solvers; an error is thrown that you can catch.*)


(* ::Input:: *)
(*Catch@CIThreshold[Graph[{1, 2, 3, 4, 5, 6, 7}, {Null, {{1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 6}, {4, 7}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {7 -> RGBColor[1, 0, 0]}}],6,p]*)


(* ::Text:: *)
(*Solving fully still works.*)


(* ::Input:: *)
(*Catch@CIThreshold[Graph[{1, 2, 3, 4, 5, 6, 7}, {Null, {{1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 6}, {4, 7}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {7 -> RGBColor[1, 0, 0]}}],6,p,False]*)


(* ::Chapter:: *)
(*Optimization Runs*)


(* ::Text:: *)
(*Everything below can be safely deleted and is not necessary for the CoffeeCode interface. It contains some helpers for graph generation.*)


(* ::Subsection:: *)
(*Rooted Graph Product/Concatenated Codes*)


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


(* ::Section:: *)
(*Two-Level Tree Graphs*)


(* ::Subsection:: *)
(*Graph Generation*)


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

