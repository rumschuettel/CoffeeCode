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


(* ::Subsection:: *)
(*Paths and Unicorns*)


(* ::Text:: *)
(*These lines HAVE TO BE MODIFIED.*)
(**)
(*Edit scripts/make-build-dir.sh to activate a valid build environment for building with a modern gcc (>8.2). You can remove all the conda... blabla lines there if the compilation environment is installed globally anyhow. The options you pass to cmake there have to match those that you would pick when setting up a solver directory manually with ccmake. This also includes setting the paths for the boost and nauty libraries (via the options -D PATH_BOOST=/path/to/boost -D PATH_NAUTY=/path/to/nauty).*)


(* ::Input::Initialization:: *)
Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]


(* ::Input::Initialization:: *)
(* fix paths in this file *)
Get["CCInterfacePaths.m",Path->CurrentDir];


(* ::Input:: *)
(*(* execute this cell to open the build script for manual editing *)*)
(*SystemDialogInput["FileOpen",PARALLELBUILDSETUPSCRIPT]*)


(* ::Subsection:: *)
(*Build Environments*)


(* ::Input::Initialization:: *)
On[Assert];(* enable sanity checks instead of failing silently *)


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
StringForm[a_,b_Parallel`Kernels`kernel]:=""



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
(*(* for iterating over tuples, we'll only need the system vertices RestrictGroup, not the full group *)*)


(* ::Input::Initialization:: *)
Clear[RestrictGroup]
RestrictGroup[group_PermutationGroup,kSys_Integer]:=RestrictGroup[group,kSys]=group/.{p_Integer/;p>kSys->Nothing}/.Cycles[{}]->Nothing


(* ::Subsection:: *)
(*Plot Environment Vertices in Red*)


(* ::Input::Initialization:: *)
Clear[EnvironmentPlot]
EnvironmentPlot[g_Graph,systemSize_,vl_:None]:=Graph[g,VertexStyle->(
Thread[VertexList[g][[-(VertexCount@g-systemSize);;]]->Red]
),VertexSize->.2,VertexLabels->vl]


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
(*RestrictGroup[GroupForGraph[%,2],2]*)
(*OrbitsHaveMultipleHairsQ[%,%%,2]*)


(* ::Section:: *)
(*CoffeeCode Link [See individual sections for tests]*)


(* ::Input::Initialization:: *)
CCCUSTOMINSTANCEFILE="cc-instance-custom.h";


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
(*RestrictGroup[GroupForGraph[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {1, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {3 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2],2]*)
(*TrivialSGSTransversalMMForm[GroupOrbits@%,2]*)


(* ::Input::Initialization:: *)
(* we use these two heads to determine the type of variable passed to CoffeeCode *)
Unprotect[Numeric,Exact];Clear[Numeric,Exact];Protect[Numeric,Exact];

Clear[ExportSymmetricCCInstance]
(* we don't re-check pType for validity here, see sanity checks in PauliAction *)
ExportSymmetricCCInstance[graph_Graph,kSys_Integer,pType_,extraStabilizedVertices_List:{},name_:"graphstate_instance"]:=Module[{
group=RestrictGroup[GroupForGraph[graph,kSys,extraStabilizedVertices],kSys]//Echo[#,"group"]&,
kTot=VertexCount@graph,
poly
},

(* pType determines polynomial type used *)
poly=If[StringQ[pType],
If[pType=="Mono",
"using Polynomial = UnivariatePolynomial",
"using Polynomial = MultivariatePolynomial"
],

"constexpr static "<>
If[Length@First@pType==1,
"UnivariateSamples",
"MultivariateSamples"
]<>
"<"<>ToString[Length@pType]<>"> SamplePoints{"<>
ToString[NumberForm[List@@pType,10]]<>
"};\n"<>
"using Polynomial = SampledPolynomial<SamplePoints>"

];

(* pick trivial or nauty solver/iterator *)
If[!OrbitsHaveMultipleHairsQ[group,graph,kSys]\[And]ProductOfSymmetricGroupsQ@group,
(* write Orbits as expected for TrivialSGS Solver *)
With[{
orbits=GroupOrbits@group,
graphCanonicalized=CanonicalizeByOrbits[graph,kSys,extraStabilizedVertices]
},

If[Length@orbits==0\[Or]Max[Length/@orbits]==1,Echo["WARN: symmetry group is 0. Not a valid symmetric problem instance."];Throw["invalid instance"]
];
Echo["Product symmetry group detected. CanonicalImage=trivial."];

"struct "<>name<>" {\n"<>
TrivialSGSTransversalMMForm[orbits,kSys]<>"\n"<>
"constexpr static size_t k_sys = "<>ToString[kSys]<>", k_env = "<>ToString[kTot-kSys]<>";\n"<>
"constexpr static AdjacencyMatrixT<"<>ToString[kTot]<>"> adjacency_matrix{"<>ToString[Normal@AdjacencyMatrix@graphCanonicalized]<>"};\n"<>
poly<>";\n"<>
"};"
]
,
(* write SGSGenerators as expected for Nauty *)
With[{},
Echo["Nested symmetry group or multiple hairs per orbit detected. CanonicalImage=nauty."];

"struct "<>name<>" {\n"<>
"constexpr static size_t k_sys = "<>ToString[kSys]<>", k_env = "<>ToString[kTot-kSys]<>";\n"<>
"constexpr static AdjacencyMatrixT<"<>ToString[kTot]<>"> adjacency_matrix{"<>ToString[Normal@AdjacencyMatrix@graph]<>"};\n"<>
poly<>";\n"<>
"};"
]
]
]


(* ::Input:: *)
(*(* numeric sampled polynomial *)*)
(*g=Graph[ConcatGraph[StarGraph[10],StarGraph[5],{2,3},1],VertexLabels->"Index"]*)
(*ExportSymmetricCCInstance[g,VertexCount@g-1,Numeric[{.11111111, .2, .3},{.33, .44, .55},{.1, .2, .3}]]*)


(* ::Input:: *)
(*(* can enforce to use the non-nauty symmetric solver by flagging certain vertices as not symmetric *)ExportSymmetricCCInstance[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,"Mono",{2}]*)


(* ::Input:: *)
(*(* multiple hairs per orbit can mean connecting to the same environment vertex. *)*)
(*ExportSymmetricCCInstance[Graph[{1, 2, 3, 4, 5}, {Null, {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 1}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {5 -> RGBColor[1, 0, 0]}}],4,"Multi"]*)


(* ::Input:: *)
(*(* it can also mean two different hairs. *)*)
(*ExportSymmetricCCInstance[Graph[{2, 3, 1, 4}, {Null, {{3, 1}, {1, 2}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {1 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2,"Mono"]*)


(* ::Input:: *)
(*(* we can also numerically sample from this polynomial; easiest is just one value for p=p1=p2=p3 *)*)
(*ExportSymmetricCCInstance[Graph[{2, 3, 1, 4}, {Null, {{3, 1}, {1, 2}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {1 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2,Numeric[{.1905}]]*)


(* ::Input:: *)
(*(* or multiple samples, but still assuming p1=p2=p3 *)*)
(*ExportSymmetricCCInstance[Graph[{2, 3, 1, 4}, {Null, {{3, 1}, {1, 2}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {1 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2,Numeric[{.1904}, {.001}, {.001}]]*)


(* ::Input:: *)
(*(* or one sample, but for different ps *)*)
(*ExportSymmetricCCInstance[Graph[{2, 3, 1, 4}, {Null, {{3, 1}, {1, 2}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {1 -> RGBColor[1, 0, 0], 4 -> RGBColor[1, 0, 0]}}],2,Numeric[{.1904,.1905,.1906}]]*)


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
MakeCC[graph_Graph,kSys_Integer,pType_,extraStabilizedVertices_List:{}]:=With[{
instanceFileContent=ExportSymmetricCCInstance[graph,kSys,pType,extraStabilizedVertices]
},

Export[ParallelBuildDirectory[]<>CCCUSTOMINSTANCEFILE,instanceFileContent,"Text"];
MakeCC[kSys,VertexCount@graph-kSys]
]


(* ::Input:: *)
(*(* this fails if the folder for the symmetric solver is not already existent with a valid cc-instance-custom.h: MakeCC[9,1] *)*)


(* ::Input:: *)
(*(* try "Mono", "Multi" *)*)
(*MakeCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,"Mono"]*)


(* ::Input:: *)
(*(* try various numerics as above *)*)
(*MakeCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,Numeric[{.1,.11, .12}]]*)


(* ::Input::Initialization:: *)
Clear[RunCC]
RunCC[graph_Graph,kSys_Integer,pType_,extraStabilizedVertices:_List:{}]:=With[{},
If[!MakeCC[graph,kSys,pType,extraStabilizedVertices],Return[False]];

With[{
ccResult=BuglessRunProcess["CoffeeCode",All,ProcessDirectory->ParallelBuildDirectory[]]
},
If[ccResult["ExitCode"]!=0\[Or]ccResult["StandardError"]!="",Echo[ccResult];Throw["run error"]];

Check[
ImportString[ccResult["StandardOutput"],"RawJSON"],
Echo[ccResult];
Throw["import error"]
]
]
]


(* ::Input:: *)
(*RunCC[Graph[{1, 2, 3, 4, 5, 6, 7, 8}, {Null, {{1, 2}, {2, 3}, {2, 4}, {1, 5}, {5, 6}, {5, 7}, {1, 8}}}, {FormatType -> TraditionalForm, ImageSize -> {100, 100}, AlignmentPoint -> Center, GraphLayout -> {"Dimension" -> 2}, ImageSize -> {100, 100}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {8 -> RGBColor[1, 0, 0]}}],7,Numeric[{.19}],{2}]*)


(* ::Subsection:: *)
(*Full Solver Interface*)


(* ::Input::Initialization:: *)
Clear[CCPath]
CCPath[kSys_Integer,kEnv_Integer,pType_String]:=CCPath[kSys,kEnv,pType]=If[pType=="Mono",CCRELEASEPATHFULLSOLVER,CCRELEASEPATHFULLSOLVERARBCH]<>"CoffeeCode."<>ToString[kSys]<>"."<>ToString[kEnv]

(* RunCC overload *)
RunCC[graph_Graph,kSys_Integer,pType_String,"Full"]:=Module[{
adjacencyMatrix=ExportAdjacencyMatrix@graph,
kEnv=VertexCount@graph-kSys
},
With[{
executable=CCPath[kSys,kEnv,pType]
},

If[!FileExistsQ[executable],
Echo[executable<>" not found"];
Throw["run error"]
];

With[{
ccResult=BuglessRunProcess[executable,All,adjacencyMatrix]
},
If[ccResult["ExitCode"]!=0\[Or]ccResult["StandardError"]!="",Echo[ccResult];Throw["run error"]];

Check[
ImportString[ccResult["StandardOutput"],"RawJSON"],
Echo[ccResult];
Throw["import error"]
]
]
]
]


(* ::Input:: *)
(*RunCC[Graph[{1, 2, 3, 4, 5, 6, 7}, {Null, {{1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 6}, {4, 7}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {7 -> RGBColor[1, 0, 0]}}],6,"Mono","Full"]*)


(* ::Input:: *)
(*RunCC[Graph[{1, 2, 3, 4, 5, 6, 7}, {Null, {{1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 6}, {4, 7}}}, {FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {None}, VertexSize -> {0.2}, VertexStyle -> {7 -> RGBColor[1, 0, 0]}}],6,"Multi","Full"]*)


(* ::Subsection:: *)
(*Importing JSON format*)


(* ::Input::Initialization:: *)
Clear[CoeffArrayToPoly,MultArrayToPoly]
CoeffArrayToPoly[{q1_,__},{coeffs__Integer}]:=Sum[{coeffs}[[i]] q1^(i-1),{i,1,Length@{coeffs}}]
CoeffArrayToPoly[{q1_,q2_,q3_},{coeffs__List}]:=(#1 q1^#2[[1]] q2^#2[[2]] q3^#2[[3]])&@@@{coeffs}//Total
MultArrayToPoly[{q0_,qs__}]:=Function[{mult,list},{
mult,
q0 CoeffArrayToPoly[{qs},list]//Simplify
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

pp=If[MatchQ[p,_Exact],
(* either of the form
Exact[p]
Exact[p/3,q/5,1-q,pq]
*)
If[!(Length@p==1\[Or]Length@p==4),
Echo["either one or four symbolic parameters for Exact[...]"];
Throw["wrong Exact format"]
];

(* create list of length four *)
If[Length@p==1,
{1-p[[1]],p[[1]]/3,p[[1]]/3,p[[1]]/3},
List@@p
],

(* or of the form of triples/singles
Numeric[{{.3, .1, .1}, {.2, .2, .15}, ...}]
*)
If[!MatchQ[p,_Numeric],
Echo["either Exact or Numeric p expected"];
Throw["wrong format for p"]
];
If[!symmetricSolver,
Echo["numeric argument for full solver not supported"];
Throw["wrong Numeric argument"]
];

If[Length@p==0,Throw["no sample points given"]];
If[!Equal@@Length/@p,Throw["not all sample points have the same parameter count"]];
If[Length@First@p!=3,Throw["sample points have to be 4-tuples"]];

(* create list of length 4 *)
List@@p/.{{p1_Real,p2_Real,p3_Real}:>{1-p1-p2-p3,p1,p2,p3}}
],
pType,

\[Lambda],\[Lambda]a,\[Lambda]m,\[Lambda]ma
},

If[kEnv>kSys\[Or]kSys<1\[Or]kTot<2,Throw["invalid params"]];
pType=If[MatchQ[p,_Numeric],
Numeric@@(({#[[2]]/#[[1]],#[[3]]/#[[1]],#[[4]]/#[[1]]})&/@pp),
If[Length@p==1,"Mono","Multi"]
];

(* get lambda and lambda_a *)
With[{
result=If[symmetricSolver,
RunCC[graph,kSys,pType,extraStabilizedVertices],
RunCC[graph,kSys,pType,"Full"]
]
},

If[MatchQ[p,_Numeric],
(* return precalculated values *)
Precalculated[
((First/@pp)^kSys Last@result["lambda"])->First@result["lambda"],
((First/@pp)^kSys Last@result["lambda_a"])-> First@result["lambda_a"]
]
,
(* calculate exact polynomial *)
With[{
toPolyFn=MultArrayToPoly[{pp[[1]]^kSys,pp[[2]]/pp[[1]],pp[[3]]/pp[[1]],pp[[4]]/pp[[1]]}]
},

{{\[Lambda]m,\[Lambda]},{\[Lambda]ma,\[Lambda]a}}={
toPolyFn@@@result["lambda"]//Transpose,
toPolyFn@@@result["lambda_a"]//Transpose
};

(* final expression in the variables given *)
Thread/@Thread[{\[Lambda],\[Lambda]a/2^(kTot-kSys)}->{\[Lambda]m,\[Lambda]ma}]
]
]
]
]


(* ::Input:: *)
(*PauliActionCC[Graph[{1, 2, 3, 4}, {Null, {{1, 2}, {2, 3}, {2, 4}}}, {FormatType -> TraditionalForm, FormatType -> TraditionalForm, GraphLayout -> {"Dimension" -> 2}, VertexLabels -> {"Index"}, VertexSize -> {0.2}, VertexStyle -> {4 -> RGBColor[1, 0, 0]}}],3,Numeric@@({Range@50/100}//Transpose//N)]*)


(* ::Chapter:: *)
(*User Interface*)


(* ::Section:: *)
(*Entropy and CI*)


(* ::Input::Initialization:: *)
(* presets for important channels *)
CHANNELGENERIC={q0,q1,q2,q3};
CHANNELDEPOLARIZING={1-p,p/3,p/3,p/3};
CHANNELBB84={1-2p+p^2,p-p^2,p^2,p-p^2};
CHANNEL2PAULI={1-p,p/2,0,p/2};


(* ::Input::Initialization:: *)
Clear[ShannonEntropyMult,CIMult]
ShannonEntropyMult[l_List]:=l/.Rule[0,_]:>0/.Rule[term_,mult_]:>-mult*term Log2[term]//Total
CIMult[{\[Lambda]_,\[Lambda]a_}]:=(ShannonEntropyMult[\[Lambda]a]-ShannonEntropyMult[\[Lambda]])/Log2@Total@(Last/@List@@\[Lambda]a)
CIMult[Precalculated[\[Lambda]_->_,\[Lambda]a_->multa_]]:=(\[Lambda]a-\[Lambda])/Log2@multa


(* ::Input::Initialization:: *)
Clear[ZeroCrossings]
ZeroCrossings[l_List]:=Module[{t,u,v},
t={Sign[l],Range[Length[l]]}//Transpose;
u=Select[t,First[#]!=0&];
v=SplitBy[u,First]; 
{Most[Max[#[[All,2]]]&/@v],Rest[Min[#[[All,2]]]&/@v]}//Transpose
]


(* ::Input::Initialization:: *)
Clear[CIThreshold]
CIThreshold[pa_List,p_,compile_:False]:=Module[{
targetFunction
},
If[!compile,
targetFunction=Function[p,CIMult@pa] (* this scope conflict is intended *)
];

FindRoot[targetFunction[p],
{p,.04,.16},
Method->"Secant",
AccuracyGoal->9,
MaxIterations->1000,
WorkingPrecision->MachinePrecision
]
]

(* this cannot be exact, of course; maybe there's a better way of doing this, for now we just return the data points surrounding wherever the CI changes sign *)
CIThreshold[pa_Precalculated,p_Numeric,range_:All]:=Module[{
data=(CIMult@pa)[[range]],
pts=(List@@p)[[range]],
crossings
},

crossings=ZeroCrossings[data];

Between[{First@pts[[#1]],First@pts[[#2]]}]&@@@crossings
]

CIThreshold[g_Graph,kSys_Integer,p_,symmetricSolver_:True,extraStabilizedVertices_List:{}]:=With[{
pa=PauliActionCC[g,kSys,p,symmetricSolver,extraStabilizedVertices]
},
CIThreshold[pa,p]
]


(* ::Input:: *)
(*ConcatGraph[StarGraph[15],StarGraph[15],{2},1]*)
(*(* for numeric input, the threshold is just a list of points where the function value crosses *)*)
(*CIThreshold[%,VertexCount@%-1,Numeric@@({Range@50/100}//N//Transpose)]*)


(* ::Section::Closed:: *)
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
(*Unique Graph Hash*)


(* ::Text:: *)
(*This graph hash can be used in filenames; it canonicalizes the graph first to ensure that for different isomorphic graphs one gets the same hash.*)
(*Note that, in principle, the filename can be decoded to the underlying graph by replacing - with /, Base64-decode, and then ImportString with Graph6 as format.*)


(* ::Input::Initialization:: *)
Clear[GraphHash]
GraphHash[g_Graph,kSys_Integer]:=StringTrim@ExportString[ExportString[CanonicalizeByOrbits[g,kSys],"Graph6"],"Base64"]//StringReplace[#,{"/"->"-"}]&


(* ::Subsection::Closed:: *)
(*Rooted Graph Product*)


(* ::Input::Initialization:: *)
Clear[RootedGraphProduct]
RootedGraphProduct[g1_Graph,g2_Graph,v1_Integer,v2_Integer]:=With[{
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

RootedGraphProduct[g1_Graph,g2_Graph,v1_List,v2_Integer]:=
Fold[f,{g1,g2,v2},v1]//.f[{gg1_,gg2_,vv2_},vv1_]:>{RootedGraphProduct[gg1,gg2,vv1,vv2],gg2,vv2}//First


(* ::Subsection:: *)
(*Concatenated Graphs*)


(* ::Text:: *)
(*Maybe you want to fill these in?*)


(* ::Section:: *)
(*Two-Level Tree Graphs*)


(* ::Subsection:: *)
(*Graph Generation*)


(* ::Input::Initialization:: *)
Clear[All2LevelGraphs]
(* all 2 level graphs with vertex count < n *)
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
"graph"->Fold[f,head,{Range@Length@legs+1,legs}\[Transpose]]//.f[g1_Graph,{i_Integer,g2_Graph}]:>RootedGraphProduct[g1,g2,{i},1]
|>
]&/@partitions
]


(* ::Subsection:: *)
(*Numeric Search Run*)


(* ::Input:: *)
(*(* get precomputed sample ranges *)*)
(*Get["p_values.m",Path->CurrentDir]*)


(* ::Input:: *)
(*(* the parameter is the resoultion *)*)
(*SAMPLES["depolarizing"][10]*)
(*SAMPLES["depolarizing"][5]*)


(* ::Input:: *)
(*Clear[MultiResolutionSamples]*)
(*MultiResolutionSamples[chan_,threshold_,n_:10]:=MultiResolutionSamples[chan]=With[{*)
(*coarse=SAMPLES[chan][10],*)
(*fine=SAMPLES[chan][100],*)
(*ultra=SAMPLES[chan][1000],*)
(*ludicrous=SAMPLES[chan][10000]*)
(*},*)
(*Sort@Union[*)
(*coarse,*)
(*Nearest[fine,threshold,n],*)
(*Nearest[ultra,threshold,n],*)
(*Nearest[ludicrous,threshold,n]*)
(*]*)
(*]*)


(* ::Input:: *)
(*(* we take finer and finer samples around the points of interest where we expect a crossing *)*)
(*MultiResolutionSamples["depolarizing",Rest@CHANNELDEPOLARIZING/.p->0.1904]*)
(*ListLogPlot[Norm[#-Rest@CHANNELDEPOLARIZING/.p->0.1904]&/@%]*)


(* ::Input:: *)
(*(* multi-resolution samples; we need to find better thresholds for BB84 and 2Pauli *)*)
(*samples={*)
(*MultiResolutionSamples["depolarizing",Rest@CHANNELDEPOLARIZING/.p->0.1904],*)
(*MultiResolutionSamples["BB84",Rest@CHANNELBB84/.p->0.1121],*)
(*MultiResolutionSamples["2Pauli",Rest@CHANNEL2PAULI/.p->0.1135]*)
(*};*)
(*sampleCounts=Length/@samples*)
(*combinedSamples=Join@@samples//N;*)
(**)


(* ::Input:: *)
(*RESULTSPATHLOCAL=RESULTSPATH<>"2level/";*)
(*FILEPREFIX="2level-numeric-";*)
(*Print["Saving output to "<>RESULTSPATHLOCAL<>FILEPREFIX<>"..."]*)


(* ::Input:: *)
(**)
(*Range[10,10]//Scan[Function[n,With[{*)
(*graphs=All2LevelGraphs[n]*)
(*},*)
(*Echo["all 2 level graphs with less than "<>ToString[n]<>" vertices: "<>ToString@Length[graphs]];*)
(*graphs//Scan[Function[g,With[{*)
(*G=HairGraph[g["graph"],{1}],*)
(*kSys=VertexCount@g["graph"]*)
(*},*)
(*Module[{*)
(*filename=RESULTSPATHLOCAL<>FILEPREFIX<>GraphHash[G,kSys]<>"."<>ToString[kSys],*)
(*prettyG=EnvironmentPlot[G,kSys],*)
(*pNumeric=Numeric@@combinedSamples ,*)
(*pa*)
(*},*)
(*If[FileExistsQ[filename],Echo["skipping "<>filename];Return[0,Module]];*)
(*Echo["running "<>filename];*)
(**)
(*pa=Catch@PauliActionCC[G,kSys,pNumeric(* numeric poly solver only works with the symmetric solver *)];*)
(*Echo["done PA for "<>filename];*)
(*If[StringQ[pa],Echo[pa];Return[]]; (* display errors *)*)
(**)
(*With[{*)
(*pdep=CIThreshold[pa,pNumeric,;;Length@samples[[1]]],*)
(*pbb84=CIThreshold[pa,pNumeric,Length@samples[[1]]+1;;Length@samples[[1]]+Length@samples[[2]]],*)
(*p2pauli=CIThreshold[pa,pNumeric,-Length@samples[[3]];;]*)
(*},*)
(* <|*)
(*"graph"->G,*)
(*"kSys"->kSys,*)
(*"pauli action"->Compress[pa],*)
(*"threshold dep"->pdep,*)
(*"threshold bb84"->pbb84,*)
(*"threshold 2pauli"->p2pauli*)
(*|>//Put[#,filename]&;*)
(**)
(*Echo["done "<>filename<>", pdep="<>ToString@pdep<>", pbb84="<>ToString@pbb84<>", p2pauli="<>ToString@p2pauli];*)
(*]]]]];*)
(*]]];*)
