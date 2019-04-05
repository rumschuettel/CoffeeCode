#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Clear[CurrentDir]
CurrentDir:=ExpandFileName[First[$ScriptCommandLine]];
Print["Working directory: ",CurrentDir];*)


(* ::Input::Initialization:: *)
(* fix paths in this file *)
(*Get["CCInterfacePaths.m",Path->CurrentDir];*)
Get["CCInterfacePaths.m"];


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


(* ::Input::Initialization:: *)
Clear[BuglessRunProcess]
(*BuglessRunProcess[exec_String,what_,input_String:"",opts:OptionsPattern[]]:=Module[{
process=StartProcess[exec,opts],
out
},
If[input\[NotEqual]"",WriteLine[process,input]];
While[ProcessStatus[process,"Running"],Pause[1]];
out=ReadString[process];
<|
"ExitCode"\[Rule]ProcessInformation[process]["ExitCode"],
"StandardError"\[Rule]"",
"StandardOutput"\[Rule]out
|>
]*)

BuglessRunProcess=RunProcess;

Unprotect[StringForm];
StringForm[a_,b_Parallel`Kernels`kernel]:=""


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
Clear[RestrictGroup]
RestrictGroup[group_PermutationGroup,kSys_Integer]:=RestrictGroup[group,kSys]=group/.{p_Integer/;p>kSys->Nothing}/.Cycles[{}]->Nothing


(* ::Input::Initialization:: *)
Clear[MaxGroupActionBase]
MaxGroupActionBase[group_PermutationGroup]:=MaxGroupActionBase[group]=List@@@GroupGenerators@group//Flatten//Max


(* ::Input::Initialization:: *)
Clear[SGSTransversal]
SGSTransversal[group_PermutationGroup,addEmptyGroupForMaxGAB_:True]:=Module[{
chain=GroupStabilizerChain[group],
sgs
},
(* the head element of the group stabilizer chain already delivers a transversal of a strong generating set *)
sgs=f@@@Partition[chain,2,1]/.f[Rule[bA_List,gA_PermutationGroup],Rule[bB_List,gB_PermutationGroup]]:>Complement[bB,bA]->PermutationGroup@If[GroupOrder@gB==1,
GroupGenerators@gA,
DeleteDuplicates@Flatten[RightCosetRepresentative[gB,#]&/@GroupElements[gA]]
]/.Cycles[{}]->Nothing;

(* we want to ensure that the transversal contains one element for the maximum index that the group acts on nontrivially, even if that element is the empty group *)
If[addEmptyGroupForMaxGAB\[And]Max@Flatten[First/@sgs]<MaxGroupActionBase[group],
sgs~Join~{{MaxGroupActionBase[group]}->PermutationGroup[{}]},
sgs
]
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
CCCUSTOMINSTANCEFILE="cc-instance-custom.h";


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
permutationIndices=If[Length@cycles>0,
PermutationList[#,kSys]&/@cycles,
{Range@kSys}
]
},

(* the pivot point is the first index for which the SGSGenerator acts nontrivially on (so one point past the points that are stabilized).
If the pivot point for the SGSGenerator is 2, that means we want to check indices 0 and 1 in C++; in MM this would correspond to indices {1, 2}, which is what \[LeftDoubleBracket];;2\[RightDoubleBracket] truncates to.
Therefore we validly don't subtract 1 from the pivot point, but do subtract 1 from all other indices *)
StringRiffle[Flatten@{
"\tSGSGenerator<"<>ToString@(#1)<>", Group<",
StringRiffle[("\t\tPermutation<"<>StringRiffle[ToString/@(##-1),","]<>">")&/@#2,",\n"],
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
group=RestrictGroup[GroupForGraph[graph,kSys,extraStabilizedVertices],kSys],
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
Echo["Product symmetry group detected. CanonicalImage=trivial."];

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
Echo["Nested symmetry group or multiple hairs per orbit detected. CanonicalImage=nauty."];

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
q0 CoeffArrayToPoly[{qs},list]//Simplify
}]


(* ::Input::Initialization:: *)
Clear[PauliActionCC]
PauliActionCC[graph_Graph,kSys_Integer,p_,symmetricSolver_:True,extraStabilizedVertices:_List:{}]:=Module[{
kTot=VertexCount@graph,
kEnv=VertexCount@graph-kSys,

pp=If[ListQ@p,p,{1-p,p/3,p/3,p/3}],
\[Lambda],\[Lambda]a,\[Lambda]m,\[Lambda]ma,

taylorF
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
Clear[GraphHash]
GraphHash[g_Graph,kSys_Integer]:=StringTrim@ExportString[ExportString[CanonicalizeByOrbits[g,kSys],"Graph6"],"Base64"]//StringReplace[#,{"/"->"-"}]&


(* ::Input:: *)
(*(* simulated annealing code *)*)


(* ::Input::Initialization:: *)
Clear[switchLeaf]
switchLeaf[G_Graph]:=Module[{
leavesIndices,
leafIndex,
leaf,

neighbour,

newG
},
(* get random vertex with degree 1. Note that vertex index is not necessarily equal to vertex name, so once we pick a random position of a vertex we still have to get the actual vertex name to delete *)
leavesIndices=Flatten@Position[VertexDegree@G,1];
leafIndex=RandomChoice@leavesIndices;
leaf=VertexList[G][[leafIndex]];

(* delete edge between leaf and neighbour *)
neighbour=First@AdjacencyList[G,leaf];
newG=EdgeDelete[G,leaf\[UndirectedEdge]neighbour];

(* add new random edge *)
newG=EdgeAdd[newG,leaf\[UndirectedEdge]RandomChoice@Complement[VertexList@G,{leaf,neighbour}]];

(* return with index labels. Vertex names should be avoided whenever possible *)
Graph[newG,VertexLabels->"Index"]
]


(* ::Input::Initialization:: *)
Clear[simAnn]
simAnn[costFunctionPolynomial_,kSys_Integer,depolProb_,initialG_Graph,maximumTemperature_Real:100,annealingSteps_Integer:100,annealingIterations_Integer:1000]:=
Module[{
curf,curG,
topf,topG
},
curG=initialG;
curf=costFunctionPolynomial@curG/.{p->depolProb};
topf=curf;
topG=curG;

Range[annealingSteps]//Map[Function[currentStep,Module[{
(* keep variables as local as possible *)
temperature=maximumTemperature/currentStep,
(*count the number of moves (accepting better function value) and hops (accepting worse function value) for diagnostics purposes*)
hopCounter=0,moveCounter=0
},

Print["Annealing step: ",currentStep,", Temp: ",N[temperature]];

Range[annealingIterations]//Map[Function[iteration,Module[{
newG=switchLeaf[curG],
newPoly,newf,hash
},
(* disallow disconnected graphs? *)
If[!ConnectedGraphQ@newG,
(*PrintTemporary["skipping disconnected"];*)
Return[Module]
];

newPoly=costFunctionPolynomial[newG];
newf=newPoly/.{p->depolProb};
hash=GraphHash[newG,kSys];
filename=RESULTSPATH<>"graph-"<>hash<>"."<>ToString@kSys;
If[!FileExistsQ[filename],
<|
"graph" -> newG,
"kSys" -> kSys,
"pauli action" -> Compress[newPoly]
|>//Put[#,filename]&;
];

(* new best found *)
If[newf>curf,
curf=newf;
curG=newG;
moveCounter++;
,

(* actually hop there only with some probability *)
If[Exp[-(curf-newf)/temperature]>RandomReal[],
curf=newf; 
curG=newG;
hopCounter++;
];
];

(* display that we found a new top choice this iteration *)
If[curf>topf,
topf=curf;topG=curG;Print["Iteration: ",iteration,", cI = ",topf];(*Print[Graph[topG,VertexLabels->"Name",VertexStyle->{VertexCount[topG]->Red}]];*)
];


]]];
Print[moveCounter," moves, ",hopCounter," hops"];
Print["Best cI: ",topf];
Print["\n"];
]]];

Graph[topG,VertexStyle->{VertexCount[topG]->Red}]
]


(* ::Input::Initialization:: *)
Clear[estimateInitialTemperature]
estimateInitialTemperature[costFunctionPolynomial_,depolProb_,initialG_Graph,numberTrials_Integer:100,initialAcceptanceProbability_Real:0.8]:=Module[{costFunction,functionEvaluations,initf},
costFunction[G_]:=costFunctionPolynomial[G]/.{p->depolProb};
initf=costFunction[initialG];
functionEvaluations=Table[costFunction[switchLeaf[initialG]]-initf,numberTrials];
-Mean[Select[functionEvaluations,0<#&]]/Log[initialAcceptanceProbability]
];


(* ::Input::Initialization:: *)
Clear[ShannonEntropyMult,CIMult]
ShannonEntropyMult[l_List]:=l/.Rule[0,_]:>0/.Rule[term_,mult_]:>-mult*term Log2[term]//Total
CIMult[\[Lambda]_,\[Lambda]a_]:=(ShannonEntropyMult[\[Lambda]a]-ShannonEntropyMult[\[Lambda]])/Log2@Total@(Last/@List@@\[Lambda]a)


(* ::Input::Initialization:: *)
(* collect parameters *)
args=Rest[$ScriptCommandLine];
pp=ToExpression[args[[1]]];
n=ToExpression[args[[2]]];
kSys=ToExpression[args[[3]]];


(* ::Input::Initialization:: *)
(* define objective function *)
cI[G_]:=cI[G]=CIMult@@PauliActionCC[G,kSys,p,False];


(* ::Input::Initialization:: *)
(* actual computation *)
H=Nest[switchLeaf,StarGraph[n],25];
T=estimateInitialTemperature[cI,pp,H,250];
Print[T];
G=simAnn[cI,kSys,pp,H,T,100,1000]
