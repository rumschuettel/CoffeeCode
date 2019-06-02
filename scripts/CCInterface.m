(* ::Package:: *)

(* ::Input::Initialization:: *)
Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]


(* ::Input::Initialization:: *)
(* fix paths in this file *)
Get["CCInterfacePaths.m",Path->CurrentDir];


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
IGBlissAutomorphismGroup[{graph,"VertexColors"->ConstantArray[0,kSys]~Join~ConstantArray[kTot^2,kTot-kSys]+extraColors}]
]


(* ::Input::Initialization:: *)
Clear[RestrictGroup]
RestrictGroup[group_PermutationGroup,kSys_Integer]:=RestrictGroup[group,kSys]=group/.{p_Integer/;p>kSys->Nothing}/.Cycles[{}]->Nothing


(* ::Input::Initialization:: *)
Clear[EnvironmentPlot]
EnvironmentPlot[g_Graph,systemSize_,vl_:None]:=Graph[g,VertexStyle->(
Thread[VertexList[g][[-(VertexCount@g-systemSize);;]]->Red]
),VertexSize->.2,VertexLabels->vl]


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


(* ::Input::Initialization:: *)
Clear[CoeffArrayToPoly,MultArrayToPoly]
CoeffArrayToPoly[{q1_,__},{coeffs__Integer}]:=Sum[{coeffs}[[i]] q1^(i-1),{i,1,Length@{coeffs}}]
CoeffArrayToPoly[{q1_,q2_,q3_},{coeffs__List}]:=If[MatchQ[q1,0]\[Or]MatchQ[q2,0]\[Or]MatchQ[q3,0],Module[{qq1,qq2,qq3},
Limit[(#1 qq1^#2[[1]] qq2^#2[[2]] qq3^#2[[3]]),{qq1->q1,qq2->q2,qq3->q3}]
],
(#1 q1^#2[[1]] q2^#2[[2]] q3^#2[[3]])
]&@@@{coeffs}//Total
MultArrayToPoly[{q0_,qs__}]:=Function[{mult,list},If[Length@list==0,Nothing,{
mult,
q0 CoeffArrayToPoly[{qs},list]//Simplify
}]]


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

(* or of the form of triples
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
If[Length@First@p!=4,Throw["sample points have to be 4-tuples"]];

(* create list of 4-tuples *)
N@List@@p
],
pType,

\[Lambda],\[Lambda]a,\[Lambda]m,\[Lambda]ma
},


If[kEnv>kSys\[Or]kSys<1\[Or]kTot<2,Throw["invalid params"]];
pType=If[MatchQ[p,_Numeric],
Numeric@@(({#[[2]]/#[[1]],#[[3]]/#[[1]],#[[4]]/#[[1]]})&/@pp),
If[Length@First@p==1,"Mono","Multi"] (* this at the moment always evaluates to Multi *)
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


(* ::Input::Initialization:: *)
Clear[ShannonEntropyMult,CIMult]
ShannonEntropyMult[l_List]:=l/.Rule[0,_]:>0/.Rule[term_,mult_]:>-mult*term Log2[term]//Total
CIMult[{\[Lambda]_,\[Lambda]a_}]:=(ShannonEntropyMult[\[Lambda]a]-ShannonEntropyMult[\[Lambda]])/Log2@Total@(Last/@List@@\[Lambda]a)
CIMult[Precalculated[\[Lambda]_->_,\[Lambda]a_->multa_]]:=(\[Lambda]a-\[Lambda])/Log2@multa


(* ::Input::Initialization:: *)
(* finds zero crossings in a list *)
Clear[ZeroCrossings]
ZeroCrossings[l_List]:=Module[{t,u,v},
t={Sign[l],Range[Length[l]]}//Transpose;
u=Select[t,First[#]!=0&];
v=SplitBy[u,First]; 
{Most[Max[#[[All,2]]]&/@v],Rest[Min[#[[All,2]]]&/@v]}//Transpose
]


(* ::Input::Initialization:: *)
Clear[CIThreshold]
CIThreshold[pa_List,q_Exact,compile_:False]:=Module[{
targetFunction,
var=First@q
},
Assert[Length@q==1\[And]var==p];

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

(* this cannot be exact, of course; for now we just numerically interpolate and find roots there; alternatively use ZeroCrossings above *)
CIThreshold[pa_Precalculated,p_Numeric,range_:All,mode_:"FindRoot"]:=Module[{
data=(CIMult@pa)[[range]],
pts=(List@@p)[[range]],
crossings,dataF,ptsF,x
},

If[mode=="ZeroCrossings",
crossings=ZeroCrossings[data];
Between[{pts[[#1]],pts[[#2]]}]&@@@crossings,

(* else *)
dataF=Interpolation[data,InterpolationOrder->2];
ptsF=ListInterpolation[#,InterpolationOrder->2]&/@Transpose[pts];

With[{xx=x/.FindRoot[dataF[x],
{x,1,Length@data,1,Length@data},
Method->Automatic,
AccuracyGoal->6,
MaxIterations->1000,
WorkingPrecision->30
]},
#[xx]&/@ptsF
]
]
]

CIThreshold[g_Graph,kSys_Integer,p_,symmetricSolver_:True,extraStabilizedVertices_List:{}]:=With[{
pa=PauliActionCC[g,kSys,p,symmetricSolver,extraStabilizedVertices]
},
CIThreshold[pa,p]
]


(* ::Input::Initialization:: *)
Clear[GraphHash]
GraphHash[g_Graph,kSys_Integer]:=StringTrim@ExportString[ExportString[CanonicalizeByOrbits[g,kSys],"Graph6"],"Base64"]//StringReplace[#,{"/"->"-"}]&


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


(* ::Input::Initialization:: *)
Clear[All2LevelGraphs]
(* all 2 level graphs with vertex count < n *)
All2LevelGraphs[n_,segmentSizes_List:{}]:=With[{
partitions=Cases[
IntegerPartitions[n-2,{2,n},If[Length@segmentSizes==0,Range[n],segmentSizes]],
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


(* ::Input::Initialization:: *)
(* get precomputed sample ranges *)
Get["p_values.m",Path->CurrentDir]


(* ::Input::Initialization:: *)
Clear[MultiResolutionSamples]
MultiResolutionSamples[chan_,threshold_,n_:10]:=MultiResolutionSamples[chan]=With[{
coarse=SAMPLES[chan][10],
fine=SAMPLES[chan][100],
ultra=SAMPLES[chan][1000],
ludicrous=SAMPLES[chan][10000]
},
Sort@Union[
coarse,
Nearest[fine,threshold,n],
Nearest[ultra,threshold,n],
Nearest[ludicrous,threshold,n]
]
]



