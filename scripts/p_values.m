(* ::Package:: *)

(* ::Text:: *)
(*Values where to sample the various channels on.*)
(*Each has to be a list of either 1- or 3- tuples; p0 will automatically be calculated.*)
(*If given a list of 1-tuples, CoffeeCode will assume that p1=p2=p3=p/3, and p0=1-p*)
(*If given a list of 3-tuples, CoffeeCode will assume {1-p1-p2-p3,p1,p2,p3}.*)
(*This can be assessed in PauliActionCC.*)
(*Note that numbers have to have head real, so be careful not to use 0 but 0.*)


(* ::Input:: *)
(*depolarizing channel : (1 - p, p/3, p/3, p/3)*)
(*threshold upper bound (antideg.) : p = 1/4*)
(*hashing point : p = 0.1893   *)


(* ::Input::Initialization:: *)
SAMPLES["depolarizing"][res_:1000]:=SAMPLES["depolarizing"][res]=Table[With[{p=0.185+(0.2-0.185)*i/res},{p/3,p/3,p/3}],{i,1,res}];


(* ::Input:: *)
(*BB84 channel: ( (1-p)^2, p-p^2,p^2,p-p^2)*)
(*threshold upper bound (antideg.): p = 0.146447*)
(*hashing point: p = 0.1100*)


(* ::Input::Initialization:: *)
SAMPLES["BB84"][res_:1000]:=SAMPLES["BB84"][res]=Table[With[{p=0.107+(0.125-0.107)*i/res},{p-p^2,p^2,p-p^2}],{i,1,res}];


(* ::Input:: *)
(*Two-Pauli channel: (1-p,p/2,0,p/2)*)
(*threshold upper bound (antideg.): p = 1/3*)
(*hashing point: p = 0.1135*)


(* ::Input::Initialization:: *)
SAMPLES["2Pauli"][res_:1000]:=SAMPLES["2Pauli"][res]=Table[With[{p=0.11+(0.25-0.11)*i/res},{p/2,0.,p/2}],{i,1,res}];
