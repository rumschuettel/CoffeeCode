(* ::Package:: *)

(* ::Input::Initialization:: *)
Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]


(* ::Input::Initialization:: *)
CCRELEASEPATHFULLSOLVER=ParentDirectory[CurrentDir]<>"/build/release_full/";
CCRELEASEPATHFULLSOLVERARBCH=ParentDirectory[CurrentDir]<>"/build/release_full_arbch/";

PARALLELBUILDPATH=ParentDirectory[CurrentDir]<>"/build/";
PARALLELBUILDSETUPSCRIPT=CurrentDir<>"/make-build-dir.sh";

RESULTSPATH="/results/jo/";


Print["Path to precompiled full solvers: "<>CCRELEASEPATHFULLSOLVER];
Print["Path to precompiled full solvers for all channels: "<>CCRELEASEPATHFULLSOLVERARBCH];
Print["Path to create parallel solver subdirectories: "<>PARALLELBUILDPATH];
Print["Script how to setup new parallel solver subdirectories: "<>PARALLELBUILDSETUPSCRIPT];
Print["Results path, can be used as RESULTSPATH (can be ignored): "<>RESULTSPATH];
