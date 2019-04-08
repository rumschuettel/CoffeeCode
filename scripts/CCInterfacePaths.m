(* ::Package:: *)

(* ::Input::Initialization:: *)
Clear[CurrentDir]
CurrentDir:=If[$InputFileName!="",DirectoryName[$InputFileName],NotebookDirectory[]]


(* ::Input::Initialization:: *)
CCRELEASEPATHFULLSOLVER=ParentDirectory[CurrentDir]<>"/build/release_full_arbch/";

PARALLELBUILDPATH=ParentDirectory[CurrentDir]<>"/build/";
PARALLELBUILDSETUPSCRIPT=CurrentDir<>"/make-build-dir.sh";
SAGEPYTHONSCRIPT=CurrentDir<>"/sage-python.sh";

RESULTSPATH=$HomeDirectory<>"/Desktop/GraphStateResults/";


Print["Path to precompiled full solvers: "<>CCRELEASEPATHFULLSOLVER];
Print["Path to create parallel solver subdirectories: "<>PARALLELBUILDPATH];
Print["Script how to setup new parallel solver subdirectories: "<>PARALLELBUILDSETUPSCRIPT];
Print["SAGE launcher script: "<>SAGEPYTHONSCRIPT];
Print["Results path, can be used as RESULTSPATH (can be ignored): "<>RESULTSPATH];
