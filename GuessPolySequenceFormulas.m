(*************************************************************************)
(*************************************************************************)
(********************* : GuessPolySequenceFormulas.m: ********************) 
(*************************************************************************)
(*************************************************************************)

BeginPackage["GuessPolySequenceFormulas`", {"GuessSequenceData`"}]

(*************************************************************************)
(**** : Package revision metadata information:                        ****)
(*************************************************************************)

PackageName = "GuessPolySequenceFormulas.m";
PackageVersion = "2014.04.28-v5";
PackageAuthor = "Maxie D. Schmidt";
PackageShortDesc = 
     "Package routines and utilities for intelligent guessing of " <> "\n" <> 
     "univariate polynomial sequence formulas."; 

(*************************************************************************)
(**** : Create default usage information for the package functions:   ****)
(**** : Provide detailed usage information for the package functions: ****)
(*************************************************************************)

GuessPolynomialSequence::usage = 
     "GuessPolynomialSequence[polyseq, polyVar]"; 

DiffSequenceFormula::usage = 
     "DiffSequenceFormula[polySeq, polyVar, polyFormula, StartIndex -> 1]"; 

(*************************************************************************)
(**** : Package error handling and messages:                          ****)
(*************************************************************************)

GPSFPkgMsgs::Errors::NonPolynomialInput = 
     "The input `1` is not a polynomial in `2`"; 

GPSFPkgMsgs::Warnings::InsufficientSeqElements = 
     "Only `1` sequence elements passed as input, " <> 
     "possibly an insufficient number to guess formulas " <> 
     "(consider passing \[GreaterEqual] `2` polynomials as input)"; 
GPSFPkgMsgs::Warnings::InsufficientFactorData = 
     "Only `1` rows of coefficient factor data specified, " <> 
     "possibly an insufficient number to correctly process formulas " <> 
     "(consider setting TriangularSequenceNumRows" <> 
     "\[LongRightArrow]`2` or more)"; 

GPSFPkgMsgs::BuildLocalSequenceData::SeqFactorsList = 
     "Invalid list of sequence factors `1`"; 
GPSFPkgMsgs::BuildLocalSequenceData::InvalidSeqID = 
     "Invalid sequence ID `1`"; 
GPSFPkgMsgs::PackageMultipleTSeqData::ListLength = 
     "Upper and Lower list lengths do not match (`1` / `2`)"; 
GPSFPkgMsgs::GuessMultipleFactorPolySequence::PolySeqDegree = 
     "Unable to compute the polynomial sequence degrees"; 

(*************************************************************************)
(**** : Package configuration and global settings: ****)
(*************************************************************************)

Begin["`PkgConfig`"]

Debugging = True;
KeepStatusMessages = True;
KeepStatusMatchMessages = True;
AllowSymbolicData = False; 

End[] (* PkgConfig *)

Begin["`Private`"]

(*************************************************************************)
(**** : Misc utilities for the package: ****)
(*************************************************************************)

PrintDebug[args___] := Module[{debugFlag, printPrefix}, 
     debugFlag = GuessPolySequenceFormulas`PkgConfig`Debugging; 
     If[debugFlag, 
          printPrefix = Style["DEBUGGING: ", Red, Bold];
          Print[printPrefix, args];
     ];
]

PrintStatus[msgs___] := Module[{keepMsgFlag, printPrefix}, 
     keepMsgFlag = GuessPolySequenceFormulas`PkgConfig`KeepStatusMessages; 
     If[keepMsgFlag, 
          printPrefix = Style["STATUS: ", Blue, Bold];
          Print[printPrefix, msgs];
     ];
]

PrintMatch[msgs___] := Module[{keepMsgFlag, printPrefix}, 
     keepMsgFlag = GuessPolySequenceFormulas`PkgConfig`KeepStatusMatchMessages; 
     If[keepMsgFlag, 
          printPrefix = Style["MATCH_DATA: ", Cyan, Bold, Underlined];
          Print[printPrefix, msgs];
     ];
]

PrintMatchInfo[msgs___] := Module[{keepMsgFlag, printPrefix}, 
     keepMsgFlag = GuessPolySequenceFormulas`PkgConfig`KeepStatusMatchMessages; 
     If[keepMsgFlag, 
          printPrefix = Style["MATCH_INFO: ", Green, Bold, Underlined];
          Print[printPrefix, msgs];
     ];
]

PrintVariable[var_] := Module[{keepMsgFlag, printPrefix}, 
     keepMsgFlag = GuessPolySequenceFormulas`PkgConfig`KeepStatusMatchMessages; 
     If[keepMsgFlag, 
          printPrefix = Style["VAR_VALUE: ", Orange, Bold, Underlined];
          Print[printPrefix, SymbolName[var], ": ", var];
     ];
]

ConfigDebuggingMessages[setOnOff_:{False,False,False}] := Block[{}, 
     GuessPolySequenceFormulas`PkgConfig`Debugging = setOnOff[[1]];
     GuessPolySequenceFormulas`PkgConfig`KeepStatusMessages = setOnOff[[2]];
     GuessPolySequenceFormulas`PkgConfig`KeepStatusMatchMessages = setOnOff[[3]];
];

EnableDebugging[] := ConfigDebuggingMessages[ConstantArray[True, 3]];
DisableDebugging[] := ConfigDebuggingMessages[ConstantArray[False, 3]];

(*% \label{scode_linelbl_GPSF.m_StartRuntimeCfgOptionsListing_v1} %*)
DefaultPkgConfigOptions = {
     StartIndex -> 1, 
     IndexMultiples -> {0, 1}, 
     IndexOffsetPairs -> Null, 
     SequenceFactors -> {"S1"},     (* See GuessSequenceData.m *)
     LimitFormulaCount -> Infinity, (* Limit the number of returned formulas *) 
     ClearLocalSequenceData -> True, 
     FactorFunction -> Null,        (* Null for local handling functions, or 
                                       possibility for a user-defined function *)
     UserGuessFunction -> DefaultUserGuessFn, 
     FSFFunction -> FindSequenceFunction, 
     FSFFnExpectedStartIndex -> 1, 
     FSFImposeTimingConstraints -> False, 
     FSFTimingConstraintSeconds -> 1.0, 
     CheckForConstantSequences -> True, 
     RequireMatchingULIndexFormula -> True, 
     RequireMatchingRemSequenceFormula -> False, 
     IssueAllWarnings -> False, 
     AllowSymbolicData -> False, 
     ProcessRemSequenceFormulas -> True, 
     ReverseRemSeqFormulaIndex -> False, 
     TriangularSequenceNumRows -> 12, 
     ExtractSequenceDataFunction -> DefaultExtractSequenceDataFn, 
     ExtractUIIndexDataFunction -> DefaultExtractUIIndexDataFn, 
     ReturnNullSequenceType -> Null, 
     RemSeqProcFnFlag -> Null, 
     GenerateNotebookOutput -> True, 
     GeneratePlaintextOnly -> False, 
     ReturnFormulasOnly -> True, 
     ReturnFullFormulaData -> False, 
     DisplayVars -> {j, i}, 
     PrintFormulas -> True, 
     PrintLaTeXFormulas -> True, 
     PrintPkgDebuggingMsgs -> {True, True, Print}, (* {T/F, VerboseQ, PrintFn} *) 
     PrintPkgStatusMsgs -> {True, True, Print},    (* {T/F, VerboseQ, PrintFn} *) 
     PrintPkgMatchMsgs -> {True, True, Print},     (* {T/F, VerboseQ, PrintFn} *)
     PrintRuntimeTimingMsgs -> {True, False, PrintTemporary}, 
     SimplifyFormulas -> True, 
     SimplifyFns -> {Simplify}, 
     FullSimplifyFormulas -> False, 
     FullSimplifyFns -> {PowerExpand, FunctionExpand, FullSimplify}, 
     DiffPolySequenceFormulas -> True, 
     EnableDebugging -> False, 
     AbortOnSymmetricRemSequence -> True, 
     RequireSeqDegreeFormula -> True, (* Otherwise, specify C[i][#1] *)
     ProcessPolyDataFn -> Null 
}; 

GetDefaultConfigOption[option_] := (option /. DefaultPkgConfigOptions); 

(*************************************************************************)
(**** : Configuration of (semi) static data:                          ****)
(*************************************************************************)

PkgNumSequenceFactors = 0;
PkgSequenceData = {}; 

GetTriangularSequenceData[numRows_, startRowIndex_, 
                          seqGenFn_, rowRangeFn_] := 
Module[{seqData, seqPkgFn, row, rowIndex, colIndexRange, nextRow}, 

     seqData = {}; 
     seqPkgFn = PackageTriangularSequenceData[#1, #2, seqGenFn[#1, #2]]&;
     For[row = 0, row < numRows, row++, 
          rowIndex = row + startRowIndex; 
          colIndexRange = Range[##]& @@ rowRangeFn[rowIndex]; 
          nextRow = Map[seqPkgFn[rowIndex, #1]&, colIndexRange];
          seqData = Union[seqData, nextRow];
     ];

     PrintDebug["GetTSeqData seqData: ", seqData]; 
     Return[seqData];

]

Options[BuildLocalSequenceData] = DefaultPkgConfigOptions; 
BuildLocalSequenceData[cfgOptions : OptionsPattern[]] := 
Module[{sequenceIDs, numRows, sindex, seqID, seqMetaData, startRowIndex, 
        seqGenFn, rowRangeFn, curSeqData}, 

     sequenceIDs = OptionValue[SequenceFactors]; 
     If[(Head[sequenceIDs] != List) || (Length[sequenceIDs] == 0), 
          Message[GPSFPkgMsgs::BuildLocalSequenceData::SeqFactorsList, 
                  sequenceIDs]; 
          Return[Null]; 
     ]; 

     (** : Reset the local sequence data storage array: **) 
     PkgNumSequenceFactors = Length[sequenceIDs]; 
     PkgSequenceData = {}; 
     
     numRows = OptionValue[TriangularSequenceNumRows]; 
     For[sindex = 1, sindex <= Length[sequenceIDs], sindex++, 

          (** : Uses the local sequences package data / implementations: **)
          seqID = sequenceIDs[[sindex]]; 
          seqMetaData = QuerySequenceMetaData[seqID]; 
          If[Length[seqMetaData] == 0, 
               Message[GPSFPkgMsgs::BuildLocalSequenceData::InvalidSeqID, 
                       seqID]; 
               PkgSequenceData = {}; 
               PkgNumSequenceFactors = 0; 
               Return[False]; 
          ]; 
          PrintDebug[seqKeys, " ", seqKey, " ", seqMetaData];

          startRowIndex = (StartRowIndex) /. seqMetaData; 
          seqGenFn = (GeneratorFunction) /. seqMetaData; 
          rowRangeFn = (RowRangeFunction) /. seqMetaData; 

          curSeqData = GetTriangularSequenceData[numRows, startRowIndex, 
                                                 seqGenFn, rowRangeFn]; 
          PkgSequenceData = Append[PkgSequenceData, curSeqData]; 

     ]; 

     Return[True]; 

]

FactorIntegerBySequences[int_] := 
Module[{firstSeqFactorData, numSequences, fullFactorData, sindex, curSeqData, 
        nextFullFactorData, eindex, prevFactors, prevInt, divCondFn, 
        curSeqFactorData, nextFactorDataFn, nextFactorData}, 

     firstSeqFactorData = FactorIntegerBySequence[int]; 
     numSequences = PkgNumSequenceFactors; 

     (** : Add brackets around the first sequence index data: **) 
     (*fullFactorData = Map[{ {ExtractTSeqULIndexData[#1]}, 
                            ExtractTSeqValue[#1] }&, firstSeqFactorData]; *)
     fullFactorData = firstSeqFactorData; 

     For[sindex = 2, sindex <= numSequences, sindex++, 

          curSeqData = PkgSequenceData[[sindex]]; 
          nextFullFactorData = {}; 
          For[eindex = 1, eindex <= Length[fullFactorData], eindex++, 

               prevFactors = ExtractTSeqULIndexData[ fullFactorData[[eindex]] ]; 
               prevInt = ExtractTSeqValue[ fullFactorData[[eindex]] ]; 
               divCondFn = Divisible[prevInt, ExtractTSeqValue[#1]]&;
               curSeqFactorData = Select[curSeqData, divCondFn];
               localReplaceFn = ReplacePart[#1, 
                                2 -> (prevInt / ExtractTSeqValue[#1])]&;
               curSeqFactorData = Map[localReplaceFn, curSeqFactorData]; 
               
               nextFactorDataFn = 
                    { Append[prevFactors, ExtractTSeqULIndexData[#1]], 
                      ExtractTSeqValue[#1] }&; 
               nextFactorData = Map[nextFactorDataFn, curSeqFactorData]; 
               nextFullFactorData = Union[nextFullFactorData, 
                                          nextFactorData]; 

          ]; 
          fullFactorData = nextFullFactorData; 

     ]; 
     
     Return[fullFactorData]; 

]

FactorIntegerBySequence[int_] := 
Module[{seqData, leadingIntCoeff, allowSymbolicData, divCondFn, 
        factorData, replaceFn, rFactorData}, 

     seqData = PkgSequenceData[[1]]; 

     (** : Handle symbolic data by extracting the integer coefficient: **) 
     leadingIntCoeff = int; 
     allowSymbolicData = GuessPolySequenceFormulas`PkgConfig`AllowSymbolicData; 
     If[allowSymbolicData && !IntegerQ[int] && 
        (Length[int] > 0) && (Head[int[[1]]] == Integer), 
          leadingIntCoeff = int[[1]]; 
     ]; 

     (** : First determine which sequence elements divide the integer: **) 
     divCondFn = Divisible[leadingIntCoeff, ExtractTSeqValue[#1]]&;
     factorData = Select[seqData, divCondFn];

     (** : Then replace the sequence elements with the remainder of int: **)
     replaceFn = ReplacePart[#1, 2 -> (int / ExtractTSeqValue[#1])]&;
     
     replaceFn = {{ExtractTSeqULIndexData[#1]}, int / ExtractTSeqValue[#1]}&; 
     rFactorData = Map[replaceFn, factorData];
     Return[rFactorData];

]

(**** : Helper routines for handling triangular sequence data: ****)
PackageTriangularSequenceData[upperi_, loweri_, value_] := 
Module[{indexData, valueData, rData}, 
     
     indexData = {upperi, loweri};
     valueData = value;
     rData = {indexData, valueData};
     Return[rData];

]

ExtractTSeqULIndexData[edata_] := Module[{}, 
     Return[ edata[[1]] ];
]

ExtractTSeqUpperIndex[edata_] := Module[{}, 
     Return[ edata[[1]][[1]] ];
]

ExtractTSeqLowerIndex[edata_] := Module[{}, 
     Return[ edata[[1]][[2]] ];
]

ExtractTSeqValue[edata_] := Module[{}, 
     Return[ edata[[2]] ];
]

PackageMultipleTSeqData[uiList_, liList_, remValue_] := 
Module[{ulIndexData, rData}, 
     
     If[Length[uiList] != Length[liList], 
          Message[GPSFPkgMsgs::PackageMultipleTSeqData::ListLength, 
                  uiList, liList]; 
          Return[{ {}, remValue}];
          Return[Null]; 
     ];
     
     If[Length[uiList] == 0, (* Return empty list for the sequence data *) 
          Return[{ {}, remValue}]; 
     ]; 

     ulIndexData = MapIndexed[{uiList[[ First[#2] ]], liList[[ First[#2] ]]}&, 
                              Range[Length[uiList]]]; 
     
     rData = {ulIndexData, remValue}; 
     
     Return[rData]; 

]

ExtractMultipleTSeqUpperIndices[edata_] := 
Module[{ulIndexData, uIndexList}, 
     
     ulIndexData = ExtractTSeqULIndexData[edata]; 
     uIndexList = Map[(#1)[[1]]&, ulIndexData]; 
     Return[uIndexList];

]

ExtractMultipleTSeqLowerIndices[edata_] := 
Module[{ulIndexData, lIndexList}, 
     
     ulIndexData = ExtractTSeqULIndexData[edata]; 
     lIndexList = Map[(#1)[[2]]&, ulIndexData]; 
     Return[lIndexList];

]

ExtractRemainingValue[edata_] := ExtractTSeqValue[edata]

DefaultUserGuessFn[polyIndex_, sumIndex_] := 1;

IsEmptySet[set_] := (Length[set] == 0);
IsConstantSequence[seq_] := Module[{}, 
     Return[!IsEmptySet[seq] && (Length[Tally[seq]] == 1)]; 
]

ConstantFn[constant_, args___] := constant;
GetConstantFn[constant_] := ConstantFn[constant, ##]&
GetNullValuedFn[] := GetConstantFn[Null]
SetAttributes[ConstantFn, Constant];

GetFixedSequenceDataFn[fixedSeqData_, startIndex_:0, 
                       genParamIndex_:0] := Which[ 
     #1 < startIndex, SEQ[genParamIndex][#1], 
     #1 >= startIndex + Length[fixedSeqData], SEQ[genParamIndex][#1], 
     _, fixedSeqData[[1 - startIndex + #1]]
]; 

PackagePolyIndexMatchData[uiIndexData_, seqData_] := {uiIndexData, seqData};
DefaultExtractUIIndexDataFn[elem_] := elem[[1]];
DefaultExtractSequenceDataFn[elem_] := elem[[2]]; 

Options[DefaultPackageFormulaDataFn] = DefaultPkgConfigOptions; 
DefaultPackageFormulaDataFn[ulIndexFnData_, remSeqFnData_, 
                            cfgOpts : OptionsPattern[]] := 
Module[{}, 

     Return[Null]; 

]; 

Options[DiffSequenceFormula] = DefaultPkgConfigOptions; 
DiffSequenceFormula[polySeq_, polyVar_, polyFormula_, 
                    cfgOpts : OptionsPattern[]] := 
Module[{seqStartIndex, indexAdjust, seqDiff, diffFn, polyDiffs}, 

     seqStartIndex = OptionValue[StartIndex]; 
     indexAdjust = seqStartIndex - 1; 
     si = Unique[si]; 
     seqDiff = (#1 - polyFormula[#2 + indexAdjust, si, polyVar])&; 
     diffFn = seqDiff[#1, #2]&; 
     polyDiffs = MapIndexed[(Simplify[Expand[diffFn[#1, First[#2]]]] === 0)&, 
                            polySeq]; 
     Return[polyDiffs]; 

] 

(**** : Begin ComputeSequenceFormula[...] helper functions : ****)

(*formatReturnDataFn[matchQ_, seqFn_] := {matchQ, seqFn}; *)
formatReturnDataFn[matchQ_, seqFn_] := Module[{}, 
     Return[{matchQ, seqFn}]; 
]

(**** : End ComputeSequenceFormula[...] helper functions : ****)

Options[ComputeSequenceFormula] = DefaultPkgConfigOptions; 
ComputeSequenceFormula[seqData_, startIndex_, 
                       cfgOpts : OptionsPattern[]] := 
Module[{FSFFn, expectedStartIndex, checkConstantSequences, initFormula, 
        constFn, indexAdjust, seqFormulaFn, simplifyQ, fullSimplifyQ, 
        simplifyFns, simpFormula, simpSeqFormulaFn}, 

     FSFFn = OptionValue[FSFFunction];
     expectedStartIndex = OptionValue[FSFFnExpectedStartIndex];
     checkConstantSequences = OptionValue[CheckForConstantSequences]; 
     
     PrintDebug["seqData: ", seqData]; 
     PrintDebug["startIndex: ", startIndex]; 
     PrintDebug["CheckConst = ", checkConstantSequences, "; expected = ", 
                expectedStartIndex];

     (** : attempt to compute a closed-form formula for the sequence: **)
     initFormula = FSFFn[seqData];

     (** : Return NULL if FSFFn[...] cannot find a formula expression: **) 
     If[Head[initFormula] == FSFFn, 
          If[checkConstantSequences && IsConstantSequence[seqData], 
               constFn = GetConstantFn[ seqData[[1]] ]; 
               Return[formatReturnDataFn[True, constFn]]; 
          ];
          Return[formatReturnDataFn[False, GetConstFn[Null]]];
     ];

     indexAdjust = startIndex - expectedStartIndex;
     PrintDebug["indexAdjust: ", indexAdjust];
     
     seqFormulaFn = initFormula /. (# -> # - indexAdjust); 
     
     (** : Attempt to simplify the formula if the options are set : **)
     (** : Note: FullSimplifyFormulas takes precedence over the   : **) 
     (** :       default simplify option SimplifyFormulas         : **)
     simplifyQ = OptionValue[SimplifyFormulas]; 
     fullSimplifyQ = OptionValue[FullSimplifyFormulas]; 
     simplifyFns = Which[fullSimplifyQ, OptionValue[FullSimplifyFns], 
                         simplifyQ, OptionValue[SimplifyFns], _, {}]; 
 
     If[(!simplifyQ && !fullSimplifyQ) || IsEmptySet[simplifyFns], 
          Return[formatReturnDataFn[True, seqFormulaFn]]; 
     ]; 
    
     (** : Apply the simplification functions in default order: **)
     PrintStatus["Pre seqFormulaFn: ", seqFormulaFn]; 
     xi = Unique[xvar]; 
     simpFormula = Last[ComposeList[simplifyFns, seqFormulaFn[xi]]]; 
     simpSeqFormulaFn = Evaluate[simpFormula /. (xi -> #1)]&;
     PrintStatus["Post seqFormulaFn: ", simpSeqFormulaFn]; 

     Return[formatReturnDataFn[True, simpSeqFormulaFn]];

]

(**** : Perform pre-processing of the input polynomial sequence terms : ****)
PreProcessPolynomialData[inputPolyForm_, polyIndex_, polyVar_, userGuessFn_, 
                         cfgOpts : OptionsPattern[], 
                         processOpts : OptionsPattern[]] := 
Module[{pxvar, cfList, processCoeffFn, modcfList, modPolySum, modPoly, 
        modTerms, modPolyMetadata, rData}, 

     (** : Perform basic error checking : **)
     PrintStatus["In pre-process ... PolyIndex = ", polyIndex];

     (** : Adjust / scale the original coefficients: **) 
     cfList = CoefficientList[inputPolyForm, polyVar]; 
     processCoeffFn = Cancel[(#1) / userGuessFn[polyIndex, #2]]&;
     modcfList = MapIndexed[processCoeffFn[#1, First[#2] - 1]&, cfList];

     (** : Get the new, modified polynomial function : **)
     modPolySum = MapIndexed[(#1) * Power[pxvar, First[#2] - 1]&, modcfList];
     modPoly = Function[pxvar, modPolySum][polyVar]; 
     modTerms = MapIndexed[(#1) * Power[polyVar, First[#2] - 1]&, modcfList];
     modPoly = Plus @@ modTerms;
     modPolyMetadata = Null; 

     (** : Setup the returned polynomial / accounting information: **) 
     rData = {ModPolyFn -> modPoly, 
              ModPolyMetadata -> modPolyMetadata}; 
     Return[rData]; 
     
]

Options[ProcessPolynomialWrapper] = DefaultPkgConfigOptions; 
ProcessPolynomialWrapper[poly_, polyIndex_, polyVar_, indexOffsets_, 
                         factorFn_, cfg : OptionsPattern[]] := 
     Module[{cflist, cfFactorData, processPolyDataFn, matchData}, 

     (** : Perform error checking: **)
     PrintStatus["In Process Fn ...  PolyIndex = ", polyIndex, " p=", poly, " v=", polyVar];
     
     (** : Get the coefficient factor data (apply user guess function): **) 
     cflist = CoefficientList[poly, polyVar]; 
     cfFactorData = Map[factorFn, CoefficientList[poly, polyVar]]; 

     (** : Compute the possible matches by calling the subroutine: **) 
     processPolyDataFn = OptionValue[ProcessPolyDataFn]; 
     matchData = processPolyDataFn[cfFactorData, indexOffsets, cfg]; 
     
     Return[matchData];

] 

ProcessSingleFactorPoly[poly_, polyIndex_, polyVar_, indexOffsets_, 
                        factorFn_, cfgOptions : OptionsPattern[]] := 
ProcessPolynomialWrapper[poly, polyIndex, polyVar, indexOffsets, factorFn, 
                         ProcessPolyDataFn -> ProcessSingleFactorPolyData, 
                         cfgOptions]; 

ProcessMultipleFactorPoly[poly_, polyIndex_, polyVar_, indexOffsets_, 
                          factorFn_, cfgOptions : OptionsPattern[]] := 
ProcessPolynomialWrapper[poly, polyIndex, polyVar, indexOffsets, factorFn, 
                         ProcessPolyDataFn -> ProcessMultipleFactorPolyData, 
                         cfgOptions]; 

getMultipleIndexOffsetFunction[uiList_, liList_, indexOffsetsList_, 
                               valData_] := 
Module[{uiOffsets, liOffsets, ulIndexFn, getOffsetFns, uiFns, liFns, 
        ulOffsetFn, mapListFn1, mapListFn2, mapListFn12, 
        intermedFn, ulOffsetFn2}, 

     (** : Error checking on the lengths of the input lists? : **) 

     PrintStatus["indexOffsetsList: ", indexOffsetsList]; 
     uiOffsets = Map[(#1)[[1]]&, indexOffsetsList];
     PrintStatus["uiOffsets: ", uiOffsets]; 

     liOffsets = Map[(#1)[[2]]&, indexOffsetsList];
     
     ulIndexFn[ulIndex_, ulOffset_] := (ulIndex + ulOffset * (#1))&; 
     getOffsetFns[indexList_, offsetList_] := 
          MapIndexed[ulIndexFn[indexList[[ First[#2] ]], #1]&, offsetList]; 
     
     uiFns = getOffsetFns[uiList, uiOffsets]; 
     PrintStatus["I: ", uiFns]; 
     
     liFns = getOffsetFns[liList, liOffsets]; 
     PrintStatus["II: ", liFns]; 
     
     ulOffsetFn = PackageMultipleTSeqData[uiFns, liFns, valData]; 
     PrintStatus["III: ", ulOffsetFn]; 

     (*x = Unique[xi]; 
     mapListFn1 = Map[((#1)[x])&, #1]&;
     mapListFn2 = Map[MapListFn1[#1]&, #1]&; 
     intermedFn = HoldForm[Evaluate[ReplacePart[ulOffsetFn, 
                  1 -> (mapListFn2[ ulOffsetFn[[1]] ])] /. (x -> #)]]; 
     ulOffsetFn2 = Function @@ intermedFn; *)
     
     (** : What these next few lines are actually used for in the code: **)
     (** : Should transform a list that looks like                      **)
     (** : {{{0 + 0 #1 &, 0 + 1 #1 &}}, _} or like                      **) 
     (** : {{{a + b #1 &, c + d #1 &}}, _} ...                          **)
     (** : into the correct function syntax of a list that looks like   **) 
     (** : {{{0 + 0 #1, 0 + 1 #1}}, _}& or like                         **)
     (** : {{{a + b #1, c + d #1}}, _}&                                 **)
     (*% \label{scode_linelbl_GPSF.m_With_usage_example_v1} %*)
     x = Unique[xi]; 
     mapListFn12 = Map[(Map[((#1)[x])&, #1]&)[#1]&, #1]&; 
     ulOffsetFn2 = With[{uloFn = ulOffsetFn, mlFn12 = mapListFn12, xvar = x}, 
                   Function @@
                   HoldForm[Evaluate[ReplacePart[uloFn, 
                   1 -> (mapListFn12[ ulOffsetFn[[1]] ])] /. (xvar -> #)]] ]; 
    
     PrintStatus["IV: ", ulOffsetFn2]; 
     
     Return[ulOffsetFn2]; 

]

Options[ProcessSingleFactorPolyData] = DefaultPkgConfigOptions; 
ProcessSingleFactorPolyData[polyFactorData_, indexOffsets_, 
                            cfg : OptionsPattern[]] := 
Module[{cfListFactors, firstcfFactors, matches, cfi, curFactorData, 
        uIndexList, lIndexList, ulOffsetFn, testMemberFormFn, testOffsets, 
        remTermPos, remSeqTerms, matchData}, 
 
     (*ConfigDebuggingMessages[{True, False, True}]; *)
     PrintStatus["In Second Process Fn ... "];

     (** : Later can pick the index with the fewest factor data terms ... : **) 
     cfListFactors = polyFactorData; 
     firstcfFactors = cfListFactors[[1]];
     cfListFactors = Drop[cfListFactors, 1];
     
     PrintDebug["firstcfFactors: ", firstcfFactors]; 
     PrintDebug["cfListFactors: ", cfListFactors]; 

     matches = {}; 
     For[cfi = 1, cfi <= Length[firstcfFactors], cfi++, 
                    
          curFactorData = firstcfFactors[[cfi]]; 
          PrintVariable[curFactorData]; 
          
          uIndexList = ExtractMultipleTSeqUpperIndices[curFactorData];
          PrintVariable[uIndexList]; 

          lIndexList = ExtractMultipleTSeqLowerIndices[curFactorData]; 
          PrintVariable[lIndexList]; 
          ulOffsetFn = getMultipleIndexOffsetFunction[uIndexList, lIndexList, 
                                                      indexOffsets, _];
          PrintStatus["after I ", ulOffsetFn]; 
         
          testMemberFormFn = MemberQ[#1, ulOffsetFn[ First[#2]] ]&; 
          testOffsets = MapIndexed[testMemberFormFn[#1, #2]&, cfListFactors];
          PrintStatus["after II"]; 

          If[!MemberQ[testOffsets, False], 
               
               (** : Potential match, get the leftover sequence terms: **)
               remTermPos = MapIndexed[Position[#1, ulOffsetFn[First[#2]]]&, 
                                       cfListFactors]; 
               PrintMatchInfo["remTermPos (pre flatten): ", remTermPos]; 
               remTermPos = Flatten[remTermPos, 1];
               PrintMatchInfo["remTermPos (post flatten): ", remTermPos]; 

               remSeqTerms = MapIndexed[ExtractTSeqValue[ 
                             cfListFactors[[ First[#2] ]][[#1]][[1]] ]&, 
                             remTermPos]; 
               remSeqTerms = Prepend[remSeqTerms, 
                                     ExtractTSeqValue[curFactorData]]; 

               (*matchData = {{uIndex, lIndex}, remSeqTerms};*)
               matchingULIndexData = ExtractTSeqULIndexData[ 
                    PackageMultipleTSeqData[uIndexList, lIndexList, Null] 
               ]; 
               matchData = PackagePolyIndexMatchData[matchingULIndexData, 
                                                     remSeqTerms]; 
               
               (*matchData = {{uIndex, lIndex}, remSeqTerms, Reverse[remSeqTerms]};*)
               PrintMatchInfo["matchData : ", matchData];
               matches = Append[matches, matchData];

          ];

     ];

     PrintMatchInfo["full matches list: ", matches];
     Return[matches]

]

(**** : Begin FilterByRemSequenceMatches[...] helper functions : ****)
     
getPermsMapFn[seqEntry_] := 
Module[{seqData, uiIndexPerms}, 
     
     seqData = seqEntry[[1]];
     uiIndexPerms = Tuples[seqEntry[[2]]];
     Return[Map[{#1, seqData}&, uiIndexPerms]];

] 

(**** : End FilterByRemSequenceMatches[...] helper functions : ****)

Options[FilterByRemSequenceMatches] = DefaultPkgConfigOptions; 
FilterByRemSequenceMatches[matchData_, preProcFn_, 
                           cfgOptions : OptionsPattern[]] := 
Module[{localFilterBySeqDataFn, lastIndex, lastMatchData, nullSequenceValue, 
        mdi, prevMatchData, seqLenDiff, seqEqualsFn, 
        diffPrevIndexRemSeqTermsFn, fullPermsData}, 

     (** : Returns a list of processed pairs of the form **)
     (** : {SeqData, {UI Index Match Pairs}}:            **)
     localFilterBySeqDataFn[mdIndex_, inputMatchData_] := 
     Module[{extractSeqDataFn, extractUIDataFn, seqSplitFn, splitResults, 
             procDataList, sindex, splitEntry, seqData, uiDataList}, 
          
          PrintStatus["inputMatchData[[", mdIndex, "]]: ", 
                      inputMatchData[[mdIndex]]];

          extractSeqDataFn = OptionValue[ExtractSequenceDataFunction];
          extractUIDataFn = OptionValue[ExtractUIIndexDataFunction];
          seqSplitFn = preProcFn[extractSeqDataFn[#1]]&;
          splitResults = SplitBy[inputMatchData[[mdIndex]], seqSplitFn];
          PrintDebug["localFilterBySeqDataFn[", mdIndex, "]: ", splitResults];
          
          procDataList = {};
          For[sindex = 1, sindex <= Length[splitResults], sindex++, 
               splitEntry = splitResults[[sindex]];
               seqData = seqSplitFn[splitEntry[[1]]];
               uiDataList = Map[extractUIDataFn[#1]&, splitEntry];
               PrintDebug["=> sindex = ", sindex, ": ", uiDataList];
               procDataList = Append[procDataList, {seqData, uiDataList}];
          ];
          
          PrintDebug["localFilterBySeqDataFn[", mdIndex, "]: ", procDataList];
          Return[procDataList];

     ]; (* localFilterBySeqDataFn *) 

     lastIndex = Length[matchData];
     lastMatchData = localFilterBySeqDataFn[lastIndex, matchData];
     
     (* Store a list of lists for the last arguments: *)
     lastMatchData = Map[{(#1)[[1]], {(#1)[[2]]}}&, lastMatchData]; 
     PrintStatus["At mdi index = lastIndex (", lastIndex, "): ", lastMatchData];

     nullSequenceValue = OptionValue[ReturnNullSequenceType];
     
     For[mdi = lastIndex - 1, 1 <= mdi, mdi--, 

          (** : Return the empty set if there are no current : **)
          (** : consistent sequence matches:                   **) 
          If[IsEmptySet[lastMatchData], 
               Return[{}]; 
          ]; 
          
          prevMatchData = localFilterBySeqDataFn[mdi, matchData]; 
          PrintStatus["At mdi index = ", mdi, " (I): ", prevMatchData];

          seqLenDiff = lastIndex - mdi; 
          seqEqualsFn = (Drop[#1, -seqLenDiff] == #2)&;
          diffPrevIndexRemSeqTermsFn[lastmd_] := 
          Module[{lastmdSeq, prevSeqMatch}, 
               
               PrintDebug["==> lastmd: ", lastmd];
               lastmdSeq = lastmd[[1]];
               prevSeqMatch = Select[prevMatchData, 
                                     seqEqualsFn[lastmdSeq, (#)[[1]] ]&];
               If[Length[prevSeqMatch] == 0, 
                    Return[nullSequenceValue];
               ];
               
               prevSeqMatch = prevSeqMatch[[1]]; (* Select returns a list *)
               PrintDebug["==> prevSeqMatch: ", prevSeqMatch];
               Return[{lastmdSeq, Prepend[lastmd[[2]], prevSeqMatch[[2]] ]}];
          
          ];
          lastMatchData = Map[diffPrevIndexRemSeqTermsFn[#1]&, lastMatchData];
          lastMatchData = Cases[lastMatchData, Except[nullSequenceValue]]; 
          PrintStatus["At mdi index = ", mdi, " (II): ", lastMatchData];

     ]; (* for mdi *) 

     PrintStatus[" ==== After mdi loop ==== "];
     PrintDebug["lastMatchData: ", lastMatchData];

     (** : Generate all possible permutations of the form :  **)
     (** : {{UI Index Data Pairs}, Processed Sequence Data}: **) 
     PrintDebug["fullPermsData (pre): ", Map[getPermsMapFn[#1]&, lastMatchData]]; 
     fullPermsData = Flatten[Map[getPermsMapFn[#1]&, lastMatchData], 1];
     PrintDebug["fullPermsData (post): ", fullPermsData]; 

     Return[fullPermsData];

]

(**** : Begin VerifyFormulaMatches[...] helper functions : ****)
     
Options[computeIndexFormulaFn] = DefaultPkgConfigOptions; 
computeIndexFormulaFn[seqData_, indexOffset_, 
                      cfgOptions : OptionsPattern[]] := 
Module[{seqStartIndex, formulaFnData, nullFn, formulaFnj, offsetFn, 
        fullIndexFn}, 

     seqStartIndex = OptionValue[StartIndex];
     
     formulaFnData = ComputeSequenceFormula[seqData, seqStartIndex, 
                                            cfgOptions]; 
     If[!formulaFnData[[1]], 
          nullFn = formulaFnData[[2]];
          Return[{False, nullFn, nullFn, nullFn}];
     ];

     formulaFnj = formulaFnData[[2]]; 
     offsetFn = (indexOffset * (#))&;
     fullIndexFn = (offsetFn[#2] + formulaFnj[#1])&; 
     
     (** : Later, change these to transformation rules, or give an : **) 
     (** : extract data function for the indices ... :               **) 
     formulaFnData = {True, formulaFnj, offsetFn, fullIndexFn}; 
     
     Return[formulaFnData]; 

] 

Options[computeRemSeqFormulaFn] = DefaultPkgConfigOptions; 
computeRemSeqFormulaFn[rseqData_, 
                       cfgOptions : OptionsPattern[]] := 
Module[{requireRemSeqFormula, processRemSeqFormulas, reverseRemSeqIndices, 
        rsFormulaData, rsInitFormula, rsFormula}, 

     requireRemSeqFormula = OptionValue[RequireMatchingRemSequenceFormula];
     processRemSeqFormulas = OptionValue[ProcessRemSequenceFormulas];
     reverseRemSeqIndices = OptionValue[ReverseRemSeqFormulaIndex];

     If[requireRemSeqFormula || processRemSeqFormulas, 

          rsFormulaData = ComputeSequenceFormula[rseqData, 0, cfgOptions];
          rsInitFormula = rsFormulaData[[2]];
          rsFormula = rsInitFormula; 
          If[!reverseRemSeqIndices, 
               rsFormula = (rsInitFormula[#2])&,     (* index i *)
               rsFormula = (rsInitFormula[#1 - #2])& (* index j-i *)
          ]; 
          
          Return[{rsFormulaData[[1]], rsFormula, rseqData}];
     
     ]; 
     
     Return[{True, GetFixedSequenceDataFn[rseqData, 0, #1]&, rseqData}]; 

] 
    
getUIIndexSeqFn[cfData_, indexPos_, factorPos_] := 
Module[{mapExtractFn}, 
     
     mapExtractFn = (#1)[[factorPos]][[indexPos]]&; 
     Return[Map[mapExtractFn[#1]&, cfData]];

]

(**** : End VerifyFormulaMatches[...] helper functions : ****)

Options[VerifyFormulaMatches] = DefaultPkgConfigOptions; 
VerifyFormulaMatches[inputMatchData_, indexOffsets_, 
                     cfgOptions : OptionsPattern[]] := 
Module[{requireULIndexFormula, requireRemSeqFormula, noProcessSymmetricRemSeq, 
        reverseRemSeqIndices, seqFactorIndex, uiIndexOffset, liIndexOffset, 
        formulaMatches, mresult, matchResultData, coeffData, 
        upperIndexSeq, uiFn, lowerIndexSeq, liFn, remSeq, rsFn, 
        formulaData, numSeqFactors, uiFnsList, liFnsList, 
        breakOnFormulaError}, 

     PrintDebug["inputMatchData: ", inputMatchData];

     (** : The primary settings for verifying the index and : **)
     (** : sequence formulas are matched:                     **)
     requireULIndexFormula = OptionValue[RequireMatchingULIndexFormula];
     requireRemSeqFormula = OptionValue[RequireMatchingRemSequenceFormula];
     noProcessSymmetricRemSeq = OptionValue[AbortOnSymmetricRemSequence]; 
     
     (** : Configure whether the remaining sequence index input is : **) 
     (** : i or j-i :                                                **)
     reverseRemSeqIndices = OptionValue[ReverseRemSeqFormulaIndex];
     SetOptions[computeRemSeqFormulaFn, 
                ReverseRemSeqFormulaIndex -> reverseRemSeqIndices]; 

     (** : Setup other configuration setting options: **)
     (*PackageFormulaFn = OptionValue[PackageIntermediateFormulaDataFunction]; *)

     (** : Setup sequence offset and position information: **) 
     numSeqFactors = Length[OptionValue[SequenceFactors]]; 

     (** : Find which of the possible matches give consistent formulas: **)
     formulaMatches = {}; 
     For[mresult = 1, mresult <= Length[inputMatchData], mresult++, 

          matchResultData = inputMatchData[[mresult]];
          PrintStatus["matchResultData (mresult = ", mresult, "): ", 
                      matchResultData];  
                    
          (** : Process the remaining sequence formula for this match: **)
          remSeq = matchResultData[[2]];
          If[noProcessSymmetricRemSeq && SameQ[remSeq, Reverse[remSeq]], 
               Continue[]; 
          ]; 
          
          rsFn = computeRemSeqFormulaFn[remSeq, cfgOptions]; 
          If[requireRemSeqFormula && !rsFn[[1]], 
               Continue[]; 
          ]; 

          PrintMatchInfo["Rem Seq -> ", 
                         rsFn[[2]] @@ Map[SymbolName[#1]&, {j,i}]]; 

          (** : Process the upper and lower index formulas for each : **) 
          (** : of the sequence factors                             : **)
          (*coeffData = matchResultData[[1]];*)
          coeffData = ExtractTSeqULIndexData[matchResultData]; 
          uiFnsList = {}; 
          liFnsList = {}; 
          breakOnFormulaError = False; 
          For[seqFactorIndex = 1, seqFactorIndex <= numSeqFactors, 
              seqFactorIndex++, 

               uiIndexOffset = indexOffsets[[seqFactorIndex]][[1]]; 
               upperIndexSeq = getUIIndexSeqFn[coeffData, 1, seqFactorIndex];
               uiFn = computeIndexFormulaFn[upperIndexSeq, uiIndexOffset, 
                                            cfgOptions]; 
               
               If[requireULIndexFormula && !uiFn[[1]], 
                    breakOnFormulaError = True; 
                    Break[]; 
               ]; 
               
               liIndexOffset = indexOffsets[[seqFactorIndex]][[2]]; 
               lowerIndexSeq = getUIIndexSeqFn[coeffData, 2, seqFactorIndex];
               liFn = computeIndexFormulaFn[lowerIndexSeq, liIndexOffset, 
                                            cfgOptions]; 

               If[requireULIndexFormula && !liFn[[1]], 
                    breakOnFormulaError = True; 
                    Break[]; 
               ]; 

               uiFnsList = Append[uiFnsList, uiFn[[4]] ]; 
               liFnsList = Append[liFnsList, liFn[[4]] ]; 
          
               PrintDebug["coeffData: ", coeffData];
               PrintDebug["upper seq (index=", seqFactorIndex, "): ", 
                          upperIndexSeq];
               PrintDebug["lower seq (index=", seqFactorIndex, "): ", 
                          lowerIndexSeq];
               PrintDebug["remSeq: ", remSeq];

               PrintMatchInfo["Upper Index (", seqFactorIndex, ") -> ", 
                              uiFn[[4]] @@ Map[SymbolName[#1]&, {j,i}]]; 
               PrintMatchInfo["Lower Index (", seqFactorIndex, ") -> ", 
                              liFn[[4]] @@ Map[SymbolName[#1]&, {j,i}]]; 

          ]; (* For seqFactorIndex *) 
          
          If[breakOnFormulaError, 
               Continue[]; 
          ]; 
          
          (*formulaData = {{{uiFn[[4]], liFn[[4]]}}, rsFn[[2]], rsFn[[3]]};*)
          ulIndexFnsData = PackageMultipleTSeqData[uiFnsList, liFnsList, Null]; 
          ulIndexFnsData = ExtractTSeqULIndexData[ulIndexFnsData]; 
          formulaData = {ulIndexFnsData, rsFn[[2]], rsFn[[3]]};
          formulaMatches = Append[formulaMatches, formulaData]; 

     ]; (* For mresult *) 

     PrintDebug["formulaMatches: ", formulaMatches]; 
     Return[formulaMatches];

]

(**** : Begin GuessPolynomialSequence[...] helper functions : ****)
     
getSequenceFactorFnGenerators[seqFactorsList_] := 
Module[{seqFactorFns, findex, seqFactorID, 
        seqFactorFn, seqMetaData, seqGenFn}, 

     seqFactorFns = {}; 
     For[findex = 1, findex <= Length[seqFactorsList], findex++, 

          seqFactorID = seqFactorsList[[findex]]; 
          seqMetaData = QuerySequenceMetaData[seqFactorID]; 
          seqGenFn = (GeneratorFunction) /. seqMetaData; 
          PrintDebug["seqGenFn: ", seqGenFn]; 
          seqFactorFns = Append[seqFactorFns, seqGenFn]; 

     ]; (* For findex *) 

     Return[seqFactorFns]; 

] 

getSequenceFactorFormulaTerm[matchData_, seqFactorFnList_] := 
Module[{seqFactorCfList, sfindex, seqFactorFn, uiInput, liInput, 
        seqFactorCf, seqFactorCfFn}, 

     seqFactorCfList = {}; 
     For[sfindex = 1, sfindex <= Length[seqFactorFnList], sfindex++, 

          seqFactorFn = seqFactorFnList[[sfindex]]; 
          uiInput = matchData[[1]][[sfindex]][[1]]; 
          liInput = matchData[[1]][[sfindex]][[2]]; 
          seqFactorCf = seqFactorFn[uiInput[#1, #2], liInput[#1, #2]]; 
          seqFactorCfList = Append[seqFactorCfList, seqFactorCf]; 

     ]; 

     seqFactorCfFn = Times @@ seqFactorCfList; 
     Return[seqFactorCfFn]; 

]

(**** : End GuessPolynomialSequence[...] helper functions : ****)

Options[GuessPolynomialSequence] = DefaultPkgConfigOptions; 
GuessPolynomialSequence[polyseq_, polyVar_, 
                        cfgOptions : OptionsPattern[]] := 
Module[{userGuessFn, seqStartIndex, seqFactorsList, factorFn, 
        procPolySeq, xvarPowMins, polyMaxDegrees, seqLength, 
        pindex, polyIndex, initPolyForm, procPolyData, newPoly, 
        minPolyDegree, maxPolyDegree, polyVarFactorPow, 
        testIndexOffsetParams, indexMultiples, offsetValues, offsetPairs, 
        matches, oindex, indexOffsets, processFnArgs, processFn, 
        mapProcessFn, fullMatchData, remSeqPreprocessFns, ppi, 
        ppFn, md, md2Matches, polyDegreeSeq, degreeFormula, 
        seqFactorFns, matchingFormulaFns, mindex, 
        seqFactorFn, uiInput, liInput, remSeqFormula, 
        seqDegreeFormula, pvar, powOffset, displayVars, seqFactorCf, 
        remSeqFn, polyPowTerm, userGuessTerm, innerSumTerms, sumBounds, 
        sumFormula, printFormulasQ, headerTextStr, printHeaderStr, 
        polyFormulaTerms, pLHS, pEQ, pRHS, polyFormulaStyle, printBullet, 
        remSeqDataList, computeSeqDiffs, seqDiffs, seqDiffsCorrect, 
        numPolySeqElements, numSeqFactorRows, clearLocalSeqData, 
        formulaCountMax, displayVarNames}, 

     (** : Pre-processing options and check other runtime options: **) 
     If[OptionValue[EnableDebugging], 
          EnableDebugging[], 
          DisableDebugging[]; 
     ]; 
     
     GuessPolySequenceFormulas`PkgConfig`AllowSymbolicData = 
          OptionValue[AllowSymbolicData]; 

     (** : Get configuration option settings: **) 
     userGuessFn = OptionValue[UserGuessFunction];
     seqStartIndex = OptionValue[StartIndex];
     seqFactorsList = OptionValue[SequenceFactors]; 
     factorFn = OptionValue[FactorFunction];
     
     (** : Setup the local sequence handling if the user did not : **) 
     (** : specifiy an alternate handling function:                **)
     If[factorFn == Null, 
          If[!BuildLocalSequenceData[cfgOptions], 
               Return[{}]; 
          ]; 
          (*factorFn = FactorIntegerBySequence; *)
          factorFn = FactorIntegerBySequences; 
     ]; 
     
     PrintStatus["In this function ... ", "StartIndex->", seqStartIndex]; 
     
     (** : Pre-process the input polynomial sequence: **) 
     procPolySeq = {}; 
     xvarPowMins = {}; 
     polyMaxDegrees = {};
     seqLength = Length[polyseq];
     For[pindex = 1, pindex <= Length[polyseq], pindex++, 
          
          polyIndex = seqStartIndex + pindex - 1;
          initPolyForm = polyseq[[pindex]];
          procPolyData = PreProcessPolynomialData[initPolyForm, polyIndex, 
                                                  polyVar, userGuessFn, 
                                                  cfgOptions]; 
          If[Length[procPolyData] == 0, 
               Return[{}]; 
          ]; 
          
          newPoly = (ModPolyFn /. procPolyData); 
          procPolySeq = Append[procPolySeq, newPoly]; 
          minPolyDegree = Exponent[newPoly, polyVar, Min]; 
          maxPolyDegree = Exponent[newPoly, polyVar]; 
          xvarPowMins = Append[xvarPowMins, minPolyDegree]; 
          polyMaxDegrees = Append[polyMaxDegrees, maxPolyDegree]; 

     ]; 

     PrintDebug["procPolySeq: ", procPolySeq]; 
     PrintDebug["min exponents list: ", xvarPowMins]; 

     (** : Figure out if can factor out a multiple of the variable x^pow: **) 
     polyVarFactorPow = Min[xvarPowMins];
     If[polyVarFactorPow > 0, 
          procPolySeq = Cancel[procPolySeq / Power[polyVar, polyVarFactorPow]];
          polyMaxDegrees = polyMaxDegrees - polyVarFactorPow; 
          PrintStatus[StringForm["Canceled factor of ``^", polyVar], 
                      polyVarFactorPow, ": " procPolySeq]; 
          PrintStatus["Call Message for this status ... "]; 
     ]; 

     (** : Get the sets of index offset parameters: **) 
     (*testIndexOffsetParams = { {{0, 1}} };*)            (* initial setting *)
     (*testIndexOffsetParams = { {{0, 1}}, {{0, -1}} };*) (* initial setting *)
     (*testIndexOffsetParams = { {{0, -1}, {0, 1}} };*)   (* initial setting *)
     testIndexOffsetParams = OptionValue[IndexOffsetPairs]; 
     If[testIndexOffsetParams == Null, 
          indexMultiples = OptionValue[IndexMultiples]; 
          offsetValues = Union[Flatten[Map[{#1, -#1}&, indexMultiples]]]; 
          offsetPairs = Tuples[offsetValues, 2]; 
          testIndexOffsetParams = Tuples[offsetPairs, Length[seqFactorsList]]; 
     ]; (* otherwise, use the user-defined setting *) 

     (** : Loop over the possible offset and input options : **) 
     matches = {}; 
     For[oindex = 1, oindex <= Length[testIndexOffsetParams], oindex++, 
          
          indexOffsets = testIndexOffsetParams[[oindex]]; 
          PrintStatus["indexOffsets: ", indexOffsets]; 

          processFnArgs = {#1, seqStartIndex + First[#2] - 1, 
                           polyVar, indexOffsets, factorFn, 
                           cfgOptions}&; 

          processFn = ProcessSingleFactorPoly[##]&;
          mapProcessFn = (processFn @@ processFnArgs[#1, #2])&; 
          fullMatchData = MapIndexed[mapProcessFn[#1, #2]&, procPolySeq]; 
          PrintMatchInfo["fullMatchData: ", fullMatchData];
          
          If[IsEmptySet[fullMatchData], 
               Continue[];
          ];
          
          (** : Process this information: **) 
          remSeqPreprocessFns = {Identity, Reverse}; 
          For[ppi = 1, ppi <= Length[remSeqPreprocessFns], ppi++, 

               ppFn = remSeqPreprocessFns[[ppi]];
               PrintStatus[StringForm["Pre matches: [oi=``, ppi=``]", 
                           oindex, ppi], matches];
               md = FilterByRemSequenceMatches[fullMatchData, ppFn, 
                                               cfgOptions];
               PrintStatus["Post matches: ", matches];
               PrintMatchInfo["After FilterByRemSequenceMatches (md): ", md];
               
               If[ppFn === Identity, 
                    SetOptions[VerifyFormulaMatches, 
                               ReverseRemSeqFormulaIndex -> False];
                    SetOptions[VerifyFormulaMatches, 
                               AbortOnSymmetricRemSequence -> False];
               ]; 
               If[ppFn === Reverse, 
                    SetOptions[VerifyFormulaMatches, 
                               ReverseRemSeqFormulaIndex -> True];
                    SetOptions[VerifyFormulaMatches, 
                               AbortOnSymmetricRemSequence -> True]; 
               ]; 
              
               md2Matches = VerifyFormulaMatches[md, indexOffsets, 
                                                 cfgOptions];
               PrintStatus["matches: ", matches];
               PrintStatus["md2Matches: ", md2Matches]; 
               matches = Join[matches, md2Matches]; 
               PrintMatchInfo["After FilterByRemSequenceMatches: ", 
                              Short[matches, 2]];

          ];

     ]; 
     PrintDebug["!!!!!!!!!!! [After main loop] ==================="];

     (** : Get the sequence of polynomial degrees: **) 
     polyDegreeSeq = Map[Exponent[#1, polyVar]&, procPolySeq]; 
     degreeFormula = ComputeSequenceFormula[polyDegreeSeq, 
                                            seqStartIndex, cfgOptions]; 
     
     If[!degreeFormula[[1]], 
          Message[GPSFPkgMsgs::GuessMultipleFactorPolySequence::PolySeqDegree]; 
          If[Length[polyseq] < 6, 
               Message[GPSFPkgMsgs::Warnings::InsufficientSeqElements, 
                       Length[polyseq], 8]; 
          ]; 
          Return[{}];
     ]; 
     PrintStatus["sequence degreeFormula: ", 
                 degreeFormula[[2]] @@ Map[SymbolName[#1]&, {j,i}]]; 

     (** : Format and print the computed formulas:                      **) 
     (** : Later, create separate formula creation routines (options) : **)
     (** : and printing functions                                     : **)
     Off[Sum::itraw]; (* temporarily disable this message for 
                         formula creation and printing *)

     (** : Get list of sequence factor function generators: **) 
     (*seqFactorFns = Map[QuerySequenceMetaDataByProperty[#1, 
                          GeneratorFunction]&, seqFactorsList]; *) 
     seqFactorFns = getSequenceFactorFnGenerators[seqFactorsList]; 
     
     matchingFormulasFns = {}; 
     formulaCountMax = Min[Length[matches], OptionValue[LimitFormulaCount]]; 
     For[mindex = 1, mindex <= formulaCountMax, mindex++, 

          (** : Get the components to assemble the : **) 
          (** : matching formula functions:        : **)
          seqFactorCf = getSequenceFactorFormulaTerm[matches[[mindex]], 
                                                     seqFactorFns]; 
          remSeqFormula = matches[[mindex]][[2]]; 
          seqDegreeFormula = degreeFormula[[2]]; 
 
          (** : Pretty printing of the formula variables without context: **)
          pvar = polyVar; 
          powOffset = polyVarFactorPow; 
          
          displayVarNames = {j, i}; 
          If[Length[OptionValue[DisplayVars]] == 2, 
               displayVarNames = OptionValue[DisplayVars]; 
          ]; 
          displayVarNames = Append[displayVarNames, pvar]; 
          (*displayVars = Map[SymbolName[#1]&, {j, i, pvar}] *)
          displayVars = Map[SymbolName[#1]&, displayVarNames]; 

          (** : Create th formula function: **) 
          (** : Input order: {j, i, x} (Poly index, Sum index, Poly Var) : **)
          remSeqFn = remSeqFormula[seqDegreeFormula[#1], #2]; 
          polyPowTerm = Power[#3, #2 + powOffset]; 
          userGuessTerm = userGuessFn[#1, #2]; 
          innerSumTerms = Times @@ {seqFactorCf, remSeqFn, userGuessTerm, 
                                    polyPowTerm}; 
          sumBounds = {#2, 0, seqDegreeFormula[#1]}; 
          (*% \label{scode_linelbl_GPSF.m_With_usage_example_v2} %*)
          sumFormula = With[{st = innerSumTerms, sb = sumBounds}, 
                            Function @@ HoldForm[Sum[st, sb]]];

          matchingFormulasFns = Append[matchingFormulasFns, sumFormula]; 

          (** : Print the formulas to the notebook (or skip option): **)
          printFormulasQ = OptionValue[PrintFormulas]; 
          If[!printFormulasQ, 
               Continue[]; 
          ]; 

          (** : Print a separate header for each match: **) 
          headerTextStr = StringForm["Found Matching Formula #`` / ``:", 
                                    mindex, Length[matches]]; 
          printHeaderStr = "============ " <> 
                           ToString[headerTextStr] <> 
                           " ============"; 
          Print[Style[printHeaderStr, White, Bold, 
                      Background -> Darker[Blue, 0.36], FontSize -> 14]]; 

          (** : Display the polynomial formula: **) 
          polyFormulaTerms = TraditionalForm[ 
               (Subscript["Poly", #1][#3]& @@ displayVars) 
               \[RightTeeArrow] (sumFormula @@ displayVars) 
          ]; 

          pLHS = TraditionalForm[Subscript["Poly", #1][#3]& @@ displayVars]; 
          pEQ = " \[RightTeeArrow] "; 
          pRHS = TraditionalForm[sumFormula @@ displayVars]; 
          polyFormulaStyle = Style[#1, Bold, FontSize -> 15]&; 
          Print[" ", polyFormulaStyle[pLHS], polyFormulaStyle[pEQ], 
                     polyFormulaStyle[pRHS]]; 

          (** : Output metadata as a bulleted list: **) 
          printBullet = Print[" ", "\[FilledRightTriangle]", " ", 
                             Style[#1, Purple, Bold, FontSize -> 12], ##2]&; 

          If[OptionValue[PrintLaTeXFormulas], 
               printBullet["Latex Formula Output: ", Null]; 
          ]; 
          
          remSeqDataList = matches[[mindex]][[3]]; 
          printBullet["Remaining Sequence Data: ", 
                      Shallow[remSeqDataList, {Infinity, 10}]]; 
 
          printBullet["User Function: ", 
               TraditionalForm[Subscript["U", "guess"][#1, #2]& 
                               @@ displayVars], " \[LongEqual] ", 
               TraditionalForm[userGuessFn[#1, #2]& @@ displayVars]]; 
          printBullet["Formula Function: ", 
               TraditionalForm[Subscript["PolyFormula", 
                               StringForm["index=``", mindex]][##]& 
                               @@ displayVars] 
          ]; 
          (*printBullet["Full Form", " \[LongRightArrow] ", 
                        FullForm[sumFormula @@ displayVars]]; *)
          
          computeSeqDiffs = OptionValue[DiffPolySequenceFormulas]; 
          If[!computeSeqDiffs, 
               Continue[]; 
          ]; 
          
          seqDiffs = DiffSequenceFormula[polyseq, polyVar, sumFormula, 
                                         cfgOptions]; 
          seqDiffsCorrect = (Count[seqDiffs, False] == 0); 
          printBullet["Sequence Formula Diffs: ", 
                      Short[seqDiffs, 0.5], 
                      Which[seqDiffsCorrect, 
                            Style[" [\[Checkmark]]", Medium, Bold, Green], 
                            !seqDiffsCorrect, 
                            Style[" [\[ScriptX]]", Medium, Bold, Red] 
                      ] 
          ]; 

     ];
     On[Sum::itraw]; (* turn this message back on *)

     (** : Issue possible warning messages if no formulas for the : **) 
     (** : sequence are found:                                      **)
     If[IsEmptySet[matchingFormulasFns] || OptionValue[IssueAllWarnings], 
          
          numPolySeqElements = Length[polyseq]; 
          numSeqFactorRows = OptionValue[TriangularSequenceNumRows]; 
          If[numPolySeqElements < 6, 
               Message[GPSFPkgMsgs::Warnings::InsufficientSeqElements, 
                       numPolySeqElements, 8]; 
          ]; 
          If[numSeqFactorRows < (numPolySeqElements + 4), 
               Message[GPSFPkgMsgs::Warnings::InsufficientFactorData, 
                       numSeqFactorRows, numPolySeqElements + 4]; 
          ]; 

     ]; (* end warning messages *) 

     (** : Clean up local sequence data storage: **)
     clearLocalSeqData = OptionValue[ClearLocalSequenceData]; 
     If[clearLocalSeqData, 
          PkgSequenceData = {}; 
          PkgNumSequenceFactors = 0; 
     ]; 

     Return[matchingFormulasFns]; 

]

End[] (* Private *)

(*************************************************************************)
(**** : Perform pretty printing of the package details and          : ****) 
(**** : revsision information if the package is loaded in a notebook: ****)
(*************************************************************************)
packageInfoStr = "" <> 
     GuessPolySequenceFormulas`PackageName <> ":\n" <> 
     GuessPolySequenceFormulas`PackageShortDesc <> "\n" <> 
     "Author: " <> GuessPolySequenceFormulas`PackageAuthor <> "\n" <> 
     "Package Version: " <> GuessPolySequenceFormulas`PackageVersion; 

If[$Notebooks,
     CellPrint[Cell[packageInfoStr, "Print", 
                    FontColor -> RGBColor[0, 0, 0], 
                    FontSize -> 12, 
                    CellFrame -> 0.5, 
                    CellFrameColor -> Purple, 
                    Background -> LightBlue, 
                    CellFrameMargins -> 14
                   ]
              ],
     Print[packageInfoStr]; 
]; 

EndPackage[]

(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
(*************************************************************************)

