(*************************************************************************)
(*************************************************************************)
(********************* : GuessSequenceData.m:         ********************) 
(*************************************************************************)
(*************************************************************************)

BeginPackage["GuessSequenceData`"]

(*************************************************************************)
(**** : Package revision, metadata information, and settings:         ****)
(*************************************************************************)

(*% \label{scode_linelbl_GSD.m_test_v2} %*)
LocalPackageName = "GuessSequenceData.m";
LocalPackageVersion = "2014.04.28-v5";
LocalPackageAuthor = "Maxie D. Schmidt";
LocalPackageShortDesc = 
     "Default package to provide several built-in sequences for the " <> 
     "GuessPolySequenceFormulas.m software package."

GSDPkg::FormattingData::NewlineString = "\n";
GSDPkg::FormattingData::ListBulletDelim = "  \[FilledRightTriangle] "; 

(*************************************************************************)
(**** : Create default usage information for the package functions:   ****)
(**** : Provide detailed usage information for the package functions: ****)
(*************************************************************************)

(**** : Utilities / sequence lookup functions in the package: ****)
LookupSequenceKey::usage = 
     "LookupSequenceKey[seqID]: (Intended for internal use by the " <> 
     "GuessPolySequenceFormulas.m and local GuessSequenceData.m " <> 
     "packages)" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Used to obtain the internal index for the lookup table that " <> 
     "stores the metadata information about the sequences " <> 
     "(as specified by the short or long seqID key input to the " <> 
     "function) supported by the GuessSequenceData subpackage." <> 
     "The sequences currently supported by thie package are listed by " 
     "running the local package function ListSupportedSequences[]."; 

QuerySequenceMetaData::usage = 
     "QuerySequenceMetaData[seqID]:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Returns the complete listing of metadata options stored " <> 
     "locally by the package for the sequence specified by the " <> 
     "(short or long) form of the seqID key " <> 
     "input to the function. A complete listing of the short and " <> 
     "long forms of the seqID keys currently supported by the " <> 
     "package can be viewed by running the local package function " <> 
     "ListSupportedSequences[]." <> 
     GSDPkg::FormattingData::NewlineString <> 
     "The sequence information returned is a listing of sequence " 
     "transformation rules corresponding to the public options " <> 
     "documented in the package source below." <> 
     "The listing of these options may also be viewed by " <> 
     "querying the usage from the Mathematica help lookup for the " <> 
     "function: ?QuerySequenceMetaDataByProperty" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Example Usage:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "QuerySequenceMetaData[\"S1\"]" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "seqFn = (GeneratorFunction) /. QuerySequenceMetaData[\"S2\"];" <> 
     "Table[seqFn[n, k], {n,0,10}, {k,0,10}] // TableForm" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "{Desc, NumInputs, OEISRefs} /. " <> 
     "QuerySequenceMetaData[\"EulerianE1\"]"; 

QuerySequenceMetaDataByProperty::usage = 
     "QuerySequenceMetaDataByProperty[seqID, mdataProp]:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "The function returns the metadata associated with the sequence " <> 
     "key seqID for the property mdataProp only. " <> 
     "Sequence Properties include: " <> 
     "LongSeqID, ShortSeqID, Desc, NumInputs, " <> 
     "GeneratorFunction, RowRangeFunction, " <> 
     "FormulaDisplay, LatexDisplay, Parameters, " <> 
     "Comments, OEISRefs, PrintRefs, WebRefs, OtherRefs, " <> 
     "RevisionInfo, Implementation, ImplNotes, OtherMetaData." <> 
     "The related package function QuerySequenceMetaData[seqID] " <> 
     "returned a complete listing of all metadata options for the " <> 
     "sequence in the form of list of replacement rules." <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Example Usage:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "descStr = QuerySequenceMetaDataByProperty[\"StirlingS2\",Desc];" <> 
     "Print[\"Description of Sequence S2: \", descStr];" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "seqFn = QuerySequenceMetaDataByProperty[\"S1\", GeneratorFunction];" <> 
     "Table[Power[-1,n-k]*seqFn[n, k] - StirlingS1[n,k], " <> 
     "{n,0,8}, {k,0,8}] // TableForm"; 

ListSequenceMetaData::usage = 
     "ListSequenceMetaData[seqID, fullInfoQ:False, numPrintRows:6]:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Provides a short summary of the sequence metadata information " <> 
     "when fullInfoQ is False, or a complete listing of the local " <> 
     "sequence information when the second parameter to the " 
     "function is set to True. The third paremeter numPrintRows " <> 
     "specifies how many rows of the triangular sequence values to " <> 
     "print in the display returned by the function." 
     GSDPkg::FormattingData::NewlineString <>    
     "Example Usage:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "ListSequenceMetaData[\"BinomialSquared\"];" <> 
     GSDPkg::FormattingData::NewlineString <> 
     GSDPkg::FormattingData::ListBulletDelim <> 
     "ListSequenceMetaData[\"E2\",False,10];"; 

ListSupportedSequences::usage = 
     "ListSupportedSequences[]:" <> 
     GSDPkg::FormattingData::NewlineString <> 
     "Displays a summary of all of the triangular sequences " <> 
     "supported and/or implemented by the " <> 
     "GuessSequenceData subpackage."; 

(**** : Declare public options / replacement rules for the : ****)
(**** : built-in sequence metadata stored in the package:    ****)
LongSeqID::usage         = "GuessSequenceData subpackage metadata option.";
ShortSeqId::usage        = "GuessSequenceData subpackage metadata option.";
Desc::usage              = "GuessSequenceData subpackage metadata option.";
NumInputs::usage         = "GuessSequenceData subpackage metadata option.";
StartRowIndex::usage     = "GuessSequenceData subpackage metadata option.";
AllowZeros::usage        = "GuessSequenceData subpackage metadata option.";
GeneratorFunction::usage = "GuessSequenceData subpackage metadata option.";
RowRangeFunction::usage  = "GuessSequenceData subpackage metadata option.";
FormulaDisplay::usage    = "GuessSequenceData subpackage metadata option.";
LatexDisplay::usage      = "GuessSequenceData subpackage metadata option.";
Parameters::usage        = "GuessSequenceData subpackage metadata option.";
Comments::usage          = "GuessSequenceData subpackage metadata option.";
OEISRefs::usage          = "GuessSequenceData subpackage metadata option.";
PrintRefs::usage         = "GuessSequenceData subpackage metadata option.";
WebRefs::usage           = "GuessSequenceData subpackage metadata option.";
OtherRefs::usage         = "GuessSequenceData subpackage metadata option.";
RevisionInfo::usage      = "GuessSequenceData subpackage metadata option.";
Implementation::usage    = "GuessSequenceData subpackage metadata option.";
ImplNotes::usage         = "GuessSequenceData subpackage metadata option.";
OtherMetaData::usage     = "GuessSequenceData subpackage metadata option.";

(*************************************************************************)
(**** : Package error handling and messages:                          ****)
(*************************************************************************)
GuessSequenceData::ErrorMsgs::InvalidSequenceID = 
     "Invalid sequence ID `1`"; 
GuessSequenceData::ErrorMsgs::InvalidMetadataOption = 
     "Invalid metadata property `1`"; 

GuessSequenceData::WarningMsgs::MultipleSeqKeyMatches = 
     "Multiple sequence match ID `1` ... using first match key `2`"; 

Begin["`Private`"] 

(**** : Built-in standard functions supported by the (sub)package:      ****)
(**** : Index Format: {LongIDKey, ShortIDKey, Lookup table ID, ValidQ}: ****) 
PkgSupportedSequenceIDKeys = { 
     {"StirlingS1",         "S1",         "SeqIDKey:S1:",        True}, 
     {"Stirlings1",         "s1",         "SeqIDKey:s1:",        True}, 
     {"StirlingS2",         "S2",         "SeqIDKey:S2:",        True}, 
     {"EulerianE1",         "E1",         "SeqIDKey:E1:",        True}, 
     {"EulerianE2",         "E2",         "SeqIDKey:E2:",        True}, 
     {"Binomial",           "Binom",      "SeqIDKey:Binom:",     True}, 
     {"BinomialSquared",    "Binom2",     "SeqIDKey:BinomPow2:", True}, 
     {"BinomialSymmetric",  "BinomSym",   "SeqIDKey:BinomSym:",  True}, 
     {"Documentation",      "ExampleFn",  "SeqIDKey:Example:",   False} 
}; 

LookupSequenceKey[seqID_] := 
Module[{posValues, seqInfoIndices}, 
     
     posValues = Position[PkgSupportedSequenceIDKeys, seqID];
     seqInfoIndices = Map[First[#1]&, posValues]; 
     Return[ Map[PkgSupportedSequenceIDKeys[[#1]][[3]]&, seqInfoIndices] ];

] 

QuerySequenceMetaData[seqID_] := 
Module[{}, (*{seqLookupKeys, seqLookupKey, mdataIndexPos, mdataIndex}, *)

     seqLookupKeys = LookupSequenceKey[seqID]; 
     If[Length[seqLookupKeys] == 0, 
          Message[GuessSequenceData::ErrorMsgs::InvalidSequenceID, seqID]; 
          Return[{}]; 
     ]; 
     
     If[Length[seqLookupKeys] > 1, 
          Message[GuessSequenceData::WarningMsgs::MultipleSeqKeyMatches, 
                  seqID, seqLookupKeys[[1]]];
     ]; 
     seqLookupKey = seqLookupKeys[[1]]; 
     
     mdataIndexPos = Position[FullSequencesMetaData, seqLookupKey, 2, 1]; 
     If[Length[mdataIndexPos] == 0, 
          Return[{}];
     ]; 
     
     (** : Format: { {index, 1} }: **) 
     mdataIndex = mdataIndexPos[[1]][[1]]; 
     Return[Drop[FullSequencesMetaData[[mdataIndex]], 1]]; 

]

QuerySequenceMetaDataByProperty[seqID_, mdataProp_] := 
Module[{}, (*{seqMetaData, rProp}, *)

     seqMetaData = QuerySequenceMetaData[seqLookupKey]; 
     If[Length[seqMetaData] == 0, 
          Return[Null]; 
     ]; 
     
     rProp = (mdataProp) /. seqMetaData; 
     If[rProp == mdataProp, 
          Message[GuessSequenceData::ErrorMsgs::InvalidMetadataOption, 
                  mdataProp]; 
          Return[Null]; 
     ]; 
     Return[rProp]; 

]

ListSequenceMetaData[seqID_, fullInfoQ_:False, 
                     numPrintRows_:6] := 
Module[{}, 
 
     seqProps = {LongSeqID, ShortSeqID, Desc, NumInputs}; 
     If[fullInfoQ, 
          remProps = {FormulaDisplay, LatexDisplay, Parameters, 
                      Comments, OEISRefs, PrintRefs, WebRefs, OtherRefs,    
                      RevisionInfo, Implementation, ImplNotes, 
                      OtherMetaData}; 
          seqProps = Append[seqProps, remProps]; 
     ]; 
     
     seqGenFn = QuerySequenceMetaDataByProperty[seqID, GeneratorFunction]; 
     rowRangeFn = QuerySequenceMetaDataByProperty[seqID, RowRangeFunction]; 
     
     Print["ListSequenceMetaData: Implementation later "]; 
     Print["Get list of sequence values ... "]; 
] 

ListSupportedSequences[overrideIsValidTag_:False] := 
Module[{}, 

     pkgSeqKeyData = PkgSupportedSequenceIDKeys; 
     Print["ListSupportedSequences[]: Implement later ... "]; 
     Print["Return a list of {LongKey, ShortKey} pairs ... "]; 

] 

(** : Previous, older versions of the Striling number triangles: **)
SequenceGeneratorS1[n_, k_] := Power[-1, n-k] * StirlingS1[n, k]; 
SequenceGeneratorS1[n_, k_] := Abs[StirlingS1[n, k]]; 
SequenceGeneratorS2[n_, k_] := StirlingS2[n, k]; 

SequenceGeneratorS1[n_, k_] := SeqFnS1[n, k]; 
SequenceGenerators1[n_, k_] := StirlingS1[n, k]
SequenceGeneratorS2[n_, k_] := SeqFnS1[n, k]; 
SequenceGeneratorE1[n_, k_] := SeqFnE1[n, k]; 
SequenceGeneratorE2[n_, k_] := SeqFnE2[n, k];  
SequenceGeneratorBinom[n_, k_] := Binomial[n, k]; 
SequenceGeneratorBinomSquared[n_, k_] := Power[Binomial[n, k], 2]; 
SequenceGeneratorBinomSymmetric[n_, m_] := Binomial[n + m, m]; 
SequenceGeneratorDefaultFn[n_, k_] := 1; 

SequenceGeneratorS2[n_, k_] := StirlingS2[n, k]; 

SequenceRowRangeDefaultTSeqV1[n_Integer] := 
     Which[n < 0, {0, -1}, n == 0, {0, 0}, n >= 1, {1, n}] 
SequenceRowRangeDefaultTSeqV2[n_Integer] := 
     Which[n < 0, {0, -1}, n == 0, {0, 0}, n >= 1, {0, n - 1}] 
SequenceRowRangeDefaultTSeqV3[n_Integer] := 
     Which[n < 0, {0, -1}, n == 0, {0, 0}, n >= 1, {0, n}] 
SequenceRowRangeDefaultFn[n_Integer] := If[n < 1, {0, -1}, {1, n}]

SequenceRowRangeS1[n_] := SequenceRowRangeDefaultTSeqV1[n] 
SequenceRowRanges1[n_] := SequenceRowRangeDefaultTSeqV1[n] 
SequenceRowRangeS2[n_] := SequenceRowRangeDefaultTSeqV1[n] 
SequenceRowRangeE1[n_] := SequenceRowRangeDefaultTSeqV2[n] 
SequenceRowRangeE2[n_] := SequenceRowRangeDefaultTSeqV2[n] 
SequenceRowRangeBinom[n_] := SequenceRowRangeDefaultTSeqV3[n] 
SequenceRowRangeBinomSquared[n_] := SequenceRowRangeDefaultTSeqV3[n] 
SequenceRowRangeBinomSymmetric[n_] := SequenceRowRangeDefaultTSeqV3[n] 

FullSequencesMetaData = { 

  {"SeqIDKey:S1:", 
    LongSeqID         -> "StirlingS1", 
    ShortSeqId        -> "S1", 
    Desc              -> "Stirling numbers of the first kind " <> 
                         "(unsigned triangle)", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorS1, 
    RowRangeFunction  -> SequenceRowRangeS1, 
    FormulaDisplay    -> "S1(``, ``)", 
    LatexDisplay      -> "\\genfrac{``}{``}", 
    Parameters        -> Null, 
    Comments          -> {"Sequence has OGF ... "}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {"See \\S 6.1 of the Concrete Mathematics " <> 
                          "book and \\S 26.8 in the NIST Handbook."}, 
    WebRefs           -> {"MathWorld: StirlingNumber.html"}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:s1:", 
    LongSeqID         -> "Stirlings1", 
    ShortSeqId        -> "s1", 
    Desc              -> "Stirling numbers of the first kind " <> 
                         "(signed triangle)", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGenerators1, 
    RowRangeFunction  -> SequenceRowRanges1, 
    FormulaDisplay    -> "s(``, ``)", 
    LatexDisplay      -> "\\genfrac{``}{``}", 
    Parameters        -> Null, 
    Comments          -> {"Sequence has OGF ... "}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "Standard Mathematica function (StrlingS1)", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:S2:", 
    LongSeqID         -> "StirlingS2", 
    ShortSeqId        -> "S2", 
    Desc              -> "Stirling numbers of the second kind", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorS2, 
    RowRangeFunction  -> SequenceRowRangeS2, 
    FormulaDisplay    -> "S(``, ``)", 
    LatexDisplay      -> "\\genfrac{``}{``}", 
    Parameters        -> Null, 
    Comments          -> {"Sequence has OGF ... "}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {"See \\S 6.1 of the Concrete Mathematics " <> 
                          "book and \\S 26.8 in the NIST Handbook."}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "Standard Mathematica function (StrlingS2)", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:E1:", 
    LongSeqID         -> "EulerianE1", 
    ShortSeqId        -> "E1", 
    Desc              -> "Triangle of the first-order Eulerian numbers", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorE1, 
    RowRangeFunction  -> SequenceRowRangeE1, 
    FormulaDisplay    -> "", 
    LatexDisplay      -> "", 
    Parameters        -> Null, 
    Comments          -> {""}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {"See \\S 6.2 of the Concrete Mathematics " <> 
                          "book and \\S 26.14 in the NIST Handbook."}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "", 
    ImplNotes         -> "See also EulerianE1 in the " <> 
                         "RISC Stirling.m package", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:E2:", 
    LongSeqID         -> "EulerianE2", 
    ShortSeqId        -> "EulerianE2", 
    Desc              -> "Triangle of the second-order Eulerian numbers", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorE2, 
    RowRangeFunction  -> SequenceRowRangeE2, 
    FormulaDisplay    -> "", 
    LatexDisplay      -> "", 
    Parameters        -> Null, 
    Comments          -> {""}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {"See \\S 6.2 of the Concrete Mathematics book."}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "", 
    ImplNotes         -> "See also EulerianE2 in the " <> 
                         "RISC Stirling.m package", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:Binom:", 
    LongSeqID         -> "Binomial", 
    ShortSeqId        -> "Binom", 
    Desc              -> "Binomial coefficients", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorBinom, 
    RowRangeFunction  -> SequenceRowRangeBinom, 
    FormulaDisplay    -> "Binom(``, ``)", 
    LatexDisplay      -> "\binom{``}{``}", 
    Parameters        -> Null, 
    Comments          -> {}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "Standard Mathematica function (Binomial)", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:BinomPow2:", 
    LongSeqID         -> "BinomialSquared", 
    ShortSeqId        -> "Binom2", 
    Desc              -> "Squares of the binomial coefficients", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorBinomSquared, 
    RowRangeFunction  -> SequenceRowRangeBinomSquared, 
    FormulaDisplay    -> "", 
    LatexDisplay      -> "", 
    Parameters        -> Null, 
    Comments          -> {}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "Standard Mathematica function (Binomial)", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:BinomSym:", 
    LongSeqID         -> "BinomialSymmetric", 
    ShortSeqId        -> "BinomSym", 
    Desc              -> "Binomial coefficients defined over " <> 
                         "symmetric index inputs. " <> 
                         "(see known generating functions)", 
    NumInputs         -> 2, 
    StartRowIndex     -> 0, 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorBinomSymmetric, 
    RowRangeFunction  -> SequenceRowRangeBinomSymmetric, 
    FormulaDisplay    -> "B(`1`+`2`, ``1)", 
    LatexDisplay      -> "\binom{`1`+`2`}{`1`}", 
    Parameters        -> Null, 
    Comments          -> {}, 
    OEISRefs          -> {}, 
    PrintRefs         -> {}, 
    WebRefs           -> {}, 
    OtherRefs         -> {}, 
    RevisionInfo      -> "", 
    Implementation    -> "Standard Mathematica function (Binomial)", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null 
  }, 

  {"SeqIDKey:Example:", 
    LongSeqID         -> "Long identifier string for the sequence", 
    ShortSeqId        -> "Short ID string", 
    Desc              -> "Description of the sequence", 
    NumInputs         -> 2, (* number of input variables *)
    StartRowIndex     -> 0, (* starting row for triangular sequences *) 
    AllowZeros        -> False, 
    GeneratorFunction -> SequenceGeneratorDefaultFn, 
    RowRangeFunction  -> SequenceRowRangeDefaultFn, 
    FormulaDisplay    -> "Display string for printing " <> 
                         "local Mathematica notebook formulas", 
    LatexDisplay      -> "Display string for LaTeX outputs", 
    Parameters        -> Null, (* parameters for the sequence *) 
    Comments          -> {}, 
    OEISRefs          -> {"List of related OEIS sequence numbers"}, 
    PrintRefs         -> {"List of relevant print references"}, 
    WebRefs           -> {"List of relevant website references"}, 
    OtherRefs         -> {"Other relevant references for the sequence"}, 
    RevisionInfo      -> "", 
    Implementation    -> "", 
    ImplNotes         -> "", 
    OtherMetaData     -> Null (* possibly for future use *)
  } 

}; 

(*************************************************************************)
(**** : Implementation of some non-standard triangular sequence     : ****)
(**** : functions defined by triangular recurrence relations:         ****)
(*************************************************************************)

Unprotect[SeqFnS1, SeqFnS2, SeqFnE1, SeqFnE2, GenTRec]; 

GenTRec[n_Integer, k_Integer, 
        alpha_, beta_, gamma_, alpha2_, beta2_, gamma2_] := 
        GenTRecLocal[n, k, alpha, beta, gamma, alpha2, beta2, gamma2]; 
     
GenTRecLocal[n_Integer, k_Integer, a_, b_, g_, a2_, b2_, g2_] := 
Which[ 
     n < 0 || k < 0, 0, 
     n == 0, KroneckerDelta[n, k], 
     n >= 0, GenTRecLocal[n, k, a, b, g, a2, b2, g2] = 
             (a n + b k + g) *  
             GenTRecLocal[n - 1, k, a, b, g, a2, b2, g2] + 
             (a2 n + b2 k + g2) * 
             GenTRecLocal[n - 1, k - 1, a, b, g, a2, b2, g2]
]; 

SeqFnS1[n_Integer, k_Integer] := GenTRec[n, k, 1, 0, -1, 0, 0, 1]; 
Format[SeqFnS1[n_, k_], TraditionalForm] := Subscript["S", 1][n, k]; 

SeqFnS2[n_Integer, k_Integer] := GenTRec[n, k, 0, 1, 0, 0, 0, 1]; 
Format[SeqFnS2[n_, k_], TraditionalForm] := Subscript["S", 2][n, k]; 

SeqFnE1[n_Integer, k_Integer] := GenTRec[n, k, 0, 1, 1, 1, -1, 0]; 
Format[SeqFnE1[n_, k_], TraditionalForm] := Subscript["E", 1][n, k]; 

SeqFnE2[n_Integer, k_Integer] := GenTRec[n, k, 0, 1, 1, 2, -1, -1]; 
Format[SeqFnE2[n_, k_], TraditionalForm] := Subscript["E", 2][n, k]; 

Protect[SeqFnS1, SeqFnS2, SeqFnE1, SeqFnE2, GenTRec]; 

End[]  (* Private *)


GenTRec::usage = 
     "Paremetrized triangular sequence definition used " <> 
     "to define several special case triangles, including the " <> 
     "Stirling numbers of the first and second kinds, and the " <> 
     "Eulerian numbers of both orders."; 

SeqFnS1::usage = 
     "Computes the unsigned Stirling numbers of the first kind " <> 
     "for n,k \[GreaterEqual] 0.\n" <> 
      "See the key \"S1\" from the GuessSequenceData package."; 
SeqFnS2::usage = 
     "Computes the Stirling numbers of the second kind " <> 
     "for n,k \[GreaterEqual] 0.\n" <> 
      "See the key \"S2\" from the GuessSequenceData package."; 
SeqFnE1::usage = 
     "Computes the first-order Eulerian numbers " <> 
     "for n,k \[GreaterEqual] 0.\n" <> 
     "See the key \"E1\" from the GuessSequenceData package."; 
SeqFnE2::usage = 
     "Computes the second-order Eulerian numbers " <> 
     "for n,k \[GreaterEqual] 0.\n" <> 
     "See the key \"E2\" from the GuessSequenceData package."; 

EndPackage[]

(*************************************************************************)
(*************************************************************************)
(*************************************************************************)
(*************************************************************************)

