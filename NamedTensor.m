(* ::Package:: *)

BeginPackage["NamedTensor`"];


Unprotect@@Names["NamedTensor`*"];
ClearAll@@Names["NamedTensor`*"];


NamedTensor::usage="NamedTensor[rowNames,colNames,data] is the general data structure for a named tensor. rowNames and colNames are associations that map the names of the indices (strings) to the indices in the tensor data.

NamedTensor[indexNames,data] is a construction shorthand. indexNames is a list that in this order defines the names of the indices, which are identical for rows and columns.";


NamedMatrix::usage="NamedMatrix[data] is a shorthand that converts a rank-2 tensor into a named tensor with empty names, one row and one column index.";
NamedColVector::usage="NamedColVector[data] is a shorthand that converts a rank-1 tensor into a named tensor with one row index of empty name.";
NamedRowVector::usage="NamedRowVector[data] is a shorthand that converts a rank-1 tensor into a named tensor with one column index of empty name.";


RenameIndices::usage="RenameIndices[tensor,replacements] renames the row and column indices in a NamedTensor object according to the rules or association specified (no patterns!)";


NamedTensorAsMatrix::usage="NamedTensorAsMatrix[tensor,order] converts a NamedTensor to a rank-2 tensor (standard List or SparseArray). By default, the levels are flattened together in the order as they appear in their associations, but a custom order can be specified.";
NamedTensorMatrixFunction::usage="NamedTensorMatrixFunction[f,tensor,order] applies a function f to the NamedTensor tensor. Since f expects a standard matrix, it is first converted into a matrix (the order may be specified, else the default order is used). After this, the matrix is converted back to the named tensor. f must not change the dimensions of the matrix.";


Transpose::usage=Transpose::usage<>"

Transpose[tensor,indices] transposes only a subset of indices.";
ConjugateTranspose::usage=ConjugateTranspose::usage<>"

ConjugateTranspose[tensor,indices] transposes only a subset of indices.";


Tr::usage=Tr::usage<>"

Tr[tensor,indices] contracts a named tensor in the given indices. By default all possible indices are contracted.";


Dot::usage=Dot::usage<>"

Dot[tensors,...] contracts named tensors by automatically contracting matching names. Before a row (column) name can be re-used as a row (column) name, it must appear as a column (row) name to contract.";


TensorProduct::usage=TensorProduct::usage<>"

TensorProduct[elements] gives a tensor product of NamedTensors. Elements must be either an association or a sequence of rules. New index names are the name of the tensor followed a dot and the old name of the index, joined with a dot. If the elements are a pure sequence of NamedTensors instead, no name rewriting is perfomed, i.e. the names of the rows and columns must be unique among the individual tensors.";


TensorContract::usage=TensorContract::usage<>"

TensorContract[elements,contractions] contracts named tensors. Expects either an association, a list of rules or a sequence of rules as first argument that associates with every name a NamedTensor object. The second parameter is a list of contractions, which may even span more than one index.";


NamedOperation::usage="NamedOperation[gates,state,contraction,newNames] performs a quantum operation acting on a state. The gates will automatically also be applied in a conjugate-transpose manner on the right-hand side of the states. Refer to TensorContract for the other parameters.";


NamedRowsQ::usage="NamedRowsQ[tensor,names] gives True iff tensor contains all rows in names.

NamedRowsQ[names] is the corresponding operator form.";
NamedColumnsQ::usage="NamedColsQ[tensor,names] gives True iff tensor contains all columns in names.

NamedColumnsQ[names] is the corresponding operator form.";
NamedIndicesQ::usage="NamedIndicesQ[tensor,names] gives True iff tensor contains all rows and all columns in names.

NamedIndicesQ[names] is the corresponding operator form.

NamedIndicesQ[tensor,namesRow,namesCol] allows to specify different names for rows and columns.

NamedIndicesQ[namesRow,namesCol] is the corresponding operator form.";
NamedRows::usage="NamedRows[tensor] lists all row index names in tensor.";
NamedColumns::usage="NamedColumns[tensor] lists all column index names in tensor.";
NamedIndices::usage="NamedIndices[tensor] lists all index names in tensor.";


Begin["`Private`"];


(* pre-11.2: use own implementation for TakeList; very reduced and simplistic variant, hence in our private context *)


If[Definition[TakeList]===Null,
  TakeList[list_,counts_List]:=With[{acounts=FoldList[Plus,0,counts]},Table[list[[acounts[[i-1]]+1;;acounts[[i]]]],{i,2,Length@acounts}]];
  SetAttributes[TakeList,NHoldRest];
  Protect[TakeList];
];


(* pre-10.3: use own implementation of UpTo; extremely limited *)


If[Definition[UpTo]===Null,
  Unprotect[Take];
  Take[a_,UpTo[n_]]:=If[Length[a]>n,Take[a,n],a];
  Protect[Take];
];


(* pre-10.2: use simple implementations of ContainsExactly and Nothing *)


If[Definition[ContainsExactly]===Null,
  ContainsExactly[a_,b_]:=Union[a]===Union[b];
  Protect[ContainsExactly];
];
If[Definition[Nothing]===Null,
  Unprotect[List];
  List[pre___,Nothing..,post___]:=List[pre,post];
  Protect[List];
];


(* pre-10.1: use simple implementation of KeyValueMap *)


If[Definition[KeyValueMap]===Null,
  KeyValueMap[f_,a_Association]:=f@@@Replace[Normal[a],Rule->List,{1},Heads->True];
  KeyValueMap[f_]:=KeyValueMap[f,#]&;
];


(* data extraction *)


NamedTensorRows[NamedTensor[rowNames_Association,colNames_Association,data_]]:=rowNames;
NamedTensorCols[NamedTensor[rowNames_Association,colNames_Association,data_]]:=colNames;
NamedTensorData[NamedTensor[rowNames_Association,colNames_Association,data_]]:=data;


(* constructors *)


NamedTensor[indexNames_List,data_]:=With[{rank=TensorRank[data]},
  Assert[EvenQ[rank]&&Length[indexNames]===rank/2];
  NamedTensor[AssociationThread[indexNames,Range[1,rank/2]],AssociationThread[indexNames,Range[rank/2 +1,rank]],data]
];
NamedTensor[rowNames_Association,colNames_Association,data_]/;Length[rowNames]+Length[colNames]!=TensorRank[data]:=$Failed;
NamedTensor[<||>, <||>, data_] := data;
NamedTensor/:MakeBoxes[tensor:NamedTensor[rowNames_Association,colNames_Association,data_],form:(StandardForm|TraditionalForm)]:=
  With[{above={{BoxForm`SummaryItem[{"Row indices: ",Keys[rowNames]}],BoxForm`SummaryItem[{"Length: ",Length[rowNames]}]},
                {BoxForm`SummaryItem[{"Column indices: ",Keys[colNames]}],BoxForm`SummaryItem[{"Length: ",Length[colNames]}]},
                {BoxForm`SummaryItem[{"Data type: ",Head[data]}],SpanFromLeft}},
    below=If[Head[data]===SparseArray,{{BoxForm`SummaryItem[{"Specified elements: ",Length@data["NonzeroPositions"]}],BoxForm`SummaryItem[{"Density: ",data["Density"]}]},
                                        {BoxForm`SummaryItem[{"Default: ",data["Background"]}],SpanFromLeft},
                                        {Pane[Row[{Style["Elements: ","SummaryItemAnnotation"],Take[ArrayRules[data],UpTo[4]]}],{300,Automatic},BaselinePosition->Baseline]}},
                                       {{Pane[Row[{Style["Elements: ","SummaryItemAnnotation"],Short[data]}],{300,Automatic},BaselinePosition->Baseline]}}]},
    BoxForm`ArrangeSummaryBox[
      NamedTensor,
      NamedTensor[rowNames,colNames,data],
      None,
      above,
      below,
      form,
      "Interpretable"->Automatic
    ]
  ];
NamedTensor[rowNames_Association,colNames_Association,data_][rowParts_Association,colParts_Association]:=
  With[{mapRow=AssociationThread[Values[rowNames],Keys[rowNames]],mapCol=AssociationThread[Values[colNames],Keys[colNames]],tensorRank=TensorRank[data]},
    With[{indexData=Reap[Table[With[{val=If[KeyExistsQ[mapRow,i],Lookup[rowParts,Key[mapRow[i]],All],Lookup[colParts,Key[mapCol[i]],All]]},If[NumberQ[val],Sow[i]];val],{i,tensorRank}]]},
      With[{deleteNames=If[Length[indexData[[2]]]>0,DeleteCases[Alternatives@@indexData[[2,1]]],Identity]},
        With[{newIndices=AssociationThread[deleteNames[Range[tensorRank]],Range[tensorRank-Length[indexData[[2,1]]]]]},
          NamedTensor[
            newIndices/@deleteNames[rowNames],
            newIndices/@deleteNames[colNames],
            data[[Sequence@@indexData[[1]]]]
          ]
        ]
      ]
    ]
  ];
NamedTensor/:MatrixForm[t_NamedTensor]:=MatrixForm[NamedTensorAsMatrix[t]];
NamedTensor[rowNames_Association,colNames_Association,data_][rowParts_List,colParts_Association]:=NamedTensor[rowNames,colNames,data][Association[rowParts],colParts];
NamedTensor[rowNames_Association,colNames_Association,data_][rowParts_Association,colParts_List]:=NamedTensor[rowNames,colNames,data][rowParts,Association[colParts]];
NamedTensor[rowNames_Association,colNames_Association,data_][rowParts_List,colParts_List]:=NamedTensor[rowNames,colNames,data][Association[rowParts],Association[colParts]];
NamedTensor[rowNames_Association,colNames_Association,data_][parts_]:=NamedTensor[rowNames,colNames,data][parts,parts];
NamedTensor[rowNames_Association,colNames_Association,data_][parts__Rule]:=With[{p=Association[parts]},NamedTensor[rowNames,colNames,data][p,p]];
Protect[NamedTensor];


NamedMatrix[data_/;TensorRank[data]===2]:=NamedTensor[<|""->1|>,<|""->2|>,data];
NamedColVector[data_/;TensorRank[data]===1]:=NamedTensor[<|""->1|>,<||>,data];
NamedRowVector[data_/;TensorRank[data]===1]:=NamedTensor[<||>,<|""->1|>,data];
Protect[NamedMatrix,NamedColVector,NamedRowVector];


RenameIndices[NamedTensor[rowNames_Association,colNames_Association,data_],newNames__]:=With[{renames=Association[newNames]},
  With[{map=KeyMap[Lookup[renames,Key[#],#]&]},
    NamedTensor[map[rowNames],map[colNames],data]
  ]
];


(* representation as matrix *)


NamedTensorAsMatrix[NamedTensor[rowNames_Association,colNames_Association,data_],order_List]/;Union[Keys[rowNames],Keys[colNames]]===Union[order]:=
  With[{rowList=Table[Lookup[rowNames,Key[rowName],Nothing],{rowName,order}],colList=Table[Lookup[colNames,Key[colName],Nothing],{colName,order}]},
    If[
      rowList==={}||colList==={},
      Flatten[data,colList],
      Flatten[data,{rowList,colList}]
    ]
  ];
NamedTensorAsMatrix[NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensorAsMatrix[NamedTensor[rowNames,colNames,data],DeleteDuplicates[Join[Keys[rowNames],Keys[colNames]]]];
Protect[NamedTensorAsMatrix];


NamedTensorMatrixFunction[f_,NamedTensor[rowNames_Association,colNames_Association,data_],order_List]/;Union[Keys[rowNames],Keys[colNames]]===Union[order]:=
  With[{rows=Length[rowNames],cols=Length[colNames],dimensions=TensorDimensions[data]},
    NamedTensor[
      AssociationThread[Select[order,KeyExistsQ[rowNames,#]&],Range[rows]],
      AssociationThread[Select[order,KeyExistsQ[colNames,#]&],Range[rows+1,rows+cols]],
      ArrayReshape[f@NamedTensorAsMatrix[NamedTensor[rowNames,colNames,data],order],dimensions[[Join[Values@rowNames,Values@colNames]]]]
    ]
  ];
NamedTensorMatrixFunction[f_,NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensorMatrixFunction[f,NamedTensor[rowNames,colNames,data],DeleteDuplicates[Join[Keys[rowNames],Keys[colNames]]]];
Protect[NamedTensorMatrixFunction];


(* extend System functions: functions that directly operate on the data *)


Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data1_],NamedTensor[rowNames_Association,colNames_Association,data2_]]:=NamedTensor[rowNames,colNames,f[data1,data2]];
    f[NamedTensor[rowNames1_Association,colNames1_Association,data1_],NamedTensor[rowNames2_Association,colNames2_Association,data2_]]/;ContainsExactly[Keys[rowNames1],Keys[rowNames2]]&&ContainsExactly[Keys[colNames1],Keys[colNames2]]:=
      With[{rowOrder2=AssociationThread[Values[rowNames2],Keys[rowNames2]],colOrder2=AssociationThread[Values[colNames2],Keys[rowNames2]]},
        NamedTensor[rowNames1,colNames1,f[data1,TensorTranspose[data2,Table[If[KeyExistsQ[rowOrder2,i],rowNames1[rowOrder2[i]],colNames1[colOrder2[i]]],{i,TensorRank[data1]}]]]]
      ];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],other_/;Head[other]=!=NamedTensor]:=NamedTensor[rowNames,colNames,f[data,other]];
    Protect[f]
  ],
  {f,{Plus,Minus,Times,Divide}}
];
Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],args___]:=NamedTensor[rowNames,colNames,f[data,args]];
    Protect[f]
  ],
  {f,{Conjugate,Simplify,FullSimplify,SparseArray,Normal}}
];


(* extend System functions: functions that operate on the interpretation of the tensor *)


Unprotect[Transpose,ConjugateTranspose];
Transpose[NamedTensor[rowNames_Association,colNames_Association,data_],indices_List]:=NamedTensor[
  Association[{KeyDrop[rowNames,indices],KeyTake[colNames,indices]}],
  Association[{KeyDrop[colNames,indices],KeyTake[rowNames,indices]}],
  data
];
Transpose[NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensor[colNames,rowNames,data];
Transpose[NamedTensor[rowNames_Association,colNames_Association,data_],index_]:=Transpose[NamedTensor[rowNames,colNames,data],{index}];
ConjugateTranspose[NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensor[colNames,rowNames,Conjugate[data]];
ConjugateTranspose[NamedTensor[rowNames_Association,colNames_Association,data_],indices_]:=Transpose[NamedTensor[rowNames,colNames,Conjugate[data]],indices];
Protect[Transpose,ConjugateTranspose];


(* extend System functions: functions that operate on the matrix representation and give scalars *)


Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],order_List]:=f[NamedTensorAsMatrix[NamedTensor[rowNames,colNames,data],order]];
    f[NamedTensor[rowNames_Association,colNames_Association,data_]]:=f[NamedTensorAsMatrix[NamedTensor[rowNames,colNames,data]]];
    Protect[f]
  ],
{f,{Det,MatrixRank}}];


(* extend System functions: functions that operate on the matrix representation and give matrices *)


Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],order_List]:=NamedTensorMatrixFunction[f,NamedTensor[rowNames,colNames,data],order];
    f[NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensorMatrixFunction[f,NamedTensor[rowNames,colNames,data]];
    Protect[f]
  ],
{f,{Inverse,MatrixPower,MatrixExp,MatrixLog}}];

Unprotect[MatrixFunction];
MatrixFunction[f_,NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensorMatrixFunction[MatrixFunction[f,#]&,NamedTensor[rowNames,colNames,data]];
MatrixFunction[f_,NamedTensor[rowNames_Association,colNames_Association,data_],order_]:=NamedTensorMatrixFunction[MatrixFunction[f,#]&,NamedTensor[rowNames,colNames,data],order];
Protect[MatrixFunction];


(* extend System functions: TensorProducts, TensorDimensions and TensorRank of NamedTensors *)


Unprotect[TensorProduct,TensorDimensions,TensorRank];
TensorProduct[Longest[elements__Rule?(AllTrue[Head[#[[2]]]===NamedTensor&])]]:=TensorProduct[Association[elements]];
TensorProduct[elements_Association?(AllTrue[Head[#]===NamedTensor&])]:=With[{tensorData=NamedTensorData/@Values[elements]},
  With[{tensorRankOffsets=AssociationThread[Keys[elements],Most[FoldList[Plus,0,TensorRank/@tensorData]]]},
    NamedTensor[
      Association[KeyValueMap[Function[{tensorName,tensorRows},KeyMap[If[#==="",tensorName,tensorName<>"."<>#]&,tensorRows]],NamedTensorRows/@elements+tensorRankOffsets]],
      Association[KeyValueMap[Function[{tensorName,tensorRows},KeyMap[If[#==="",tensorName,tensorName<>"."<>#]&,tensorRows]],NamedTensorCols/@elements+tensorRankOffsets]],
      TensorProduct@@tensorData
    ]
  ]
];
TensorProduct[Longest[elements__NamedTensor]]:=With[{tensors={elements}},
  With[{tensorData=NamedTensorData/@tensors},
    With[{tensorRankOffsets=Most[FoldList[Plus,0,TensorRank/@tensorData]]},
      NamedTensor[
        Association[NamedTensorRows/@tensors+tensorRankOffsets],
        Association[NamedTensorCols/@tensors+tensorRankOffsets],
        TensorProduct@@tensorData
      ]
    ]
  ]
];
TensorDimensions[NamedTensor[rowNames_Association,colNames_Association,data_]]:=With[{dim=TensorDimensions[data]},Merge[{rowNames,colNames},Map[dim[[#]]&]]];
TensorRank[NamedTensor[rowNames_Association,colNames_Association,data_]]:=TensorRank[data];
Protect[TensorProduct,TensorDimensions,TensorRank];


(* contractions: Trace and partial trace *)


Unprotect[Tr];
Tr[NamedTensor[rowNames_Association,colNames_Association,data_],indices_List]:=With[{contractionNames=Intersection[Keys[rowNames],Keys[colNames],indices]},
  With[{contractionIndices=Table[{rowNames[i],colNames[i]},{i,contractionNames}]},
    With[{contractionIndicesFlat=Flatten@contractionIndices},
      With[{remainingIndices=Complement[Range[TensorRank[data]],contractionIndicesFlat]},
        With[{newIndices=AssociationThread[remainingIndices,Range[Length[remainingIndices]]]},
          NamedTensor[
            newIndices/@KeyDrop[rowNames,contractionNames],
            newIndices/@KeyDrop[colNames,contractionNames],
            TensorContract[data,contractionIndices]
          ]
        ]
      ]
    ]
  ]
];
Tr[NamedTensor[rowNames_Association,colNames_Association,data_],indices__]:=Tr[NamedTensor[rowNames,colNames,data],{indices}];
Tr[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Tr[NamedTensor[rowNames,colNames,data],Intersection[Keys[rowNames],Keys[colNames]]];
Protect[Tr];


(* contractions: Dot, i.e. automatically contract same names *)


Unprotect[Dot];
Dot[Longest[elements__NamedTensor]]:=With[{tensors={elements}},
  With[{tensorRanks=TensorRank[NamedTensorData[#]]&/@tensors,tensorCount=Length@tensors},
    With[{totalRankRange=Range@Total@tensorRanks},
      With[{indices=TakeList[totalRankRange,tensorRanks]},
        With[{contractionIndices=Table[
          Sequence@@KeyValueMap[
            Function[{colName,colIndex},
              With[{rightTensorIndex=leftTensorIndex+LengthWhile[tensors[[leftTensorIndex+1;;]],!KeyExistsQ[NamedTensorRows[#],colName]&]+1},
                If[rightTensorIndex>tensorCount,
                  Nothing,
                  {indices[[leftTensorIndex,colIndex]],indices[[rightTensorIndex,NamedTensorRows[tensors[[rightTensorIndex]]][colName]]]}
                ]
              ]
            ],NamedTensorCols[tensors[[leftTensorIndex]]]],{leftTensorIndex,Length@tensors-1}]},
           With[{contractionIndicesFlat=Flatten@contractionIndices},
             With[{remainingIndices=Complement[totalRankRange,contractionIndicesFlat]},
               With[{newIndices=AssociationThread[remainingIndices,Range[Length[remainingIndices]]]},
                 With[{rowIndices=MapIndexed[Function[{tensor,tensorPosition},KeyValueMap[If[KeyExistsQ[newIndices,indices[[First[tensorPosition],#2]]],#1->newIndices[indices[[First[tensorPosition],#2]]],Nothing]&,NamedTensorRows[tensor]]],tensors],
                   colIndices=MapIndexed[Function[{tensor,tensorPosition},KeyValueMap[If[KeyExistsQ[newIndices,indices[[First[tensorPosition],#2]]],#1->newIndices[indices[[First[tensorPosition],#2]]],Nothing]&,NamedTensorCols[tensor]]],tensors]},
                  NamedTensor[
                    Association[rowIndices],Association[colIndices],
                    Activate@TensorContract[Inactive[TensorProduct]@@NamedTensorData/@tensors,Sort@contractionIndices]
                  ]
                ]
              ]
            ]
          ]
        ]
      ]
    ]
  ]
];
Protect[Dot];


(* contractions: full command, specify all contractions manually *)


Unprotect[TensorContract];
ClearAll[`TensorContract`tensorName,`TensorContract`indices];
`TensorContract`extract=KeyValueMap[If[#1==="",`TensorContract`tensorName,`TensorContract`tensorName<>"."<>#1]->`TensorContract`indices[`TensorContract`tensorName][[#2]]&];
`TensorContract`replace[index_->name_]:=name->`TensorContract`newIndices[index];
TensorContract[elements_Association?(AllTrue[Head[#]===NamedTensor&]),contractions_List?(AllTrue[VectorQ])]:=With[{tensorRanks=TensorRank[NamedTensorData[#]]&/@Values[elements],tensorNames=Keys[elements]},
  With[{totalRankRange=Range@Total@tensorRanks},
    Block[{`TensorContract`indices=AssociationThread[tensorNames,TakeList[totalRankRange,tensorRanks]]},
      With[{rowIndices=Association@Flatten[Table[`TensorContract`extract[NamedTensorRows[elements[`TensorContract`tensorName]]],{`TensorContract`tensorName,tensorNames}]],
        colIndices=Association@Flatten[Table[`TensorContract`extract[NamedTensorCols[elements[`TensorContract`tensorName]]],{`TensorContract`tensorName,tensorNames}]]},
        With[{contractionIndices=Table[
          Sequence@@Table[{colIndices[`TensorContract`contraction[[i]]],rowIndices[`TensorContract`contraction[[i+1]]]},{i,1,Length[`TensorContract`contraction]-1}],
          {`TensorContract`contraction,contractions}],
           rowIndexNames=AssociationThread[Values[rowIndices],Keys[rowIndices]],colIndexNames=AssociationThread[Values[colIndices],Keys[colIndices]]},
          With[{contractionIndicesFlat=Flatten@contractionIndices},
            With[{remainingIndices=Complement[totalRankRange,contractionIndicesFlat]},
              Block[{`TensorContract`newIndices=AssociationThread[remainingIndices,Range[Length[remainingIndices]]]},
                NamedTensor[
                  AssociationMap[`TensorContract`replace,KeyDrop[rowIndexNames,Flatten@contractionIndices]],
                  AssociationMap[`TensorContract`replace,KeyDrop[colIndexNames,Flatten@contractionIndices]],
                  Activate@TensorContract[Inactive[TensorProduct]@@NamedTensorData/@Values[elements],Sort@contractionIndices]
                ]
              ]
            ]
          ]
        ]
      ]
    ]
  ]
];
TensorContract[elements_Association,contraction_List?VectorQ]:=TensorContract[elements,{contraction}];
TensorContract[elements_List?(AllTrue[Head[#]===NamedTensor&]),contractions_]:=TensorContract[Association[elements],contractions];
TensorContract[elements__Rule?(AllTrue[#[[2]]===NamedTensor&]),contractions_]:=TensorContract[Association[elements],contractions];
TensorContract::invctr="Invalid contraction `1`; head must be List; elements must be valid names";
Protect[TensorContract];


`NamedOperation`adjointKeys=KeyMap[#<>"\[Dagger]"&];
`NamedOperation`newName[s_String]:=With[{split=StringSplit[s,".",2]},If[Length[split]===2,split[[1]]<>"\[Dagger]."<>split[[2]],split[[1]]<>"\[Dagger]"]];
`NamedOperation`extractStateOp[s_String]:=First[StringSplit[s,".",2]];
`NamedOperation`adjointContractions[states_Association]:=Map[
  Sequence@@Table[
    If[KeyExistsQ[states,`NamedOperation`extractStateOp[#[[i]]]],
      If[KeyExistsQ[states,`NamedOperation`extractStateOp[#[[i+1]]]],
        Nothing,
        Message[NamedOperation::invctr,#];Nothing
      ],
      If[KeyExistsQ[states,`NamedOperation`extractStateOp[#[[i+1]]]],
        {#[[i+1]],`NamedOperation`newName[#[[i]]]},
        {`NamedOperation`newName[#[[i+1]]],`NamedOperation`newName[#[[i]]]}
      ]
    ],
    {i,Length[#]-1}
  ]&
];
NamedOperation::invctr="Invalid contraction for a quantum operation: Operator must be to the left of state, but got `1`.";
`NamedOperation`adjointRenames[gateName_->NamedTensor[rowIndices_,colIndices_,data_]]:=Sequence@@(If[#==="",gateName<>"\[Dagger]"->gateName,gateName<>"\[Dagger]."<>#->gateName<>"."<>#]&/@Union[Keys[rowIndices],Keys[colIndices]]);
NamedOperation[gates_Association?(AllTrue[Head[#]===NamedTensor&]),states_Association?(AllTrue[Head[#]===NamedTensor&]),contractions_List?(AllTrue[VectorQ])]/;Intersection[Keys@gates,Keys@states]==={}:=RenameIndices[TensorContract[
  Association[gates,states,ConjugateTranspose/@`NamedOperation`adjointKeys[gates]],
  Join[contractions,`NamedOperation`adjointContractions[states][contractions]]
], AssociationMap[`NamedOperation`adjointRenames,gates]];
NamedOperation[gates_Association,states_Association,contractions_]/;Intersection[Keys@gates,Keys@states]=!={}:=$Failed;
NamedOperation[gates_Association,states_Association,contraction_List?VectorQ]:=NamedOperation[gates,states,{contraction}];
NamedOperation[gates_List,states_Association,contractions_]:=NamedOperation[Association[gates],states,contractions];
NamedOperation[gates_Association,states_List,contractions_]:=NamedOperation[gates,Association[states],contractions];
NamedOperation[gates_List,states_List,contractions_]:=NamedOperation[Association[gates],Association[states],contractions];
NamedOperation[gates__Rule,On,states__Rule,contractions_]:=NamedOperation[Association[gates],Association[states],contractions];
Protect[NamedOperation];


(* some conditionals *)


NamedRowsQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_List]:=Complement[names,Keys[rowNames]]==={};
NamedRowsQ[NamedTensor[rowNames_Association,colNames_Association,data_],name_String]:=KeyExistsQ[rowNames,name];
NamedRowsQ[names_]:=NamedRowsQ[#,names]&;
NamedColumnsQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_List]:=Complement[names,Keys[colNames]]==={};
NamedColumnsQ[NamedTensor[rowNames_Association,colNames_Association,data_],name_String]:=KeyExistsQ[colNames,name];
NamedColumnsQ[names_]:=NamedColumnsQ[#,names]&;
NamedIndicesQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_]:=NamedRowsQ[NamedTensor[rowNames,colNames,data],names]&&NamedColumnsQ[NamedTensor[rowNames,colNames,data],names];
NamedIndicesQ[NamedTensor[rowNames_Association,colNames_Association,data_],namesRow_,namesCol_]:=NamedRowsQ[NamedTensor[rowNames,colNames,data],namesRow]&&NamedColumnsQ[NamedTensor[rowNames,colNames,data],namesCol];
NamedIndicesQ[namesRow_,namesCol_]:=NamedIndicesQ[#,namesRow,namesCol]&;
NamedIndicesQ[names_]:=NamedIndicesQ[#,names,names]&;
Protect[NamedRowsQ,NamedColumnsQ,NamedIndicesQ];


NamedRows[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Keys[rowNames];
NamedColumns[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Keys[colNames];
NamedIndices[NamedTensor[rowNames_Association,colNames_Association,data_]]:=DeleteDuplicates[Join[Keys[rowNames],Keys[colNames]]];
Protect[NamedRows,NamedColumns,NamedIndices];


Unprotect[VectorQ,MatrixQ,ArrayQ,TensorQ];
VectorQ[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Length[rowNames]+Length[colNames]===1;
MatrixQ[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Length[rowNames]===Length[colNames]===1;
ArrayQ[NamedTensor[rowNames_Association,colNames_Association,data_]]:=True;
TensorQ[NamedTensor[rowNames_Association,colNames_Association,data_]]:=True;
Protect[VectorQ,MatrixQ,ArrayQ,TensorQ];


End[];


EndPackage[]
