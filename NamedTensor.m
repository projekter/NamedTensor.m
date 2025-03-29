(* ::Package:: *)

BeginPackage["NamedTensor`"];


Unprotect@@Names["NamedTensor`*"];
ClearAll@@Names["NamedTensor`*"];


NamedTensor::usage="NamedTensor[rowNames,colNames,data] is the general data structure for a named tensor. rowNames and colNames are associations that map the names of the indices (strings) to the indices in the tensor data.

NamedTensor[indexNames,data] is a construction shorthand. indexNames is a list that in this order defines the names of the indices, which are identical for rows and columns.";


NamedMatrix::usage="NamedMatrix[data] is a shorthand that converts a rank-2 tensor into a named tensor with empty names, one row and one column index.";
NamedColumnVector::usage="NamedColumnVector[data] is a shorthand that converts a rank-1 tensor into a named tensor with one row index of empty name.";
NamedRowVector::usage="NamedRowVector[data] is a shorthand that converts a rank-1 tensor into a named tensor with one column index of empty name.";
NamedProjection::usage="NamedProjection[row or column vector] makes vectors into their projections (no normalization is performed!).

NamedProjection[indexNames,data] is a construction shorthand for NamedProjection[NamedTensor[indexNames,{},data]].

NamedProjection[...,Reals] does not conjugate.";


NamedCondition::usage="NamedCondition[posProj,posOp,negProj,negOp] gives posProj\[TensorProduct]posOp+negProj\[TensorProduct]negOp. If either posProj or negProj are Automatic, they are inferred from the other. If either posOp or negOp are Automatic, they are identities. If they are lists or associations, they are wrapped in TensorProduct first.

NamedCondition[posProj,posOp] is a shorthand for NamedCondition[posProj,posOp,Automatic,Automatic]."


RenameIndices::usage="RenameIndices[tensor,replacements] renames the row and column indices in a NamedTensor object according to the rules or association specified (no patterns!)";


NamedTensorAsMatrix::usage="NamedTensorAsMatrix[tensor,order] converts a NamedTensor to a rank-2 tensor (standard List or SparseArray). By default, the levels are flattened together in the order as they appear in their associations, but a custom order can be specified.";
NamedTensorMatrixFunction::usage="NamedTensorMatrixFunction[f,tensor,order] applies a function f to the NamedTensor tensor. Since f expects a standard matrix, it is first converted into a matrix (the order may be specified, else the default order is used). After this, the matrix is converted back to the named tensor. f must not change the dimensions of the matrix.";


MergeRowIndices::usage="MergeRowIndices[tensor,indices->new,...] merges all row indices specified in `indices` into a new row index with name `new`. If `new` is omitted, the name of the first index specified in `indices` is taken instead.";
MergeColumnIndices::usage="MergeColumnIndices[tensor,indices->new,...] merges all column indices specified in `indices` into a new column index with name `new`. If `new` is omitted, the name of the first index specified in `indices` is taken instead.";
MergeIndices::usage="MergeIndices[tensor,indices->new,...] merges all indices specified in `indices` into a new index with name `new`. If `new` is omitted, the name of the first index specified in `indices` is taken instead. All indices must exist both as row and column indices.";
(*SplitRowIndex::usage="SplitRowIndex[tensor,index->{{name1,dimension1},...},...] splits the row index specified in `index` into multiple new indices with given dimensions.";
SplitColumnIndex::usage="SplitColumnIndex[tensor,index->{{name1,dimension1},...},...] splits the columm index specified in `index` into multiple new indices with given dimensions.";
SplitIndex::usage="SplitIndex[tensor,index->{{name1,dimension1},...},...] splits the index specified in `index` into multiple new indices with given dimensions.";*)


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
NamedRowsExactlyQ::usage="NamedRowsExactlyQ[tensor,names] gives True iff tensor contains all rows in names, and only those.

NamedRowsExactlyQ[names] is the corresponding operator form.";
NamedColumnsExactlyQ::usage="NamedColsQ[tensor,names] gives True iff tensor contains all columns in names, and only those.

NamedColumnsExactlyQ[names] is the corresponding operator form.";
NamedIndicesExactlyQ::usage="NamedIndicesExactlyQ[tensor,names] gives True iff tensor contains all rows and all columns in names, and only those.

NamedIndicesExactlyQ[names] is the corresponding operator form.

NamedIndicesExactlyQ[tensor,namesRow,namesCol] allows to specify different names for rows and columns.

NamedIndicesExactlyQ[namesRow,namesCol] is the corresponding operator form.";
NamedRows::usage="NamedRows[tensor] lists all row index names in tensor.";
NamedColumns::usage="NamedColumns[tensor] lists all column index names in tensor.";
NamedIndices::usage="NamedIndices[tensor] lists all index names in tensor.";


Begin["`Private`"];


(* pre-11.2: use own implementation for TakeList; very reduced and simplistic variant, hence in our private context *)


If[ToString@Definition[TakeList]==="Null",
  TakeList[list_,counts_List]:=With[{acounts=FoldList[Plus,0,counts]},Table[list[[acounts[[i-1]]+1;;acounts[[i]]]],{i,2,Length@acounts}]];
  SetAttributes[TakeList,NHoldRest];
  Protect[TakeList];
];


(* pre-10.3: use own implementation of UpTo; extremely limited *)


If[ToString@Definition[UpTo]==="Null",
  Unprotect[Take];
  Take[a_,UpTo[n_]]:=If[Length[a]>n,Take[a,n],a];
  Protect[Take];
];


(* pre-10.2: use simple implementations of ContainsExactly, ContainsAll and Nothing *)


If[ToString@Definition[ContainsExactly]==="Null",
  ContainsExactly[a_,b_]:=Union[a]===Union[b];
  Protect[ContainsExactly];
];
If[ToString@Definition[ContainsAll]==="Null",
  ContainsAll[a_,b_]:=Complement[b,a]==={};
  Protect[ContainsAll];
];
If[ToString@Definition[Nothing]==="Null",
  Unprotect[List];
  List[pre___,Nothing..,post___]:=List[pre,post];
  Protect[List];
];


(* pre-10.1: use simple implementation of KeyValueMap *)


If[ToString@Definition[KeyValueMap]==="Null",
  KeyValueMap[f_,a_Association]:=f@@@Replace[Normal[a],Rule->List,{1},Heads->True];
  KeyValueMap[f_]:=KeyValueMap[f,#]&;
];


(* this function is applied on the indices before Mathematica's TensorContract is called - you may want to change it to Echo... *)


`$ContractionMap=Identity;


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
NamedTensor[<||>, <||>, data_]:=data;
NamedTensor/:MakeBoxes[tensor:NamedTensor[rowNames_Association,colNames_Association,data_],form:(StandardForm|TraditionalForm)]:=
  With[{dims=TensorDimensions[data]},
    With[{above={{BoxForm`SummaryItem[{"Row dimensions: ",AssociationThread[Keys[rowNames],dims[[Values[rowNames]]]]}],BoxForm`SummaryItem[{"Length: ",Length[rowNames]}]},
                  {BoxForm`SummaryItem[{"Column dimensions: ",AssociationThread[Keys[colNames],dims[[Values[colNames]]]]}],BoxForm`SummaryItem[{"Length: ",Length[colNames]}]},
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


NamedMatrix[data_/;TensorRank[data]===2]:=NamedTensor[<|""->1|>,<|""->2|>,data];
NamedColumnVector[data_/;TensorRank[data]===1]:=NamedTensor[<|""->1|>,<||>,data];
NamedRowVector[data_/;TensorRank[data]===1]:=NamedTensor[<||>,<|""->1|>,data];
NamedProjection[NamedTensor[rowNames_Association,<||>,data_]]:=NamedTensor[rowNames,rowNames+Length[rowNames],TensorProduct[data,Conjugate[data]]];
NamedProjection[NamedTensor[<||>,colNames_Association,data_]]:=NamedTensor[colNames,colNames+Length[colNames],TensorProduct[Conjugate[data],data]];
NamedProjection[NamedTensor[rowNames_Association,<||>,data_],Reals]:=NamedTensor[rowNames,rowNames+Length[rowNames],TensorProduct[data,data]];
NamedProjection[NamedTensor[<||>,colNames_Association,data_],Reals]:=NamedTensor[colNames,colNames+Length[colNames],TensorProduct[data,data]];
NamedProjection[indexNames_List,data_,rest___]:=With[{rank=TensorRank[data]},
  Assert[Length[indexNames]===rank];
  NamedProjection[NamedTensor[AssociationThread[indexNames,Range[rank]],<||>,data],rest]
];
Protect[NamedMatrix,NamedColumnVector,NamedRowVector,NamedProjection];


NamedCondition[posProj:NamedTensor[posProjRowNames_Association,posProjColNames_Association,_],posOp:NamedTensor[posOpRowNames_Association,posOpColNames_Association,_],
  negProj:NamedTensor[negProjRowNames_,negProjColNames_,_],negOp:NamedTensor[negOpRowNames_,negOpColNames_,_]]/;
  ContainsExactly[Keys[posProjRowNames],Keys[negProjRowNames]]&&ContainsExactly[Keys[posProjColNames],Keys[negProjColNames]]&&
  ContainsExactly[Keys[posOpRowNames],Keys[negOpRowNames]]&&ContainsExactly[Keys[posOpColNames],Keys[negOpColNames]]:=posProj\[TensorProduct]posOp+negProj\[TensorProduct]negOp;
NamedCondition[posProj_,posOp:NamedTensor[posOpRowNames_Association,posOpColNames_Association,posOpData_],negProj_,Automatic]:=
  With[{shift=Length[posOpRowNames],dims=TensorDimensions[posOpData]},
    NamedCondition[posProj,posOp,negProj,
      NamedTensor[
        AssociationThread[Keys[posOpRowNames],Range[shift]],
        AssociationThread[Keys[posOpColNames],Range[shift+1,shift+Length[posOpColNames]]],
        ArrayReshape[IdentityMatrix[Sqrt[Times@@dims],If[Head[posOpData]===SparseArray,SparseArray,List]],dims]
      ]
    ]
  ];
NamedCondition[posProj:NamedTensor[posProjRowNames_Association,posProjColNames_Association,posProjData_],posOp_,Automatic,negProj_]:=
  With[{shift=Length[posProjRowNames],dims=TensorDimensions[posProjData]},
    NamedCondition[posProj,posOp,
      NamedTensor[
        AssociationThread[Keys[posProjRowNames],Range[shift]],
        AssociationThread[Keys[posProjColNames],Range[shift+1,shift+Length[posProjColNames]]],
        ArrayReshape[IdentityMatrix[Sqrt[Times@@dims],If[Head[posProjData]===SparseArray,SparseArray,List]],dims]
      ]-posProj,negProj
    ]
  ];
NamedCondition[posProj_,Automatic,negProj_,negOp_NamedTensor]:=NamedCondition[negProj,negOp,posProj,Automatic];
NamedCondition[Automatic,posOp_,negProj_NamedTensor,negOp_]:=NamedCondition[negProj,negOp,Automatic,posOp];
NamedCondition[posProj_,posOp_,negProj_,negOp_]:=With[{heads=Union[Head/@{posProj,posOp,negProj,negOp}]},
  NamedCondition[TensorProduct[posProj],TensorProduct[posOp],TensorProduct[negProj],TensorProduct[negOp]]/;ContainsAll[{Association,List,NamedTensor},heads]&&heads=!={NamedTensor}
]
NamedCondition[posProj_,posOp_]/;ContainsAll[{Association,List,NamedTensor},Union[Head/@{posProj,posOp}]]:=
  NamedCondition[TensorProduct[posProj],TensorProduct[posOp],Automatic,Automatic];
Protect[NamedCondition];


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
      Flatten[data,Join[rowList,colList]],
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


(* reshaping the tensor *)


MergeRowIndices[NamedTensor[rowNames_Association,colNames_Association,data_],merges_List/;VectorQ[merges,Head[#]===Rule&]]:=
  With[{mergeLists=merges[[All,1]],mergeInto=merges[[All,2]],allMerges=Flatten[merges[[All,1]]]},
    If[ContainsAll[Keys@rowNames,allMerges]&&DuplicateFreeQ[allMerges]&&ContainsNone[KeyDrop[rowNames,allMerges],mergeInto],
      With[{allMergeValues=Values[rowNames[[Key/@allMerges]]],mergeLen=Length[merges]},
        NamedTensor[
          Association@KeyValueMap[
            With[{mergeListPosition=FirstPosition[mergeLists,#1]},
              Which[
                MissingQ[mergeListPosition],
                #1->#2-Count[allMergeValues,_?(LessThan[#2])]+mergeLen,(* this index was not merged; keep it - Flatten puts the merged levels to the beginning *)
                mergeListPosition[[2]]==1,
                mergeInto[[mergeListPosition[[1]]]]->mergeListPosition[[1]],(* this was the first index mentioned in the merge sublist, we keep it and rename it *)
                True,
                Nothing(* all other indices are deleted *)
              ]
            ]&,
            rowNames
          ],
          #-Count[allMergeValues,_?(LessThan[#])]+mergeLen&/@colNames,
          Flatten[data,Map[rowNames,mergeLists,{2}]]
        ]
      ],
      $Failed
    ]
  ];
MergeRowIndices[t_NamedTensor,merges_List/;VectorQ[merges,Head[#]===String&]]:=MergeRowIndices[t,{merges->First[merges]}];
MergeRowIndices[t_NamedTensor,merges__]:=MergeRowIndices[t,Switch[Head[#],Rule,#,List,#->First[#]]&/@{merges}];
MergeColumnIndices[NamedTensor[rowNames_Association,colNames_Association,data_],merges_List/;VectorQ[merges,Head[#]===Rule&]]:=
  With[{mergeLists=merges[[All,1]],mergeInto=merges[[All,2]],allMerges=Flatten[merges[[All,1]]]},
    If[ContainsAll[Keys@colNames,allMerges]&&DuplicateFreeQ[allMerges]&&ContainsNone[KeyDrop[colNames,allMerges],mergeInto],
      With[{allMergeValues=Values[colNames[[Key/@allMerges]]],mergeLen=Length[merges]},
        NamedTensor[
          #-Count[allMergeValues,_?(LessThan[#])]+mergeLen&/@rowNames,
          Association@KeyValueMap[
            With[{mergeListPosition=FirstPosition[mergeLists,#1]},
              Which[
                MissingQ[mergeListPosition],
                #1->#2-Count[allMergeValues,_?(LessThan[#2])]+mergeLen,(* this index was not merged; keep it - Flatten puts the merged levels to the beginning *)
                mergeListPosition[[2]]==1,
                mergeInto[[mergeListPosition[[1]]]]->mergeListPosition[[1]],(* this was the first index mentioned in the merge sublist, we keep it and rename it *)
                True,
                Nothing(* all other indices are deleted *)
              ]
            ]&,
            colNames
          ],
          Flatten[data,Map[colNames,mergeLists,{2}]]
        ]
      ],
      $Failed
    ]
  ];
MergeColumnIndices[t_NamedTensor,merges_List/;VectorQ[merges,Head[#]===String&]]:=MergeColumnIndices[t,{merges->First[merges]}];
MergeColumnIndices[t_NamedTensor,merges__]:=MergeColumnIndices[t,Switch[Head[#],Rule,#,List,#->First[#]]&/@{merges}];
MergeIndices[t_NamedTensor,merges__]:=MergeRowIndices[MergeColumnIndices[t,merges],merges];
Protect[MergeRowIndices,MergeColumnIndices,MergeIndices];


(* extend System functions: functions that directly operate on the data *)


Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data1_],NamedTensor[rowNames_Association,colNames_Association,data2_]]:=NamedTensor[rowNames,colNames,f[data1,data2]];
    f[NamedTensor[rowNames1_Association,colNames1_Association,data1_],NamedTensor[rowNames2_Association,colNames2_Association,data2_]]/;ContainsExactly[Keys[rowNames1],Keys[rowNames2]]&&ContainsExactly[Keys[colNames1],Keys[colNames2]]:=
      With[{rowOrder2=AssociationThread[Values[rowNames2],Keys[rowNames2]],colOrder2=AssociationThread[Values[colNames2],Keys[colNames2]]},
        NamedTensor[rowNames1,colNames1,f[data1,TensorTranspose[data2,Table[If[KeyExistsQ[rowOrder2,i],rowNames1[rowOrder2[i]],colNames1[colOrder2[i]]],{i,TensorRank[data1]}]]]]
      ];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],other_/;Head[other]=!=NamedTensor]:=NamedTensor[rowNames,colNames,f[data,other]];
    f[nt1:NamedTensor[rowNames1_Association,colNames1_Association,data1_],nt2:NamedTensor[rowNames2_Association,colNames2_Association,data2_]]:=
      With[{dim1=Dimensions[data1],dim2=Dimensions[data2],missing1=Complement[Keys[rowNames2],Keys[rowNames1]],missing2=Complement[Keys[rowNames1],Keys[rowNames2]]},
        With[{missing1Dimensions=dim2[[Values[rowNames2[[Key/@missing1]]]]],missing2Dimensions=dim1[[Values[rowNames1[[Key/@missing2]]]]]},
          f[If[missing1==={},nt1,TensorProduct[nt1,NamedTensor[missing1,ArrayReshape[IdentityMatrix[Times@@missing1Dimensions,Head[data1]],Join[missing1Dimensions,missing1Dimensions]]]]],
            If[missing2==={},nt2,TensorProduct[nt2,NamedTensor[missing2,ArrayReshape[IdentityMatrix[Times@@missing2Dimensions,Head[data2]],Join[missing2Dimensions,missing2Dimensions]]]]]]
          /;missing1Dimensions===dim2[[Values[colNames2[[Key/@Complement[Keys[colNames2],Keys[colNames1]]]]]]]&&missing2Dimensions===dim1[[Values[colNames1[[Key/@Complement[Keys[colNames1],Keys[colNames2]]]]]]]
        ]
      ];
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
  {f,{Simplify,FullSimplify,SparseArray,Normal,N,Chop,SetPrecision,SetAccuracy,
      IntegerPart,FractionalPart,Round,Floor,Ceiling,Rationalize,
      RealSign,UnitStep,RealAbs,Clip,Rescale,
      Re,Im,Conjugate,Abs,Arg,Sign,
      Boole}}
];
(* for SparseArrays, Rationalize doesn't work as expected... *)
Unprotect[Rationalize];
Rationalize[SparseArray[Automatic,dim_,background_,{1,csr_,data_}]]:=SparseArray[Automatic,dim,Rationalize[background],{1,csr,Rationalize[data]}];
Rationalize[SparseArray[Automatic,dim_,background_,{1,csr_,data_}],dx_]:=SparseArray[Automatic,dim,Rationalize[background,dx],{1,csr,Rationalize[data,dx]}];
Protect[Rationalize];


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
Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data1_],NamedTensor[rowNames_Association,colNames_Association,data2_]]:=f[data1,data2];
    f[NamedTensor[rowNames1_Association,colNames1_Association,data1_],NamedTensor[rowNames2_Association,colNames2_Association,data2_]]/;ContainsExactly[Keys[rowNames1],Keys[rowNames2]]&&ContainsExactly[Keys[colNames1],Keys[colNames2]]:=
      With[{rowOrder2=AssociationThread[Values[rowNames2],Keys[rowNames2]],colOrder2=AssociationThread[Values[colNames2],Keys[colNames2]]},
        f[data1,TensorTranspose[data2,Table[If[KeyExistsQ[rowOrder2,i],rowNames1[rowOrder2[i]],colNames1[colOrder2[i]]],{i,TensorRank[data1]}]]]
      ];
    f[NamedTensor[rowNames1_Association,colNames1_Association,data1_],NamedTensor[rowNames2_Association,colNames2_Association,data2_]]/;!ContainsExactly[Keys[rowNames1],Keys[rowNames2]]||!ContainsExactly[Keys[colNames1],Keys[colNames2]]=False;
    Protect[f]
  ],
  {f,{Equal,SameQ}}
];


(* extend System functions: functions that operate on the matrix representation and give matrices *)


Unprotect[Inverse];
Inverse[NamedTensor[rowNames_Association,colNames_Association,data_],order_List]/;Union[Keys[rowNames],Keys[colNames]]===Union[order]:=
  With[{rows=Length[rowNames],cols=Length[colNames],dimensions=TensorDimensions[data]},
    NamedTensor[
      AssociationThread[Select[order,KeyExistsQ[colNames,#]&],Range[cols]],
      AssociationThread[Select[order,KeyExistsQ[rowNames,#]&],Range[cols+1,cols+rows]],
      ArrayReshape[Inverse@NamedTensorAsMatrix[NamedTensor[rowNames,colNames,data],order],dimensions[[Join[Values@colNames,Values@rowNames]]]]
    ]
  ];
Inverse[NamedTensor[rowNames_Association,colNames_Association,data_]]:=Inverse[NamedTensor[rowNames,colNames,data],DeleteDuplicates[Join[Keys[rowNames],Keys[colNames]]]];
Protect[Inverse];
Do[
  With[{f=f},
    Unprotect[f];
    f[NamedTensor[rowNames_Association,colNames_Association,data_],order_List]:=NamedTensorMatrixFunction[f,NamedTensor[rowNames,colNames,data],order];
    f[NamedTensor[rowNames_Association,colNames_Association,data_]]:=NamedTensorMatrixFunction[f,NamedTensor[rowNames,colNames,data]];
    Protect[f]
  ],
{f,{MatrixPower,MatrixExp,MatrixLog}}];

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
            TensorContract[data,`$ContractionMap[contractionIndices]]
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
                    Activate@TensorContract[Inactive[TensorProduct]@@NamedTensorData/@tensors,`$ContractionMap[Sort@contractionIndices]]
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
                  Activate@TensorContract[Inactive[TensorProduct]@@NamedTensorData/@Values[elements],`$ContractionMap[Sort@contractionIndices]]
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
NamedRowsExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_List]:=ContainsExactly[Keys[rowNames],names];
NamedRowsExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],name_String]:=Keys[rowNames]==={name};
NamedRowsExactlyQ[names_]:=NamedRowsExactlyQ[#,names]&;
NamedColumnsExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_List]:=ContainsExactly[Keys[colNames],names];
NamedColumnsExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],name_String]:=Keys[colNames]==={name};
NamedColumnsExactlyQ[names_]:=NamedColumnsExactlyQ[#,names]&;
NamedIndicesExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],names_]:=NamedRowsExactlyQ[NamedTensor[rowNames,colNames,data],names]&&NamedColumnsExactlyQ[NamedTensor[rowNames,colNames,data],names];
NamedIndicesExactlyQ[NamedTensor[rowNames_Association,colNames_Association,data_],namesRow_,namesCol_]:=NamedRowsExactlyQ[NamedTensor[rowNames,colNames,data],namesRow]&&NamedColumnsExactlyQ[NamedTensor[rowNames,colNames,data],namesCol];
NamedIndicesExactlyQ[namesRow_,namesCol_]:=NamedIndicesExactlyQ[#,namesRow,namesCol]&;
NamedIndicesExactlyQ[names_]:=NamedIndicesExactlyQ[#,names,names]&;
Protect[NamedRowsQ,NamedColumnsQ,NamedIndicesQ,NamedRowsExactlyQ,NamedColsExactlyQ,NamedIndicesExactlyQ];


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


Protect[NamedTensor];


End[];


EndPackage[]
