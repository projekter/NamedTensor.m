# NamedTensor Mathematica package

Requires at least Mathematica 10.0.

License: LPPL 1.3c

Support the development:
- [![PayPal](https://img.shields.io/badge/donate-via%20PayPal-blue.svg?style=flat)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=UTR3MRBYJ825A&source=url)
- ![Bitcoin](https://img.shields.io/badge/donate-BTC-blue.svg?style=flat) 3KBFpoJuA4eSPLGXEf3jicqaV1czhK36fH
- ![Ethereum](https://img.shields.io/badge/donate-ETH-blue.svg?style=flat) 0xE0F774221290b1E41ea62c2dd9af5dbD3df7c685

## Purpose
This package introduces _named_ tensor objects. While Mathematica fully supports tensors in the form of nested lists since version 9.0, performing contractions of large tensors is cumbersome: each index has to be addressed by its position. This is a time-consuming and error-prone process. This package solves the issue by introducing an object that assigns names to every index. In automatically distinguishes row and column indices to allow automatic contractions.

## Creating named tensors
We start by including the package:
```Mathematica
<<NamedTensor`
```
If we already have ordinary matrices, row, or column vectors, those can be automatically turned into named objects:
```Mathematica
nMat = NamedMatrix[myMatrix];
nColVec = NamedColVector[myVec];
nRowVec = NamedRowVector[myVec];
```
Here, `myMatrix` was an existing tensor of rank 2 (a list of lists), `myVec` was a rank-1 tensor (a list). Both might as well have been sparse objects. In fact, it is possible to convert between sparse and full form in the usual manner:
```Mathematica
nMatSparse = SparseArray@nMat;
nColVecDense = Normal@nColVec;
```
The wrappers above turn their arguments into a `NamedTensor` object. Since there is at most one row and one column in these cases, no names have to be specified.

Now we want to create a named tensor that is not a matrix or a vector. For example, let us define a CNot gate:
```Mathematica
nCNot = NamedTensor[{"C", "T"}, 
   SparseArray@{{{{1, 0}, {0, 0}}, {{0, 1}, {0, 0}}}, {{{0, 0}, {0, 1}}, {{0, 0}, {1, 0}}}}];
```
Here, we used another wrapper function: `NamedTensor[list_List, data_]` will create a named tensor with the indices named as given in the first parameter, and it will create identical names for rows and columns. The `data` parameter must be a tensor with all row indices first followed by all column indices.

The most explicit way of constructing the named tensor object would be to specify all indices manually. Equivalent to the last statement is
```Mathematica
nCNot = NamedTensor[<|"C" -> 1, "T" -> 2|>, <|"C" -> 3, "T" -> 4|>,
   SparseArray@{{{{1, 0}, {0, 0}}, {{0, 1}, {0, 0}}}, {{{0, 0}, {0, 1}}, {{0, 0}, {1, 0}}}}];
```
This is the standard form for any `NamedTensor` object: `NamedTensor[rowIndices_Association, colIndices_Association, data_]`. For every row and every column index, we have to specify a name and the corresponding index. Of course, no index must occur more than once. All indices must for a contiguous range that spans from 1 to `TensorRank[data]`.

Note that *tensor names must always be strings* and they must not contain the dot `"."`!

## Querying named tensors
You can extract information from `NamedTensor` objects by using standard Mathematica functions: `TensorDimensions`, `TensorRank`, `VectorQ`, and `MatrixQ` are overloaded to work on those objects. `ArrayQ` and the undocumented `TensorQ` always give  `True`. Note that `TensorDimensions` returns an association that for every index name as key contains a list with the dimensions of its row and its column dimension.

Additionally, the conditionals `NamedRowsQ[_NamedTensor, names_]`, `NamedColumnsQ[_NamedTensor, names_]`, `NamedIndicesQ[_NamedTensor, namesRowAndCol_]`, and `NamedIndicesQ[_NamedTensor, namesRow_, namesCol_]` give `True` if and only if the given names exist in the tensor. `names` may either be a single string or a list of strings. All conditionals also exist in operator form, e.g., `NamedRowsQ[names_]` to be applied to a `NamedTensor`.

The list of row index names can be obtained by `NamedRows[_NamedTensor]`, likewise `NamedColumns[_NamedTensor]` and `NamedIndices[_NamedTensor]` for columns and all indices. The order of those lists correspond to the order with which the tensors were created. `NamedIndices` lists all rows first, then all columns that did not occur as a row already.

## Working with named tensors
Many of the functions commonly available in Mathematica also work with `NamedTensor` objects.

### Threadable functions
The functions `Plus`, `Minus`, `Times` and `Divide` can be used either with two `NamedTensor` objects of the same structure. Here, same structure means that the objects must at least have the same row and column index names. If their internal order does not match, they are automatically reordered.

The functions `Conjugate`, `Simplify`, `FullSimplify`, `SparseArray`, and `Normal` directly work on the data associated with the tensor.

### Structure functions

The functions `Transpose` and `ConjugateTranspose` swap row and column indices with the same names (without any operation on the data). They may additionally be supplied with a list of index names to be swapped; if none is given, all possible swappings are performed. Note that `TensorTranspose` does not work on `NamedTensor` objects, and in fact it does not make sense, as the actual index order is irrelevant.

The helper function `RenameIndices[_NamedTensor, newNames__]`, where `newNames` is either an association, a list of rules or a single rule allows to rename indices. This is particularly useful after contractions.

Finally, we can also take parts of a `NamedTensor`, though `Part` itself does not work. But similarly to associations, the `NamedTensor` itself can be invoked as a function with two arguments. Those arguments must be associations (or lists of rules) that specify filers: For every index, only the given elements are taken. It is possible to specify slices or keywords such as All. If a single entry is taken, this index is removed from the output tensor. The first argument specifies row indices, the second one column indices. If only one is given (or a sequence of rules), both rows and column are targeted.

### Matrix functions
Lots of functions are defined for matrices only. However, since `NamedTensor` objects know about row and column indices, it is easy to change from the tensor to the associated matrix representation. The function `NamedTensorAsMatrix[_NamedTensor, order_List]` converts the data into a matrix or a vector, i.e., a tensor of rank at most 2. By default, the order initially specified at creation is used (which can be queried by `NamedIndices`), but a different order can be set. If the `MatrixForm` of a `NamedTensor` is requested, `NamedTensorAsMatrix` is internally called.

Since it happens quite often that algorithms to apply a function only exist for matrix forms, `NamedTensorMatrixFunction[f_, _NamedTensor, order_List]` first converts the tensor into a matrix (using a custom index order, if desired), then applies `f` to the matrix and converts the result - which must not change the shape of the matrix - back to a named tensor with the same indices.

Those two functions are internally used to make `Det`, `MatrixRank`, `Inverse`, `MatrixPower`<sup>1</sup>, `MatrixExp`<sup>1</sup>, `MatrixLog` and `MatrixFunction` work on `NamedTensor`s.

<sup>1</sup> Without vector as additional argument.

### Contractions
Contractions are the most important argument for named tensors. The simplest contraction is the one without summation, i.e., a pure tensor product. `TensorProduct` is overloaded in the following ways:

- `TensorProduct[_Association]` expects an association of string keys with `NamedTensor` values. It forms the tensor product, gives back a new `NamedTensor` object, and the indices in this new object will have the names they had before, prefixed with the name of their tensor (which is the key in the association) and connected with a dot<sup>1</sup>.
- `TensorProduct[__Rule]` is the very same thing, but allows to drop the `<||>` association markers.
- `TensorProduct[__NamedTensor]` is a third form that can only be used if the indices in all tensors have unique names (row and column indices separately). The resulting object will contain unprefixed indices.

The next contraction is the trace: all row and column indices of the same name are contracted. Indeed, `Tr` is overloaded to achieve exactly this. Additionally, as a second parameter, a list or sequence of indices may be given, which are then the only indices to be traced over (partial trace).

For matrices, the most common contraction is `Dot`. It is overloaded so that it is possible to perform the `Dot` operation with `NamedTensor` objects. A column index will be automatically contracted with the next available row index of the same name. After this, the name is free to appear again, but it is not permitted to have two row (column) indices of the same name with no matching column (row) index in between.

While in principle, `Dot` and `RenameIndices` together sufficient to perform arbitrary contractions, this is still not comfortable. Most power is provided by `TensorContract`. It is overloaded accept an association (or a list or sequence of rules) as first parameter, which, as with `TensorProduct`, defines names for the individual tensors - they will work as prefixes<sup>1</sup> for the index names. The second parameter, as in Mathematica's version of `TensorContract`, is a list of lists. The individual sub-lists must now contain the index names to be contracted. This is different from Mathematica's `TensorContract` in multiple respects:

- Of course, now the indices are strings and not numbers.
- More than two indices may be specified to perform multiple contractions. This is done by splitting the list into overlapping pairs of two.
- The order is important: The first index of a pair is a column index, the second index of a pair is a row index. Consequently, if for example three indices are given, the first column index is contracted with the second row index; and the second column index is contracted with the third row index.

<sup>1</sup> As a special case, the empty index name can (only) be accessed by the name of the tensor itself without a dot.

### Quantum operations

In quantum physics, a very common operation is `A ϱ A†`: One or multiple operators act from the left on a density matrix, and their adjoints act on the right. To simplify this process `NamedOperation` is provided. It is very similar to `TensorContract`, apart from the fact that the first parameter of `TensorContract` is now split into two parameters: one for the operators, one for the density tensors. The adjoints of the operators do not need to be specified. `NamedOperation` will automatically construct both the adjoints and their corresponding contraction instruction, then call `TensorContract` and finally rename the remaining indices of the adjoints so that they have the same name as the original operators. For this to work, every column index of the operators must appear in a contraction, so that only row indices are left over.

## Examples
### Tensor product of Pauli matrices
We define a function `nPauli[a__]` that creates a named tensor of Pauli matrices (here in sparse form). The names are `"1"`, `"2"`, ... until the number of parameters specified.
```Mathematica
nPauli[a__] := 
  With[{count = Length@{a}}, With[{names = ToString /@ Range[count]},
    NamedTensor[
      AssociationThread[names, Range[1, 2 count, 2]], 
      AssociationThread[names, Range[2, 2 count, 2]], 
      TensorProduct @@ SparseArray @* PauliMatrix /@ {a}]
    ]
  ];
```

### Projection onto a pure state
A pure state can be represented as a column vector. The corresponding projection can be written as the tensor product with the conjugate transpose:
```Mathematica
nProjection[state_NamedTensor /; NamedColumns[state] === {}] := TensorProduct[state, ConjugateTranspose[state]]
```

With this, we can first define the four Bell states as pure states and then build the density matrix of a Bell-diagonal state:
```Mathematica
nBell[0] = NamedTensor[<|"A" -> 1, "B" -> 2|>, <||>, {{1, 0}, {0, 1}} / Sqrt[2]];
nBell[1] = NamedTensor[<|"A" -> 1, "B" -> 2|>, <||>, {{1, 0}, {0, -1}} / Sqrt[2]];
nBell[2] = NamedTensor[<|"A" -> 1, "B" -> 2|>, <||>, {{0, 1}, {1, 0}} / Sqrt[2]];
nBell[3] = NamedTensor[<|"A" -> 1, "B" -> 2|>, <||>, {{0, 1}, {-1, 0}} / Sqrt[2]];
nBellDiagonal[coefficients_List /; Length[coefficients] === 4] :=
  Sum[coefficients[[i]] nProjection[nBell[i -1]], {i, 4}]
```

### Partial transposition and negativity
We can easily calculate the negativity of a state. For this, we need to perform its partial transposition with respect to one of the two subsystems, and we also need to calculate the nuclear norm of this partially transposed state.
```Mathematica
nPartialTranspose[state_NamedTensor] := Block[{subsystems},
  Transpose[state, subsystems[[1]]]
  /; Length[subsystems = NamedIndices[state]] === 2
]

nNuclearNorm[state_NamedTensor] := Total@SingularValueList@NamedTensorAsMatrix[state]
```
Since `SingularValueList` is not among the overloaded functions, we first need to convert to a matrix. Here, we queried the names of the indices used in `nPartialTranspose` so that this function works regardless of which names are chosen for the subsystems.

With this, we can calculate the relevant part of the negativity of a Bell diagonal state:
```Mathematica
Assuming[a ∈ Reals && b ∈ Reals && c ∈ Reals && d ∈ Reals,
  Simplify@nNuclearNorm@nPartialTranspose@nBellDiagonal@{a,b,c,d}
]
```
which will give
```Mathematica
1/2 (Abs[a + b + c - d] + Abs[a + b - c + d] + Abs[a - b + c + d] + Abs[-a + b + c + d])
```

### Entanglement distillation
We will implement the DEJMPS entanglement distillation protocol. It first performs a bilateral rotation on both Alice and Bob's side, then local CNot and finally coincidence measurements on the targets of the CNots. First, we define the local rotation:
```Mathematica
nDEJMPS$rotation = NamedMatrix[{{1, I}, {I, 1}} / Sqrt[2]];
```
Alice will directly apply this, Bob will apply its conjugate.

Then we directly implement the DEJMPS scheme as sketched above.
```Mathematica
nDEJMPS[ϱAB1_NamedTensor?(NamedIndicesQ[{"A", "B"}]), ϱAB2_NamedTensor?(NamedIndicesQ[{"A", "B"}])] :=
  RenameIndices[
    Sum[
      NamedOperation[
        <|
          "rA1" -> nDEJMPS$rotation, "rA2" -> nDEJMPS$rotation,
          "rB1" -> Conjugate@nDEJMPS$rotation, "rB2" -> Conjugate@nDEJMPS$rotation,
          "cnotA" -> nCNot, "cnotB" -> nCNot
        |>,
        <|"ϱ1" -> ϱAB1, "ϱ2" -> ϱAB2 |>,
        {{"cnotA.C", "rA1", "ϱ1.A"}, {"cnotA.T", "rA2", "ϱ2.A"},
         {"cnotB.C", "rB1", "ϱ1.B"}, {"cnotB.T", "rB2", "ϱ2.B"}}
      ]["cnotA.T" -> result, "cnotB.T" -> result],
      {result, 2}
    ],
    "cnotA.C" -> "A", "cnotB.C" -> "B"
  ]
```
Without writing one index, we were directly able to implement the plain text description.

- We define a function `nDEJMPS` which takes two `NamedTensor` objects as parameters. To make the interface clear, we specify that these objects must have both row and column indices that are called `"A"` and `"B"`.
- Then we start reading from the inner. We apply a quantum operation that consists of two rotations, two conjugate rotations, and also two CNots (defined at the beginning of this Readme).
- The quantum operation acts on two states to which we assign the names `"ϱ1"` and `"ϱ2"`.
- The contractions are as follows - and now we must keep in mind that this is the common contraction order, i.e., the temporal order is reversed -: Alice's part of `"ϱ1"` is rotated; the output of the rotation is fed into the `"C"` (control) index of a CNot gate. Also Alice's part of `"ϱ2"` is rotated and passed on to the target index of the same CNot gate. We do the analogous actions on Bob's side.
- We subselect on the target indices having the same value.
- We sum over the results (here, an optimization would be to just take one result and multiply by two, since both summation entries give the same tensor anyways.
- Finally, we don't want the output to give any conclusion about the internals. For this reason, we rename the remaining indices, `"cnotA.c"` and `"cnotB.C"` into `"A"` and `"B"`.

We can check this out:
```Mathematica
result = Simplify[nDEJMPS@@ConstantArray[nBellDiagonal[{a, b, c, d}], 2]];
Tr[result]
MatrixForm[result]
```
This will give the results `(b + c)^2 + (a + d)^2` as success probability and the resulting density matrix in the computational basis.

## Notes
The use of these objects introduces some overhead; however, this should be barely noticeable. Most of the time should be spent in doing the actual operations. The internal calls to Mathematica's `TensorContract` use the efficient inactivated variants of `TensorProduct` and sort the indices for optimum performance. Compared to the time usually spent on figuring out indices or bugtracking them, the overhead is negative.

As the last example showed, not all of Mathematica's internal functions are overloaded to work with `NamedTensor`. However, given the wrappers `NamedTensorAsMatrix` and `NamedTensorFunction`, it should be easy to use them. If support for a particular function is desired, please file a feature request.

Currently, there is only marginal error handling. This means that you will know if you performed something illegal, as the output is either `$Invalid` or lots of rather unhelpful error messages appear. Check the spelling of your indices and make sure that there are no duplicate indices where they should not be! Make sure you always put the index names in quotation marks, else the prefixing will not work.

This package was created and tested with Mathematica 11.3. Based on a function analysis, the following functions are used that do not natively exist in Mathematica 10.0:

- `KeyValueMap` (since 10.1)
- `ContainsExactly` (since 10.2)
- `Nothing` (since 10.2)
- `UpTo` (since 10.3)
- `TakeList` (since 11.2)

For all these functions, the package provides trimmed-down replacements that _should_ fulfill their purpose and make it work with any version since 10.0, in which support for Associations was added (so we cannot [or do not want to] go before 10.0). However, these replacement may badly affect performance and they are untested.