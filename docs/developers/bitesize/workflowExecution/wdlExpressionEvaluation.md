### WDL Expression Evaluation

#### Expressions in WOM

Expressions in WOM expose the following methods:

* List the names of inputs which the expression will need in order to evaluate.
* List the files which would need to be available in order to evaluate.
* Evaluate the type of value which the expression will evaluate to.
* Evaluate the expression.  

#### How WDL expressions become WomExpressions

Relating back to the [WDL parsing](../workflowParsing/wdlParsingOverview.md) process:

* WDL 1.0 engine functions become contextless WDLOM `ExpressionElement`s during the transliteration phase
* Reference resolution between WDLOM elements occurs during the linking phase. 
* WDLOM `ExpressionElement`s are wrapped into `WdlomWomExpression`s during the graph building phase.

During the graph construction phase static expression elements are mixed together with evaluation functions
to produce the final WOM expressions, and any differences between language versions are baked in. 

#### Where evaluation functions are coded

The four evaluators types described in "Expressions in WOM" are coded in four package objects for each version. They live in:

* WDL `version 1.0`: [`wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/linking/expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/linking/expression)
* WDL `version development`: [`wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/linking/expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/biscayne/src/test/scala/wdl/transforms/biscayne/linking/expression)

These evaluators are really just long pattern matches from various WDLOM element types to the relevant evaluator for that element.

Because the evaluations are largely the same for now between WDL versions, you'll notice that the imported evaluators mostly come from files 
in [`wdl.transforms.base.linking.expression`](https://github.com/broadinstitute/cromwell/tree/develop/wdl/transforms/new-base/src/main/scala/wdl/transforms/base/linking/expression). However, in WDL development, there are a number of imported functions from biscayne-specific directories.
If future WDL versions diverge more starkly from the 1.0 base, it is likely that the imports into these pattern match evaluators will come from a
wider variety of origins.

#### How Evaluation Functions Build WOM Expressions

The top-level Graph construction functions for WDL 1.0 and the development version are:

* WDL `version 1.0`: [`wdl.draft3.transforms.wdlom2wom.workflowDefinitionElementToWomWorkflowDefinition`](https://github.com/broadinstitute/cromwell/blob/develop/wdl/transforms/draft3/src/main/scala/wdl/draft3/transforms/wdlom2wom/package.scala)
* WDL `version development`: [`wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/wdlom2wom/package`](https://github.com/broadinstitute/cromwell/blob/develop/wdl/transforms/biscayne/src/main/scala/wdl/transforms/biscayne/wdlom2wom/package.scala)

In each case:

* The package objects import the evaluation functions (from "Where evaluation functions are coded")
    * These imports fill in the implicit parameters required by `WorkflowDefinitionElementToWomWorkflowDefinition.convert`.
    * The imported conversion functions are used by the `convert` function every time it needs to make a WomExpression from WDLOM.
* The package objects create a language-specific conversion, `workflowDefinitionElementToWomWorkflowDefinition`.
    * This is used in the graph construction phase mentioned in "How WDL expressions become WomExpressions".
* The language-specific conversion is used by the appropriate LanguageFactory when it gets asked to construct WOM from WDL. 
