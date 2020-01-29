### WDL Expression Evaluation

#### Expressions in WOM

Expressions in the WOM model all expose the following methods:

* List the names of inputs which the expression will need in order to evaluate.
* List the files which would need to be available in order to evaluate.
* Evaluate the type of value which the expression will evaluate to.
* Evaluate the expression.  

#### How WDL expressions become WomExpressions

Relating back to the [WDL parsing](../workflowParsing/wdlParsingOverview.md) process:

* WDL 1.0 engine functions become WDLOM `ExpressionElement`s during the transliteration phase
* WDLOM `ExpressionElement`s are wrapped into `WdlomWomExpression`s during the graph building phase.

It is during the graph construction phase that raw/static expression elements are mixed together with evaluation functions
to produce the final WOM expressions, and any differences between language versions are created. 

#### Where evaluation functions come from

