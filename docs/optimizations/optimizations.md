# Optimization

These optimizations are *Cromwell-specific* functionality which can be triggered from within your workflow descriptions.

You can think of them as ways of telling Cromwell that the task or workflow has a certain property which allows Cromwell to do something clever. 
 

## Portability Warning

Note that these optimizations are *outside* of the language specifications and so not all workflow engines will respect them.
In order to maintain portability of workflows, write defensively with respect to these optimizations: 

  * Remember that a Cromwell instance might have your optimization turned off.
  * Remember that your workflow might need to run on a version of Cromwell which predates the optimization.
  * Remember that to share your WDL most widely, it will need to be able to run on engines other than Cromwell - and those engines won't necessarily respect these optimizations.
