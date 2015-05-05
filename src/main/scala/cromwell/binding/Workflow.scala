package cromwell.binding

import cromwell.binding.types.WdlType

/**
 * Represents a `workflow` definition in WDL which currently
 * only supports a set of `call` declarations and a name for
 * the workflow
 *
 * @param name The name of the workflow
 * @param calls The set of `call` declarations
 */
case class Workflow(name: String, calls: Set[Call]) extends Scope {
    /** Parent node for this workflow.  Since we do not support nested
      * workflows currently, this is always `None`
      */
    val parent: Option[Scope] = None

    /**
     * All inputs for this workflow and their associated types.
     *
     * @return a Map[FullyQualifiedName, WdlType] representing the
     *         inputs that the user needs to provide to this workflow
     */
    def inputs: Map[FullyQualifiedName, WdlType] = {
        val inputs = for(call <- calls; input <- call.unsatisfiedInputs) yield (s"${call.fullyQualifiedName}.${input._1}", input._2)
        inputs.toMap
    }

    /**
     * All outputs for this workflow and their associated types
     *
     * @return a Map[FullyQualifiedName, WdlType] representing the union
     *         of all outputs from all `call`s within this workflow
     */
    def outputs: Map[FullyQualifiedName, WdlType] = {
        val outputs = for(call <- calls; output <- call.task.outputs) yield (s"${call.fullyQualifiedName}.${output.name}", output.wdlType)
        outputs.toMap
    }
}
