package cromwell.engine.workflow.lifecycle.execution.job.preparation

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import cromwell.core.CallKey
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import common.mock.MockSugar
import shapeless._
import wom.callable.Callable.RequiredInputDefinition
import wom.expression.{NoIoFunctionSet, WomExpression}
import wom.graph.CallNode.{InputDefinitionMappings, InputDefinitionPointer}
import wom.graph.CommandCallNode
import wom.types.WomSingleFileType
import wom.values.WomString

class CallPreparationSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockSugar {
  it should "disallow empty Strings being input as Files" in {
    val callKey = mock[CallKey]

    val inputExpressionPointer: InputDefinitionPointer =
      Coproduct[InputDefinitionPointer](WomString("").asWomExpression: WomExpression)
    val inputs: InputDefinitionMappings = List(
      (RequiredInputDefinition("inputVal", WomSingleFileType), inputExpressionPointer)
    )
    val node = CommandCallNode(null, null, null, inputs, Set.empty, null, None)

    callKey.node returns node

    val valueStore = ValueStore.empty

    CallPreparation.resolveAndEvaluateInputs(callKey, NoIoFunctionSet, valueStore) shouldBeInvalid
      """Failed to evaluate input 'inputVal' (reason 1 of 1): Cannot coerce the empty String value "" into a File."""
  }
}
