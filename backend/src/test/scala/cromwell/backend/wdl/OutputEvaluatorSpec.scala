package cromwell.backend.wdl

import org.scalatest.{FlatSpec, Matchers}

class OutputEvaluatorSpec extends FlatSpec with Matchers {
  behavior of "OutputEvaluator"
  
  it should "return an InvalidJobOutputs if errors are found during evaluation" in {
    OutputEvaluator.evaluateOutputs()
  }
}
