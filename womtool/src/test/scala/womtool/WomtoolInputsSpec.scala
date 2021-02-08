package womtool

import womtool.WomtoolJsonCommandSpec.TestDefinition

class WomtoolInputsSpec extends WomtoolJsonCommandSpec {
  override val womtoolCommand: String = "inputs"
  override val testDefinitions: Seq[TestDefinition] = Seq(
    TestDefinition("required inputs", Seq("inputs", "-o", "false"), "required.inputs.json"),
    TestDefinition("all inputs", Seq("inputs"), "all.inputs.json")
  )

  defineTests()
}
