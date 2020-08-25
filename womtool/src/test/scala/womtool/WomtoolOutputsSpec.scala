package womtool

import womtool.WomtoolJsonCommandSpec.TestDefinition

class WomtoolOutputsSpec extends WomtoolJsonCommandSpec {
  override val womtoolCommand: String = "outputs"
  override val testDefinitions: Seq[TestDefinition] = Seq(
    TestDefinition("all outputs", Seq("outputs"), "outputs.json")
  )

  defineTests()
}
