package cromwell.core.callcaching

import org.scalatest.{FlatSpec, Matchers}

class HashKeySpec extends FlatSpec with Matchers {

  "HashKey" should "produce consistent key value" in {
    val keys = Set(
      HashKey("command template"),
      HashKey("backend name"),
      HashKey("input count"),
      HashKey("output count"),
      HashKey("runtime attribute", "failOnStderr"),
      HashKey(checkForHitOrMiss = false, "runtime attribute", "cpu"),
      HashKey("runtime attribute", "continueOnReturnCode"),
      HashKey("input", "String stringInput"),
      HashKey("output", "String myOutput"),
      HashKey("runtime attribute", "docker")
    )
    
    keys map { _.key } should contain theSameElementsAs Set(
      "command template",
      "backend name",
      "input count",
      "output count",
      "runtime attribute: failOnStderr",
      "runtime attribute: cpu",
      "runtime attribute: continueOnReturnCode",
      "input: String stringInput",
      "output: String myOutput",
      "runtime attribute: docker"
    )
  }
  
}
