package cromwell.engine.backend.jes

import cromwell.binding.values.{WdlFile, WdlString}
import cromwell.binding.{CallInputs, Call}
import org.scalatest.{Matchers, FlatSpec}
import org.scalatest.mock.MockitoSugar
import org.mockito._

class JesBackendSpec extends FlatSpec with Matchers with MockitoSugar {

  "adjustInputPaths" should "map GCS paths and *only* GCS paths to local" in {
    val ignoredCall = mock[Call]
    val stringKey = "abc"
    val stringVal = WdlString("abc")
    val localFileKey = "lf"
    val localFileVal = WdlFile("/blah/abc")
    val gcsFileKey = "gcsf"
    val gcsFileVal = WdlFile("gs://blah/abc")


    val inputs: CallInputs = Map(
      (stringKey, stringVal),
      (localFileKey, localFileVal),
      (gcsFileKey, gcsFileVal)
    )

    val mappedInputs: CallInputs  = new JesBackend().adjustInputPaths(ignoredCall, inputs)

    mappedInputs.get(stringKey).get match {
      case WdlString(v) => assert(v.equalsIgnoreCase(stringVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(localFileKey).get match {
      case WdlFile(v) => assert(v.equalsIgnoreCase(localFileVal.value))
      case _ => fail("test setup error")
    }

    mappedInputs.get(gcsFileKey).get match {
      case WdlFile(v) => assert(v.equalsIgnoreCase("/cromwell_root/blah/abc"))
      case _ => fail("test setup error")
    }
  }
}
