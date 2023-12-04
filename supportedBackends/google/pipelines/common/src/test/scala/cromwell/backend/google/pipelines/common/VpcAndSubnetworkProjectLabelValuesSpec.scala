package cromwell.backend.google.pipelines.common

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class VpcAndSubnetworkProjectLabelValuesSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "VpcAndSubnetworkProjectLabelValues"

  private val myProjectId = "my-project"

  private val labelsTests = Table(
    ("description", "network", "subnetOption", "networkName", "subnetNameOption"),
    ("a network with a slash", "slash/net", None, "slash/net", None),
    ("a network without a slash", "net", None, "projects/my-project/global/networks/net/", None),
    ("a subnet with a slash", "slashed/net", Option("slashed/sub"), "slashed/net", Option("slashed/sub")),
    ("a subnet without a slash", "slashed/net", Option("sub"), "slashed/net", Option("sub")),
    (
      "a network with a project token",
      s"slashed/$${projectId}/net",
      None,
      "slashed/my-project/net",
      None
    ),
    (
      "a subnet with a project token",
      "slashed/net",
      Option(s"slashed/$${projectId}/sub"),
      "slashed/net",
      Option("slashed/my-project/sub")
    )
  )

  forAll(labelsTests) { (description, network, subnetOption, networkName, subnetNameOption) =>
    it should description in {
      val labels = VpcAndSubnetworkProjectLabelValues(network, subnetOption)
      labels.networkName(myProjectId) should be(networkName)
      labels.subnetNameOption(myProjectId) should be(subnetNameOption)
    }
  }
}
