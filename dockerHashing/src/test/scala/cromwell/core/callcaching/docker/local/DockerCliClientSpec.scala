package cromwell.core.callcaching.docker.local

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}

import scala.util.{Failure, Success}

class DockerCliClientSpec extends FlatSpecLike with Matchers with TableDrivenPropertyChecks {
  behavior of "DockerCliClient"

  private val lookupSuccessStdout = Seq(
    "<none>\t<none>\t<none>",
    "fauxbuntu\tlatest\tsha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
    "fauxbuntu\tmytag\tsha256:00001111222233334444555566667777888899990000aaaabbbbccccddddeeee")

  private val lookupSuccessHashes = Table(
    ("dockerCliKey", "hashValue"),
    (DockerCliKey("fauxbuntu", "latest"), "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"),
    (DockerCliKey("fauxbuntu", "mytag"), "00001111222233334444555566667777888899990000aaaabbbbccccddddeeee"))

  forAll(lookupSuccessHashes) { (dockerCliKey, hashValue) =>
    it should s"successfully lookup simulated hash for ${dockerCliKey.fullName}" in {
      val client = new TestDockerCliClient(DockerCliResult(0, lookupSuccessStdout, Nil))
      val response = client.lookupHash(dockerCliKey)
      response should be(Success(Option(s"sha256:$hashValue")))
    }
  }

  val lookupFailStderr = Seq("Error response from daemon: Bad response from Docker engine")

  it should "return not found for simulated hash fauxbuntu:unknowntag" in {
    val client = new TestDockerCliClient(DockerCliResult(0, lookupSuccessStdout, Nil))
    val dockerCliKey = DockerCliKey("fauxbuntu", "unknowntag")
    val response = client.lookupHash(dockerCliKey)
    response should be(Success(None))
  }

  it should "return expected error messages" in {
    val client = new TestDockerCliClient(DockerCliResult(1, Nil, lookupFailStderr))
    val dockerCliKey = DockerCliKey("fauxbuntu", "dockernotrunning")
    val response = client.lookupHash(dockerCliKey)
    response should be(a[Failure[_]])
    val exception = response.asInstanceOf[Failure[_]].exception
    exception.getMessage should be(
      """|Error running: docker images --digests --format {{printf "%s\t%s\t%s" .Repository .Tag .Digest}}
         |Exit code: 1
         |Error response from daemon: Bad response from Docker engine
         |""".stripMargin)
  }
}

class TestDockerCliClient(dockerCliResult: DockerCliResult) extends DockerCliClient {
  override private[local] def run(cmd: Seq[String]) = dockerCliResult
}
