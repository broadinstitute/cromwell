package cromwell.core

import java.util.Base64

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class DockerCredentialsSpec extends AnyFlatSpec with Matchers {
  behavior of "DockerCredentialsUsernameAndPassword"

  val cases = List(
    ("foouser:foopass", "foouser", "foopass"),
    ("foouser:foo:pass", "foouser", "foo:pass")
  )

  cases foreach { case (tokenString, expectedUsername, expectedPassword) =>
    it should s"decompose a token formed from ${tokenString} into $expectedUsername and $expectedPassword" in {
      val credentials: Any = new DockerCredentials(Base64.getEncoder.encodeToString(tokenString.getBytes()), None, None)

      credentials match {
        case DockerCredentialUsernameAndPassword(u, p) => {
          u should be(expectedUsername)
          p should be(expectedPassword)
        }
        case _ => fail(s"Expected to decompose ${tokenString} into username=$expectedPassword and password=$expectedPassword")
      }
    }
  }

  it should "not match strings without colons" in {
    new DockerCredentials(Base64.getEncoder.encodeToString("no-colons".getBytes), None, None) match {
      case DockerCredentialUsernameAndPassword(huh, what) => fail(s"expected no decomposition but got $huh and $what")
      case _ => // success!
    }
  }
}
