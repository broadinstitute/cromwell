package cloud.nio.impl.ftp

import org.apache.commons.net.ftp.FTPClient
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito

class FtpCredentialsSpec extends FlatSpec with Matchers with Mockito {

  behavior of "FtpCredentialsSpec"

  it should "login when using authenticated credentials" in {
    var loggedInWithAccount: Boolean = false
    var loggedInWithoutAccount: Boolean = false
    val client = mock[FTPClient]
    client.login(anyString, anyString).responds(_ => {
      loggedInWithoutAccount = true
      true
    })
    client.login(anyString, anyString, anyString).responds(_ => {
      loggedInWithAccount = true
      true
    })
    
    FtpAuthenticatedCredentials("user", "password", None).login(client)
    loggedInWithoutAccount shouldBe true
    loggedInWithAccount shouldBe false

    // reset
    loggedInWithoutAccount= false
    
    FtpAuthenticatedCredentials("user", "password", Option("account")).login(client)
    loggedInWithAccount shouldBe true
    loggedInWithoutAccount shouldBe false
  }

}
