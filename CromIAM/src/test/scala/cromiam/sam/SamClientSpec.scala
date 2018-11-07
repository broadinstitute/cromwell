package cromiam.sam

import akka.actor.ActorSystem
import akka.http.scaladsl.model.HttpRequest
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.stream.ActorMaterializer
import cromiam.auth.{Collection, User}
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamDenialException, SamRegisterCollectionException}
import cromiam.webservice.MockSamClient
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import org.scalatest.{AsyncFlatSpec, BeforeAndAfterAll, Matchers}

import scala.concurrent.ExecutionContextExecutor

class SamClientSpec extends AsyncFlatSpec with Matchers with BeforeAndAfterAll  {

  implicit val actorSystem: ActorSystem = ActorSystem("SamClientSpec")
  implicit val ece: ExecutionContextExecutor = actorSystem.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  val samClient = new MockSamClient()

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val authorizedUserWithCollection: User = User(WorkbenchUserId(samClient.authorizedUserCollectionStr), authorization)
  val unauthorizedUserWithNoCollection: User = User(WorkbenchUserId(samClient.unauthorizedUserCollectionStr), authorization)
  val notWhitelistedUser: User = User(WorkbenchUserId(samClient.notWhitelistedUser), authorization)

  val authorizedCollection: Collection = Collection(samClient.authorizedUserCollectionStr)
  val unauthorizedCollection: Collection = Collection(samClient.unauthorizedUserCollectionStr)
  val authorizedCollectionRequest: CollectionAuthorizationRequest = CollectionAuthorizationRequest(authorizedUserWithCollection, authorizedCollection, "add")
  val unauthorizedCollectionRequest: CollectionAuthorizationRequest = CollectionAuthorizationRequest(unauthorizedUserWithNoCollection, unauthorizedCollection, "add")

  val emptyHttpRequest: HttpRequest = HttpRequest()

  override protected def afterAll(): Unit = {
    actorSystem.terminate()
    super.afterAll()
  }


  behavior of "SamClient"
  it should "return true if user is whitelisted" in {
    samClient.isSubmitWhitelisted(authorizedUserWithCollection, emptyHttpRequest).map(v => assert(v))
  }

  it should "return false if user is not whitelisted" in {
    samClient.isSubmitWhitelisted(notWhitelistedUser, emptyHttpRequest).map(v => assert(!v))
  }

  it should "eventually return the collection(s) of user" in {
    samClient.collectionsForUser(authorizedUserWithCollection, emptyHttpRequest).map(collectionList =>
      assert(collectionList == samClient.userCollectionList)
    )
  }

  it should "fail if user doesn't have any collections" in {
    recoverToExceptionIf[Exception] {
      samClient.collectionsForUser(unauthorizedUserWithNoCollection, emptyHttpRequest)
    } map(exception =>
      assert(exception.getMessage == s"Unable to look up collections for user ${unauthorizedUserWithNoCollection.userId.value}!")
    )
  }

  it should "return true if user is authorized to perform action on collection" in {
    samClient.requestAuth(authorizedCollectionRequest, emptyHttpRequest).map(_ => succeed)
  }

  it should "throw SamDenialException if user is not authorized to perform action on collection" in {
    recoverToExceptionIf[SamDenialException] {
      samClient.requestAuth(unauthorizedCollectionRequest, emptyHttpRequest)
    } map { exception =>
      assert(exception.getMessage == "Access Denied")
    }
  }

  it should "register collection to Sam if user has authorization to create/add to collection" in {
    samClient.requestSubmission(authorizedUserWithCollection, authorizedCollection, emptyHttpRequest).map(_ => succeed)
  }

  it should "throw SamRegisterCollectionException if user doesn't have authorization to create/add to collection" in {
    recoverToExceptionIf[SamRegisterCollectionException] {
      samClient.requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest).map(_ => succeed)
    } map { exception =>
      assert(exception.getMessage == "Can't register collection with Sam. Status code: 400 Bad Request")
    }
  }
}
