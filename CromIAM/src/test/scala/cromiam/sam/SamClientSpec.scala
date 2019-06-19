package cromiam.sam

import akka.actor.ActorSystem
import akka.http.scaladsl.model.headers.{Authorization, OAuth2BearerToken}
import akka.http.scaladsl.model.{HttpEntity, HttpRequest, HttpResponse, StatusCodes}
import akka.stream.ActorMaterializer
import cats.Monad
import cromiam.auth.{Collection, User}
import cromiam.sam.SamClient.{CollectionAuthorizationRequest, SamDenialException, SamRegisterCollectionException}
import cromiam.webservice.MockSamClient._
import cromiam.webservice._
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import cromwell.api.model._
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId
import org.scalatest.{AsyncFlatSpec, BeforeAndAfterAll, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContextExecutor

class SamClientSpec extends AsyncFlatSpec with Matchers with BeforeAndAfterAll with Mockito {

  implicit val actorSystem: ActorSystem = ActorSystem("SamClientSpec")
  implicit val ece: ExecutionContextExecutor = actorSystem.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  private val expectedErrorResponse =
    HttpResponse(StatusCodes.InternalServerError, entity = HttpEntity("expected error"))

  val authorization = Authorization(OAuth2BearerToken("my-token"))
  val authorizedUserWithCollection = User(WorkbenchUserId(MockSamClient.AuthorizedUserCollectionStr), authorization)
  val unauthorizedUserWithNoCollection =
    User(WorkbenchUserId(MockSamClient.UnauthorizedUserCollectionStr), authorization)
  val notWhitelistedUser = User(WorkbenchUserId(MockSamClient.NotWhitelistedUser), authorization)

  val authorizedCollection = Collection(MockSamClient.AuthorizedUserCollectionStr)
  val unauthorizedCollection = Collection(MockSamClient.UnauthorizedUserCollectionStr)
  val authorizedCollectionRequest =
    CollectionAuthorizationRequest(authorizedUserWithCollection, authorizedCollection, "add")
  val unauthorizedCollectionRequest =
    CollectionAuthorizationRequest(unauthorizedUserWithNoCollection, unauthorizedCollection, "add")

  val emptyHttpRequest: HttpRequest = HttpRequest()

  override protected def afterAll(): Unit = {
    actorSystem.terminate()
    super.afterAll()
  }


  behavior of "SamClient"
  it should "return true if user is whitelisted" in {
    val samClient = new MockSamClient()
    samClient.isSubmitWhitelisted(authorizedUserWithCollection, emptyHttpRequest).map(v => assert(v))
      .asIo.unsafeToFuture()
  }

  it should "return false if user is not whitelisted" in {
    val samClient = new MockSamClient()
    samClient.isSubmitWhitelisted(notWhitelistedUser, emptyHttpRequest).map(v => assert(!v))
      .asIo.unsafeToFuture()
  }

  it should "return sam errors while checking is whitelisted" in {
    val samClient = new MockSamClient() {
      override def isSubmitWhitelistedSam(user: User, cromiamRequest: HttpRequest): FailureResponseOrT[Boolean] = {
        MockSamClient.returnResponse(expectedErrorResponse)
      }
    }
    samClient.isSubmitWhitelisted(notWhitelistedUser, emptyHttpRequest).value.unsafeToFuture() map {
      _ should be(Left(expectedErrorResponse))
    }
  }

  it should "eventually return the collection(s) of user" in {
    val samClient = new MockSamClient()
    samClient.collectionsForUser(authorizedUserWithCollection, emptyHttpRequest).map(collectionList =>
      assert(collectionList == MockSamClient.UserCollectionList)
    ).asIo.unsafeToFuture()
  }

  it should "fail if user doesn't have any collections" in {
    val samClient = new MockSamClient()
    recoverToExceptionIf[Exception] {
      samClient.collectionsForUser(unauthorizedUserWithNoCollection, emptyHttpRequest)
        .asIo.unsafeToFuture()
    } map(exception =>
      assert(exception.getMessage == s"Unable to look up collections for user ${unauthorizedUserWithNoCollection.userId.value}!")
    )
  }

  it should "return true if user is authorized to perform action on collection" in {
    val samClient = new MockSamClient()
    samClient.requestAuth(authorizedCollectionRequest, emptyHttpRequest).map(_ => succeed)
      .asIo.unsafeToFuture()
  }

  it should "throw SamDenialException if user is not authorized to perform action on collection" in {
    val samClient = new MockSamClient()
    recoverToExceptionIf[SamDenialException] {
      samClient.requestAuth(unauthorizedCollectionRequest, emptyHttpRequest)
        .asIo.unsafeToFuture()
    } map { exception =>
      assert(exception.getMessage == "Access Denied")
    }
  }

  it should "register collection to Sam if user has authorization to create/add to collection" in {
    val samClient = new MockSamClient()
    samClient.requestSubmission(authorizedUserWithCollection, authorizedCollection, emptyHttpRequest).map(_ => succeed)
      .asIo.unsafeToFuture()
  }

  it should "throw SamRegisterCollectionException if user doesn't have authorization to create/add to collection" in {
    val samClient = new MockSamClient()
    recoverToExceptionIf[SamRegisterCollectionException] {
      samClient.requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest).map(_ => succeed)
        .asIo.unsafeToFuture()
    } map { exception =>
      assert(exception.getMessage == "Can't register collection with Sam. Status code: 400 Bad Request")
    }
  }

  it should "add the user when create returns a conflict" in {
    val samClient = new BaseMockSamClient() {
      override protected def registerCreation(user: User,
                                              collection: Collection,
                                              cromiamRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
        val conflictResponse = HttpResponse(StatusCodes.Conflict, entity = HttpEntity("expected conflict"))
        returnResponse(conflictResponse)
      }

      override def requestAuth(authorizationRequest: CollectionAuthorizationRequest,
                               cromiamRequest: HttpRequest): FailureResponseOrT[Unit] = {
        Monad[FailureResponseOrT].unit
      }
    }
    samClient
      .requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest)
      .map(_ => succeed)
      .asIo
      .unsafeToFuture()
  }

  it should "fail to add the user when create returns a conflict but then adding returns an error" in {
    val samClient = new BaseMockSamClient() {
      override protected def registerCreation(user: User,
                                              collection: Collection,
                                              cromiamRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
        val conflictResponse = HttpResponse(StatusCodes.Conflict, entity = HttpEntity("expected conflict"))
        returnResponse(conflictResponse)
      }

      override def requestAuth(authorizationRequest: CollectionAuthorizationRequest,
                               cromiamRequest: HttpRequest): FailureResponseOrT[Unit] = {
        returnResponse(expectedErrorResponse)
      }
    }
    recoverToExceptionIf[UnsuccessfulRequestException] {
      samClient.requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest)
        .asIo.unsafeToFuture()
    } map { exception =>
      assert(exception.getMessage == "expected error")
    }
  }

  it should "fail to add the user when create returns an unexpected successful response" in {
    val samClient = new BaseMockSamClient() {
      override protected def registerCreation(user: User,
                                              collection: Collection,
                                              cromiamRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
        val unexpectedOkResponse = HttpResponse(StatusCodes.OK, entity = HttpEntity("elided ok message"))
        returnResponse(unexpectedOkResponse)
      }
    }
    recoverToExceptionIf[SamRegisterCollectionException] {
      samClient.requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest)
        .asIo.unsafeToFuture()
    } map { exception =>
      exception.getMessage should be("Can't register collection with Sam. Status code: 200 OK")
    }
  }

  it should "fail to add the user when create returns an unexpected failure response" in {
    val samClient = new BaseMockSamClient() {
      override protected def registerCreation(user: User,
                                              collection: Collection,
                                              cromiamRequest: HttpRequest): FailureResponseOrT[HttpResponse] = {
        val unexpectedFailureResponse = HttpResponse(StatusCodes.ImATeapot, entity = HttpEntity("elided error message"))
        returnResponse(unexpectedFailureResponse)
      }
    }
    recoverToExceptionIf[SamRegisterCollectionException] {
      samClient.requestSubmission(unauthorizedUserWithNoCollection, unauthorizedCollection, emptyHttpRequest)
        .asIo.unsafeToFuture()
    } map { exception =>
      exception.getMessage should be("Can't register collection with Sam. Status code: 418 I'm a teapot")
    }
  }
}
