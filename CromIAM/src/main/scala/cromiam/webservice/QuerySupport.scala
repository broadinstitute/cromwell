package cromiam.webservice

import akka.event.LoggingAdapter
import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import akka.stream.ActorMaterializer
import cats.data.NonEmptyList
import cromiam.auth.Collection.CollectionLabelName
import cromiam.auth.{Collection, User}
import cromiam.cromwell.CromwellClient
import cromiam.sam.SamClient
import cromiam.webservice.QuerySupport._
import cromwell.api.model._

import scala.concurrent.ExecutionContextExecutor
import scala.util.{Failure, Success, Try}

trait QuerySupport extends RequestSupport {
  val cromwellClient: CromwellClient
  val samClient: SamClient

  val log: LoggingAdapter

  implicit def executor: ExecutionContextExecutor
  implicit val materializer: ActorMaterializer

  def queryGetRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    get {
      preprocessQuery { (user, collections, request) =>
        processLabelsForGetQuery(user, collections) { uri =>
          val requestToForward = HttpRequest(HttpMethods.GET, uri, request.headers)
          complete {
            cromwellClient.forwardToCromwell(requestToForward).asHttpResponse
          }
        }
      }
    }
  }

  def queryPostRoute: Route = path("api" / "workflows" / Segment / "query") { _ =>
    post {
      preprocessQuery { (user, collections, request) =>
        processLabelsForPostQuery(user, collections) { entity =>
          complete { cromwellClient.forwardToCromwell(request.withEntity(entity)).asHttpResponse }
        }
      }
    }
  }

  /**
    * Handles shared stuff between the POST & GET Query requests. Pulls out the user information, logs the request,
    * retrieves the collections for the user, grabs the underlying HttpRequest and forwards it on to the specific
    * directive
    */
  private def preprocessQuery: Directive[(User, List[Collection], HttpRequest)] = {
    extractUserAndRequest tflatMap { case (user, cromIamRequest) =>
      log.info("Received query " + cromIamRequest.method.value + " request for user " + user.userId)

      onComplete(samClient.collectionsForUser(user, cromIamRequest).value.unsafeToFuture()) flatMap {
        case Success(Left(httpResponse)) =>
          complete(httpResponse)
        case Success(Right(collections)) =>
          toStrictEntity(Timeout) tflatMap { _ =>
            extractStrictRequest flatMap { request =>
              tprovide((user, collections, request))
            }
          }
        case Failure(e) =>
          throw new RuntimeException(s"Unable to look up collections for user ${user.userId}: ${e.getMessage}", e)
      }
    }
  }

  /**
    * Will verify that none of the GET query parameters are specifying the collection label, and then tack
    * on query parameters for the user's collections on to the query URI
    */
  private def processLabelsForGetQuery(user: User, collections: List[Collection]): Directive1[Uri] = {
    extractUri flatMap { uri =>
      val query = uri.query()

      val labelOrs: Seq[String] = query collect {
        case (key, value) if key.equalsIgnoreCase(LabelOrKey) => value
      }
      // DO NOT REMOVE THE NEXT LINE WITHOUT READING THE SCALADOC ON ensureNoLabelOrs
      ensureNoLabelOrs(user, labelOrs)

      val newQueryBuilder = query.newBuilder
      newQueryBuilder ++= query

      val collectionLabels = userCollectionLabels(user, collections)
      collectionLabels.toList foreach { collectionLabel => newQueryBuilder += (LabelOrKey -> collectionLabel) }

      provide(uri.withQuery(newQueryBuilder.result()))
    }
  }

  /**
    * Will verify that none of the POSTed query parameters are specifying the collection label, and then tack
    * on query parameters w/ a label for each of the user's collections. All of the query parameters are rebuilt
    * into a HttpEntity
    */
  private def processLabelsForPostQuery(user: User, collections: List[Collection]): Directive1[HttpEntity.Strict] = {
    import spray.json._
    import DefaultJsonProtocol._

    entity(as[String]) flatMap { paramString =>
      val params = Try(JsonParser(paramString).asInstanceOf[JsArray].elements map { _.asJsObject })
      params match {
        case Success(originalParams) =>
          val labelOrs: Vector[String] = (
            originalParams collect {
              case jsObject if jsObject.fields.keySet.exists(key => key.equalsIgnoreCase(LabelOrKey)) =>
                jsObject.fields.values.map(_.convertTo[String])
            }
            ).flatten
          // DO NOT REMOVE THE NEXT LINE WITHOUT READING THE SCALADOC ON ensureNoLabelOrs
          ensureNoLabelOrs(user, labelOrs)

          val collectionLabels = userCollectionLabels(user, collections)
          val collectionLabelParams = collectionLabels map { collectionLabel =>
            JsObject(LabelOrKey -> JsString(collectionLabel))
          }

          val newParams = JsArray(originalParams ++ collectionLabelParams.toList)
          provide(HttpEntity(ContentTypes.`application/json`, newParams.compactPrint))
        case Failure(e) => throw InvalidQueryException(e)
      }
    }
  }

  /**
    * Ensures no provided Label Ors are provided.
    *
    * Without this limitation, one could externally query ANY workflow as long as the collection name or user ID were
    * known.
    *
    * An example query to caas that would return more than expected, assuming "labelor" _were_ allowed:
    * https://caas/query?labelor=Some_Label_From_A_Hidden_Workflow:Some_Value
    *
    * Even if we built up caas' internal LabelOr queries, any workflow with "Some_Label_From_A_Hidden_Workflow"
    * would be returned.
    *
    * The permafix for this requires cromwell to support something like RQL, RSQL, or similar. Then caas could take
    * the original query, wrap it in a paren/context/group, then AND it to return only workflows from the user's
    * collections.
    *
    * Ex: (< original query here >) AND (collection IN (users_collection_a, users_collection_b, etc))
    *
    * - https://github.com/persvr/rql#rql-rules
    * - https://github.com/jirutka/rsql-parser#grammar-and-semantic
    */
  protected[this] def ensureNoLabelOrs(user: User, labelOrs: Iterable[String]): Unit = {
    labelOrs.toList match {
      case Nil => ()
      case head :: tail => throw new LabelContainsOrException(user, NonEmptyList(head, tail))
    }
  }

  /**
    * Returns the user's collections as a set of labels
    */
  protected[this] def userCollectionLabels(user: User, collections: List[Collection]): NonEmptyList[String] = {
    val userCollections = NonEmptyList(Collection.forUser(user), collections)
    userCollections.map(c => s"$CollectionLabelName:${c.name}")
  }
}

object QuerySupport {
  final case class InvalidQueryException(e: Throwable) extends
    Exception(s"Invalid JSON in query POST body: ${e.getMessage}", e)

  final class LabelContainsOrException(val user: User, val labelOrs: NonEmptyList[String]) extends
    Exception(s"User ${user.userId} submitted a labels query containing an OR which CromIAM is blocking: " +
      labelOrs.toList.mkString("LABELS CONTAIN '", "' OR LABELS CONTAIN '", "'"))

  val LabelAndKey = "label"
  val LabelOrKey = "labelor"
}
