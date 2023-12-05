package cromiam.auth

import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import spray.json._

import scala.util.{Success, Try}

final case class Collection(name: String) extends AnyVal

object Collection {

  /**
    * Parses a raw JSON string to make sure it fits the standard pattern (see below) for labels,
    * performs some CromIAM-specific checking to ensure the user isn't attempting to manipulate the
    * protected collection name label, and then returns.
    *
    * re the Option[String], I haven't been able to figure out the Akka Http fu to get this to work as just
    * a String and the calling point in SubmissionSupport. Will loop around in a future refactor
    */
  def validateLabels(labelsJson: Option[String]): Directive1[Option[Map[String, JsValue]]] = {

    val labels = labelsJson map { l =>
      Try(l.parseJson) match {
        case Success(JsObject(json)) if json.keySet.contains(CollectionLabelName) =>
          throw new LabelContainsCollectionException
        case Success(JsObject(json)) => json
        case _ => throw InvalidLabelsException(l)
      }
    }

    provide(labels)
  }

  val CollectionLabelName = "caas-collection-name"
  val LabelsKey = "labels"

  // LabelContainsCollectionException is a class because of ScalaTest, some of the constructs don't play well w/ case objects
  final class LabelContainsCollectionException
      extends Exception(s"Submitted labels contain the key $CollectionLabelName, which is not allowed\n")
  final case class InvalidLabelsException(labels: String)
      extends Exception(s"Labels must be a valid JSON object, received: $labels\n")

  /**
    * Returns the default collection for a user.
    */
  def forUser(user: User): Collection =
    Collection(user.userId.value)

  implicit val collectionJsonReader = new JsonReader[Collection] {
    import spray.json.DefaultJsonProtocol._
    override def read(json: JsValue): Collection = Collection(json.convertTo[String])
  }
}
