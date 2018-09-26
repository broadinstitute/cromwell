package wes2cromwell

import akka.http.scaladsl.model.{HttpEntity, MediaTypes, MessageEntity, Multipart}
import cromiam.auth.Collection.LabelsKey
import cromiam.webservice.SubmissionSupport._

final case class WesSubmission(workflowParams: String,
                               workflowType: String,
                               workflowTypeVersion: String,
                               tags: Option[String],
                               workflowEngineParameters: Option[String],
                               workflowUrl: String,
                               workflowAttachment: Iterable[String]
                              ) {
  val entity: MessageEntity = {
    /*
      FIXME:

      Super oversimplification going on here as Cromwell's API is expected to normalize w/ WES a bit over the course
      of the quarter and it's not worth doing it here and throwing it out later (there's close to a 0% chance this will
      be used between now and then)

      At the moment, just taking the head of workflowAttachment (if it exists) and dropping that into the workflow source,
      and we're ignoring the possible existence of a zip bundle. Eventually the way this will work is that there'll be
      an Iterable of workflow files and workflowUrl will point to which one is the source and the rest goes into the bundle.

      NB: I think we've already lost too much information for the above to happen as there's an optional use of
      Content-Disposition headers on each of these files which can be used to describe directory structure and such
      for relative import resolution
     */
    val sourcePart = workflowAttachment.headOption map { a => Multipart.FormData.BodyPart(WorkflowSourceKey, HttpEntity(MediaTypes.`application/json`, a)) }

    val urlPart = Multipart.FormData.BodyPart(WorkflowUrlKey, HttpEntity(MediaTypes.`application/json`, workflowUrl))

    val typePart = Multipart.FormData.BodyPart(WorkflowTypeKey, HttpEntity(MediaTypes.`application/json`, workflowType))
    val typeVersionPart = Multipart.FormData.BodyPart(WorkflowTypeVersionKey, HttpEntity(MediaTypes.`application/json`, workflowTypeVersion))
    val inputsPart = Multipart.FormData.BodyPart(WorkflowInputsKey, HttpEntity(MediaTypes.`application/json`, workflowParams))
    val optionsPart = workflowEngineParameters map { o => Multipart.FormData.BodyPart(WorkflowOptionsKey, HttpEntity(MediaTypes.`application/json`, o)) }
    val labelsPart = tags map { t => Multipart.FormData.BodyPart(LabelsKey, HttpEntity(MediaTypes.`application/json`, t)) }

    val parts = List(sourcePart, Option(urlPart), Option(typePart), Option(typeVersionPart), Option(inputsPart), optionsPart, labelsPart).flatten

    Multipart.FormData(parts: _*).toEntity()
  }
}
