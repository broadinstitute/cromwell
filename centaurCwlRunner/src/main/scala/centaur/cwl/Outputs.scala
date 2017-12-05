package centaur.cwl

import cats.effect.IO
import centaur.api.CentaurCromwellClient
import cromwell.api.model.SubmittedWorkflow
import cwl.{CwlDecoder, CwlType, MyriadOutputType, File => CwlFile}
import shapeless.Inl
import spray.json.{JsObject, JsString, JsValue}

import io.circe.syntax._
import io.circe.shapes._
import io.circe.generic.auto._
import io.circe.refined._
import io.circe.literal._

object Outputs {
  import cwl.Implicits._
  import cwl.Implicits.enumerationEncoder

  def handleOutput(submittedWorkflow: SubmittedWorkflow): String = {
    import scalaz.syntax.std.map._
    val metadata: Map[String, JsValue] = CentaurCromwellClient.metadata(submittedWorkflow).get.value

    //parse the expected outputs
    Thread.sleep(100)
    metadata.get("submittedFiles.workflow") match {
      case Some(JsString(some)) =>

          val cwl = CwlDecoder.decodeTopLevelCwl(some)

          cwl.value.attempt.unsafeRunSync() match {
            case Right(Right(cwl)) =>
              val typeMap: Map[String, MyriadOutputType] = cwl.fold(CwlOutputsFold)

              CentaurCromwellClient.outputs(submittedWorkflow).get.outputs match {
                case JsObject(map) =>
                  val mungeOutputMap = map.mapKeys(key => key.substring(key.lastIndexOf(".") + 1, key.length))
                  val mungeTypeMap = typeMap.mapKeys(key => key.substring(key.lastIndexOf("#") + 1, key.length))
                  mungeOutputMap.intersectWith(mungeTypeMap)({
                    case (JsString(metadata), Inl(CwlType.File)) =>


                      val fileName = better.files.File(metadata)

                      val file = fs2.io.file.readAll[IO](fileName.path, 65536)
                      val hash = fs2.hash.sha1[IO]

                      val bytes = (file through hash).runLog.unsafeRunSync

                      val hex = bytes.map{
                        b => String.format("%02X", new java.lang.Integer(b & 0xff))
                      }.mkString.toLowerCase



                      CwlFile(
                        location = Option(fileName.name),
                        checksum = Option("sha1$" + hex),

                        size = Option(fileName.size)

                      )
                    case _ => ??? //bad

                  }).asJson.pretty(io.circe.Printer.spaces2.copy(dropNullValues = true))
                case _ => ??? //bad

              }


            case _ => println(s"couldn't parse $some"); ??? //TODO: do something bad

          }
        case _ => ???
      }
  }
}
