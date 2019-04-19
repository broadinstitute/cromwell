//package cromwell.docker.registryv2.flows.aliyuncr
//
//import akka.stream.ActorMaterializer
//
//import AliyunCr._
//import cromwell.docker.DockerImageIdentifierWithoutHash
//import cromwell.docker.registryv2.flows.aliyuncr.AliyunCrAbstractFlow.HttpDockerFlow
//import mouse.all._
//
//import scala.concurrent.ExecutionContext
//
//class AliyunCrFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer)
//  extends AliyunCrAbstractFlow(httpClientFlow)(ec, materializer) {
//
//
//  override def accepts(dockerImageIdentifierWithoutHash: DockerImageIdentifierWithoutHash): Boolean =
//    dockerImageIdentifierWithoutHash.host |> isValidAliyunCrHost
//}
