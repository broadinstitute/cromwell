package cromwell.docker.registryv2.flows.gcr

import akka.actor.Scheduler
import akka.stream.ActorMaterializer
import cromwell.docker.registryv2.DockerRegistryV2AbstractFlow.HttpDockerFlow

import scala.concurrent.ExecutionContext

private class GcrFlow(httpClientFlow: HttpDockerFlow)(implicit ec: ExecutionContext, materializer: ActorMaterializer, scheduler: Scheduler) extends GcrAbstractFlow(httpClientFlow, "gcr.io")
