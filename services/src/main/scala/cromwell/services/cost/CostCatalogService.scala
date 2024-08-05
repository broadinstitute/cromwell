package cromwell.services.cost

import akka.actor.ActorRef
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport
import com.google.api.client.json.gson.GsonFactory
import com.typesafe.scalalogging.LazyLogging
import cromwell.services.instrumentation.CromwellInstrumentation
import com.google.api.services.cloudbilling.Cloudbilling
import com.google.api.services.cloudbilling.model.Sku
import scala.jdk.CollectionConverters._

class CostCatalogService(serviceRegistry: ActorRef) extends LazyLogging
  with CromwellInstrumentation {

  override def serviceRegistryActor: ActorRef = serviceRegistry
  val JsonFactory = GsonFactory.getDefaultInstance
  val HttpTransport = GoogleNetHttpTransport.newTrustedTransport
  val billingClient = constructBillingClient

  def fetchPublicCostCatalog: List[Sku] = {
    val res = billingClient.services().skus().list("parent").execute()
    res.getSkus.asScala.toList
  }

  private def constructBillingClient: Cloudbilling = {
    new Cloudbilling.Builder(
      HttpTransport,
      JsonFactory,
      null
    )
      .setApplicationName("Cromwell Cloud Billing Client")
      .build()
  }
}
