package cromwell.services.instrumentation.impl.stackdriver

import akka.actor.ActorRef
import com.google.monitoring.v3.CreateTimeSeriesRequest
import com.typesafe.config.Config
import scala.collection.JavaConverters._

class TestStackdriverInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends StackdriverInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor) {

  var metricsReceived = List[TimeSeriesRequest]()

  override def sendTimeSeriesToStackdriver(timeSeriesRequest: CreateTimeSeriesRequest) = {
    val timeSeries = timeSeriesRequest.getTimeSeries(0)
    val metric = timeSeries.getMetric

    metricsReceived = metricsReceived :+ TimeSeriesRequest(metric.getType,
      timeSeries.getPoints(0).getValue.getDoubleValue,
      timeSeries.getResource.getLabelsMap.asScala.toMap,
      metric.getLabelsMap.asScala.toMap)
  }
}


case class TimeSeriesRequest(metricPath: String,
                             metricValue: Double,
                             resourceLabels: Map[String, String],
                             metricLabels: Map[String, String])
