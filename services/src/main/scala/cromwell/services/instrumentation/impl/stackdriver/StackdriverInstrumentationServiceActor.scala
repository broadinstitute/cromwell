package cromwell.services.instrumentation.impl.stackdriver

import java.util

import com.google.api.gax.core.FixedCredentialsProvider
import com.google.api.{Metric, MonitoredResource}
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.monitoring.v3.{MetricServiceClient, MetricServiceSettings}
import com.google.monitoring.v3._
import com.google.protobuf.util.Timestamps

import scala.collection.JavaConverters._


class StackdriverInstrumentationServiceActor {

  def pushMetrics() = {
    val projectId = "broad-dsde-cromwell-dev"

    // Instantiates a client// Instantiates a client

    val appDefaultCred = GoogleCredentials.getApplicationDefault

    val metricServiceSettings = MetricServiceSettings.newBuilder.setCredentialsProvider(FixedCredentialsProvider.create(appDefaultCred)).build
    val metricServiceClient = MetricServiceClient.create(metricServiceSettings)

    // Prepares an individual data point
    val interval = TimeInterval.newBuilder.setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
    val value = TypedValue.newBuilder.setDoubleValue(123.45).build
    val point: Point = Point.newBuilder.setInterval(interval).setValue(value).build

    val pointList = List[Point](point)

    val name = ProjectName.of(projectId)

    // Prepares the metric descriptor
    val metricLabels = Map[String, String](("store_id", "Pittsburg"))
    val metric = Metric.newBuilder.setType("custom.googleapis.com/stores/daily_sales").putAllLabels(metricLabels.asJava).build

    // Prepares the monitored resource descriptor
    val resourceLabels = Map[String, String](("project_id", projectId))
    val resource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build

    // Prepares the time series request
    val timeSeries = TimeSeries.newBuilder.setMetric(metric).setResource(resource).addAllPoints(pointList.asJava).build
    val timeSeriesList = new util.ArrayList[TimeSeries]
    timeSeriesList.add(timeSeries)

    val request = CreateTimeSeriesRequest.newBuilder.setName(name.toString).addAllTimeSeries(timeSeriesList).build

    // Writes time series data
    metricServiceClient.createTimeSeries(request)

    System.out.printf("Done writing time series data.%n")

    metricServiceClient.close()
  }

}
