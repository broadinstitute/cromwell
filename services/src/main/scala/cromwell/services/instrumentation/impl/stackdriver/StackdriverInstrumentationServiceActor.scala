package cromwell.services.instrumentation.impl.stackdriver

import java.util

import akka.actor.{Actor, ActorRef, Props}
import com.google.api.gax.core.FixedCredentialsProvider
import com.google.api.{Metric, MetricDescriptor, MonitoredResource}
import com.google.cloud.monitoring.v3.{MetricServiceClient, MetricServiceSettings}
import com.google.monitoring.v3._
import com.google.protobuf.util.Timestamps
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation._
import cromwell.services.instrumentation.impl.stackdriver.StackdriverConfig._
import cromwell.services.instrumentation.impl.stackdriver.StackdriverInstrumentationServiceActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.util.Try


class StackdriverInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with StrictLogging {
  implicit lazy val executionContext = context.dispatcher

  val stackdriverConfig = StackdriverConfig(serviceConfig, globalConfig)

  lazy val projectName: ProjectName = ProjectName.of(stackdriverConfig.googleProject)
  val credentials = stackdriverConfig.auth.credentials(List(MonitoringScope))
  lazy val metricLabelsMap = generateMetricLabels()

  var metricsMap = Map.empty[StackdriverMetric, Vector[Double]]

  // Instantiates a client
  val metricServiceSettings = MetricServiceSettings.newBuilder.setCredentialsProvider(FixedCredentialsProvider.create(credentials)).build
  final val metricServiceClient = MetricServiceClient.create(metricServiceSettings)

  // Prepares the monitored resource descriptor
  lazy val resourceLabels = Map[String, String](("project_id", stackdriverConfig.googleProject))
  lazy val monitoredResource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build

  // Give the actor time to warm up, then start sending the metrics to Stackdriver at an interval
  context.system.scheduler.schedule(InitialDelay, stackdriverConfig.flushRate, self, SendStackdriverMetricCommand)


  override def receive = {
    case SendStackdriverMetricCommand => sendMetricData()
    case InstrumentationServiceMessage(cromwellMetric) => cromwellMetric match {
      case CromwellTiming(bucket, value, _) => updateMetricMap(bucket, value.toMillis.toDouble, StackdriverGauge)
      case CromwellGauge(bucket, value) => updateMetricMap(bucket, value.toDouble, StackdriverGauge)
      case CromwellCount(bucket, value, _) => updateMetricMap(bucket, value.toDouble, StackdriverCumulative)
      case CromwellIncrement(bucket) => updateMetricMap(bucket, metricValue = 1D, metricKind = StackdriverCumulative)
    }
    case ShutdownCommand =>
      // flush out metrics (if any) before shut down
      sendMetricData()
      metricServiceClient.close()
      context stop self
  }


  private def generateMetricLabels(): Map[String, String] = {
    def labelFromConfig(op: StackdriverConfig => Option[String], key: String): Option[(String, String)] = {
      op(stackdriverConfig).map(v => (key.replace("-", "_"), v))
    }

    labelFromConfig(_.cromwellInstanceIdentifier, CromwellInstanceIdentifier).toMap ++
      labelFromConfig(_.cromwellInstanceRole, CromwellInstanceRole) ++
      labelFromConfig(_.cromwellPerfTestCase, CromwellPerfTest)
  }


  private def updateMetricMap(bucket: CromwellBucket, metricValue: Double, metricKind: StackdriverMetricKind): Unit = {
    val metricObj = StackdriverMetric(bucket.toStackdriverString, metricKind)

    if (metricsMap.contains(metricObj)) {
      val valueVector: Vector[Double] = metricsMap(metricObj) :+ metricValue
      metricsMap += metricObj -> valueVector
    }
    else metricsMap += metricObj -> Vector(metricValue)
  }


  private def sendMetricData(): Unit = {
    metricsMap.foreach { case (key, value: Vector[Double]) =>
      val dataPointVectorSum = value.sum
      val dataPoint = key.kind match {
        case StackdriverGauge => dataPointVectorSum / value.length
        case StackdriverCumulative => dataPointVectorSum
      }

      writeMetrics(key, dataPoint) recover {
        case e => logger.error(s"Failed to send metrics to Stackdriver API for metric ${key.name} with value $dataPoint.", e)
      }
    }

    metricsMap = Map.empty[StackdriverMetric, Vector[Double]]
  }


  private def writeMetrics(metricObj: StackdriverMetric, value: Double): Try[Unit] = {
    def timeInterval(metricKind: StackdriverMetricKind): TimeInterval = {
      metricKind match {
        case StackdriverGauge => TimeInterval.newBuilder.setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
        case StackdriverCumulative => TimeInterval.newBuilder.setStartTime(ActorCreationTime).setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
      }
    }

    def createTimeSeries(metricKind: StackdriverMetricKind, metric: Metric, resource: MonitoredResource, dataPointList: util.List[Point]): TimeSeries = {
      metricKind match {
        case StackdriverGauge => TimeSeries.newBuilder.setMetric(metric).setResource(resource).addAllPoints(dataPointList).build
        case StackdriverCumulative => TimeSeries.newBuilder.setMetric(metric).setResource(resource).setMetricKind(MetricDescriptor.MetricKind.CUMULATIVE).addAllPoints(dataPointList).build
      }
    }

    // Prepares an individual data point
    val interval = timeInterval(metricObj.kind)
    val pointValue = TypedValue.newBuilder().setDoubleValue(value).build()
    val dataPoint: Point = Point.newBuilder.setInterval(interval).setValue(pointValue).build
    val dataPointList: List[Point] = List[Point](dataPoint)

    // Prepares the metric descriptor
    val metric: Metric = Metric.newBuilder.setType(s"$CustomMetricDomain/${metricObj.name}").putAllLabels(metricLabelsMap.asJava).build

    // Prepares the time series request
    val timeSeries = createTimeSeries(metricObj.kind, metric, monitoredResource, dataPointList.asJava)
    val timeSeriesList = List[TimeSeries](timeSeries)

    val timeSeriesRequest = CreateTimeSeriesRequest.newBuilder.setName(projectName.toString).addAllTimeSeries(timeSeriesList.asJava).build

    // Writes time series data
    Try(sendTimeSeriesToStackdriver(timeSeriesRequest))
  }


  // This single line of code is a separate function to help with StackdriverInstrumentationActor tests
  def sendTimeSeriesToStackdriver(timeSeriesRequest: CreateTimeSeriesRequest): Unit = {
    metricServiceClient.createTimeSeries(timeSeriesRequest)
  }
}


object StackdriverInstrumentationServiceActor {
  val CromwellMetricPrefix = List("cromwell")

  val ActorCreationTime = Timestamps.fromMillis(System.currentTimeMillis())

  // Custom metrics must begin with this domain
  val CustomMetricDomain = "custom.googleapis.com"

  /**
    * Scope to write metrics to Stackdriver Monitoring API.
    * Used by the monitoring action.
    *
    * For some reason we couldn't find this scope within Google libraries
    */
  val MonitoringScope = "https://www.googleapis.com/auth/monitoring"

  val InitialDelay = 1.minute

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new StackdriverInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor))

  implicit class CromwellBucketEnhanced(val cromwellBucket: CromwellBucket) extends AnyVal {
    /**
      * Transforms a CromwellBucket to a Stackdriver path
      */
    def toStackdriverString = (CromwellMetricPrefix ++ cromwellBucket.prefix ++
        cromwellBucket.path.toList).mkString("/").replaceAll(" ", "_").replaceAll("\\[|\\]", "")
  }
}


sealed trait StackdriverMetricKind
object StackdriverGauge extends StackdriverMetricKind
object StackdriverCumulative extends StackdriverMetricKind


case class StackdriverMetric(name: String, kind: StackdriverMetricKind)


object SendStackdriverMetricCommand
