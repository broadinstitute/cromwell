package statsdproxy

import java.net.{InetSocketAddress, URI}
import java.nio.file.{Paths, StandardOpenOption}
import java.util.concurrent.Executors

import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits._
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.StrictLogging
import fs2.io.file.pulls.{fromPath, writeAllToFileHandle}
import fs2.io.udp.{AsynchronousSocketGroup, Packet, Socket}
import fs2.{Pull, Stream}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.ExecutionContext

/**
  * Proxies incoming udp connection to another udp socket, and write all packets to a file.
  */
object StatsDProxy extends IOApp with StrictLogging {
  implicit private val ag = AsynchronousSocketGroup()
  private val blockingEc = ExecutionContext.fromExecutor(Executors.newSingleThreadExecutor)
  private val config = ConfigFactory.load()
  private val proxyHost = config.as[String]("proxy.host")
  private val proxyPort = config.as[Int]("proxy.port")
  private val statsdHost = config.as[String]("statsd.host")
  private val statsdPort = config.as[Int]("statsd.port")
  private val filePath = config.as[String]("output-file")
  
  private val proxyAddress = new InetSocketAddress(proxyHost, proxyPort)
  private val redirectAddress = new InetSocketAddress(statsdHost, statsdPort)

  logger.info(s"Redirecting statsd packets to ${redirectAddress.getHostString}:${redirectAddress.getPort}")
  logger.info(s"Writing statsd packets to $filePath")

  private val byteStream: Stream[IO, Byte] = for {
    serverSocket <- Stream.resource(Socket[IO](proxyAddress))
    _ = logger.info(s"Proxy listening at ${proxyAddress.getHostString}:${proxyAddress.getPort}")
    clientSocket <- Stream.resource(Socket[IO]())
    receivedPacket <- serverSocket.reads()
    // Re-wrap the packet in the redirect address, send it, and then unpack it to create a stream of bytes that can be written to a file
    pipedToOutbound <- Stream(Packet(redirectAddress, receivedPacket.bytes))
      .covary[IO]
      .through(clientSocket.writes())
      .map(_ => receivedPacket.bytes.toVector :+ '\n'.toByte)
      .flatMap(Stream.emits(_))
  } yield pipedToOutbound

  private val server: Pull[IO, Nothing, Unit] = for {
    filePull <- fromPath[IO](Paths.get(URI.create(filePath)), blockingEc, List(StandardOpenOption.CREATE, StandardOpenOption.WRITE))
    _ <- writeAllToFileHandle(byteStream, filePull.resource)
  } yield ()

  override def run(args: List[String]) = {
    server.stream.compile.drain.as(ExitCode.Success)
  }
}
