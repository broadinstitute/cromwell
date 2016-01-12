import sbt._
import Keys._

object CromwellBuild extends Build {
    lazy val root = Project(id = "Cromwell",
                            base = file(".")) aggregate(CromwellBackend)

    lazy val CromwellBackend = Project(id = "Cromwell-backend",
                           base = file("cromwell-backend"))
}
