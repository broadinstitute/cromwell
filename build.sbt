import Settings._

lazy val wom = (project in file("wom"))
  .settings(womSettings: _*)

lazy val wdl = (project in file("wdl"))
  .settings(wdlSettings: _*)
  .dependsOn(wom)

lazy val cwl = (project in file("cwl"))
  .settings(cwlSettings: _*)
  .dependsOn(wom)

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .aggregate(wom)
  .aggregate(wdl)
  .aggregate(cwl)
  .dependsOn(wdl)
  .dependsOn(cwl)
