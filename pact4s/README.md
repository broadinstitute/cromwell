c# pact4s [Under construction]

pact4s is used for contract testing.

# Dependencies

```scala
  val pact4sDependencies = Seq(
    pact4sScalaTest,
    pact4sCirce,
    pact4sSpray
    http4sEmberClient,
    http4sDsl,
    http4sEmberServer,
    http4sCirce,
    circeCore,
    typelevelCat,
    scalaTest
  )

lazy val pact4s = project.in(file("pact4s"))
  .settings(pact4sSettings)
  .dependsOn(http % "test->test;compile->compile")
```

## Building and running contract tests
Clone the repo.
```
$ git clone https://github.com/broadinstitute/cromwell.git 
$ cd cromwell
```

If you are already using OpenJDK 11, run the following command. 
```
$ sbt "project pact4s" clean test  
```

Otherwise, you can run the command inside a docker container with OpenJDK 11 installed. 
This is especially useful when automating contract tests in a GitHub Action runner which does not guarantee the correct OpenJDK version.
```
docker run --rm -v $PWD:/working \
                -v jar-cache:/root/.ivy \
                -v jar-cache:/root/.ivy2 \
                -w /working \
                sbtscala/scala-sbt:openjdk-11.0.16_1.8.1_2.13.10 \
                sbt "project pact4s" clean test
```

The generated contracts can be found in the `./target/pacts` folder
- `cromwell-drshub.json`
- `cromwell-cbas.json`

