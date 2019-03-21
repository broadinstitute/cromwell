export runner="DataflowRunner"
#export runner="DirectRunner"
sbt "runMain example.WordCount --project=broad-dsde-cromwell-perf --runner=$runner --zone=us-central1-a"
