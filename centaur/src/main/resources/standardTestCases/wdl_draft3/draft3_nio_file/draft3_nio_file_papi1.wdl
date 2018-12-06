version 1.0

workflow draft3_nio_file {
  call mk_files { input: salt = 0 }
  call mk_files as mk_more_files { input: salt = 1 }

  call nio_file { input:
    ready = true,
    f = mk_files.f,
    h = mk_files.h,
    x = mk_files.x,
    y = mk_files.y
  }

  call nio_file as cc_nio_file { input:
    ready = nio_file.done,
    f = mk_files.f,
    h = mk_files.h,
    x = mk_files.x,
    y = mk_files.y
  }

  call nio_file as non_cc_nio_file { input:
    ready = nio_file.done,
    f = mk_more_files.f,
    h = mk_more_files.h,
    x = mk_more_files.x,
    y = mk_more_files.y
  }

  output {
    String f_path_prefix = nio_file.result[0]
    String g_path_prefix = nio_file.result[1]
    String h_path_prefix = nio_file.result[2]
    String x_path_prefix = nio_file.result[3]
    String y_path_prefix = nio_file.result[4]

    String errors = nio_file.errors
  }
}

struct FileBox {
  File val0
  File val1
}

task mk_files {

  parameter_meta {
    salt: "Shakes things up a little bit!"
  }

  input {
    Int salt
  }

  command {
    echo "f~{salt}" > f
    echo "h~{salt}" > h
    echo "x0~{salt}" > x0
    echo "x1~{salt}" > x1
    echo "y0~{salt}" > y0
    echo "y1~{salt}" > y1
  }

  runtime {
    docker: "ubuntu:latest"
  }

  output {
    File f = "f"
    File h = "h"
    Array[File] x = ["x0", "x1"]
    FileBox y = object { val0: "y0", val1: "y1" }
  }
}

task nio_file {

  meta {
    description: "Analyzes whether (a) the interpolated paths for NIO files start with 'gs://' and (b) the inputs are localized"
  }
  parameter_meta {
    ready: "Allows us to delay until a previous task is done"
    f: { localization_optional: true }
    g: { localization_optional: true }
    h: {
      description: "Only here to check that we can have fields before the 'nio'",
      localization_optional: true,
      after: "... and after it..."
    }
    x: { localization_optional: true }
    y: { localization_optional: true }
    done: "Always true, can be chained into the 'ready' of a subsequent invocation"
  }

  input {
    Boolean ready
    File f
    File g = f
    File? h
    Array[File] x
    FileBox y
  }

  command {

    # Check that the paths provided were cloud:
    echo ~{f} | cut -c 1-5
    echo ~{g} | cut -c 1-5
    echo ~{h} | cut -c 1-5
    echo ~{x[0]} | cut -c 1-5
    echo ~{y.val1} | cut -c 1-5

    # Check that the NIO files were not localized
    INSTANCE=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
    TOKEN=$(gcloud auth application-default print-access-token)
    INSTANCE_INFO=$(curl "https://www.googleapis.com/compute/v1/projects/broad-dsde-cromwell-dev/zones/us-central1-c/instances/$INSTANCE" -H "Authorization: Bearer $TOKEN" -H 'Accept: application/json')
    OPERATION_LINE=$(grep Operation <<< $INSTANCE_INFO)
    OPERATION_ID=$(sed 's/.*Operation: \([^ ]*\)".*/\1/' <<< $OPERATION_LINE)
    gcloud components install alpha
    PAPI_METADATA=$(gcloud --quiet alpha genomics operations describe operations/$OPERATION_ID)
    echo $PAPI_METADATA > papi_metadata

    touch errors.txt
    grep -q "draft3_nio_file.nio_file.f" <<< $PAPI_METADATA && echo "f was incorrectly localized" >> errors.txt
    grep -q "draft3_nio_file.nio_file.g" <<< $PAPI_METADATA && echo "g was incorrectly localized" >> errors.txt
    grep -q "draft3_nio_file.nio_file.__g" <<< $PAPI_METADATA && echo "__g was incorrectly localized" >> errors.txt
    grep -q "draft3_nio_file.nio_file.h" <<< $PAPI_METADATA && echo "h was incorrectly localized" >> errors.txt
    grep -q "draft3_nio_file.nio_file.x" <<< $PAPI_METADATA && echo "x was incorrectly localized" >> errors.txt
    grep -q "draft3_nio_file.nio_file.y" <<< $PAPI_METADATA && echo "y was incorrectly localized" >> errors.txt
  }
  runtime {
    docker: "google/cloud-sdk:slim"
    zones: ["us-central1-c"]

    # Depending on the final 'grep', the return code is probably going to be '1'... which is fine!
    continueOnReturnCode: true
  }
  output {
    Array[String] result = read_lines(stdout())
    String errors = read_string("errors.txt")
    Boolean done = true
  }
}
