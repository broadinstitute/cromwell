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
  
  # Unfortunately we can't use the same trick as in draft3_nio_file_papi1 because in papi 2 the instance metadata doesn't seem to contain
  # the PAPI operation ID. So instead simply check that the input files are not there. This is not as neat as the PAPI1 version
  # but should check the right thing as long as the files are localized (in general, but not here) somewhere in the same directory as where the command
  # runs. If that changes then the directory where "find" searches should be updated 
  command {

    # Check that the paths provided were cloud:
    echo ~{f} | cut -c 1-5
    echo ~{g} | cut -c 1-5
    echo ~{h} | cut -c 1-5
    echo ~{x[0]} | cut -c 1-5
    echo ~{y.val1} | cut -c 1-5

    # Check that the NIO files were not localized
    touch errors.txt
    find . -name f >> errors.txt
    find . -name h >> errors.txt
    find . -name x >> errors.txt
    find . -name y >> errors.txt
  }
  runtime {
    docker: "google/cloud-sdk:slim"
    zones: ["us-central1-c"]
  }
  output {
    Array[String] result = read_lines(stdout())
    String errors = read_string("errors.txt")
    Boolean done = true
  }
}
