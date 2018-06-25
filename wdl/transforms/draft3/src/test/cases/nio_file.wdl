version 1.0

task nio_file {
  parameter_meta {
    f: { localization_optional: true }
    g: { localization_optional: true }
    h: { localization_optional: true }
  }

  input {
    File f
    File g = f
    File? h
  }

  command {
    echo ~{f} | cut -c 1-5
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Int i = 5
  }
}
