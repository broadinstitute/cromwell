workflow read_dollared_string {
  call dollars_in_strings
}


task dollars_in_strings {
  String dollar = "$"
  command <<<
    cat > foo.txt << 'EOF'
      oops ${dollar}{BLAH}
    EOF
  >>>

  output {
    File x = "foo.txt"
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
