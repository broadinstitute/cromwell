workflow read_dollared_strings {

  call dollars_in_strings

  String dollar = "$"

  output {
    String s1 = "${dollar}{BLAH}"
    String s2 = s1

    String s3 = dollars_in_strings.s3
  }
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
    String s3 = read_string(x)
  }

  runtime {
    docker: "ubuntu:latest"
  }
}
