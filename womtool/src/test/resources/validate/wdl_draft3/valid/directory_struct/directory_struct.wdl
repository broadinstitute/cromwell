##
## WARNING!
##
## If you've stumbled across the file at random, please DO NOT use this as an example of a good struct name!
## I'm testing that 'Directory' would be a valid struct name in WDL 1.0 to make sure Cromwell can handle an edge case...
##
## But beware!
## In future versions, the 'Directory' is instead a primitive type like File or String so I wouldn't recommend
## using a struct named Directory like this, or else you'll have serious headaches when you try to upgrade...!
##

version 1.0

## TL;DR: DO NOT copy this struct name, it won't work past WDL 1.0...
struct Directory {
  Int i
}

workflow use_directory {
  call make_directory
  call read_from_directory { input: d = make_directory.d }

  output {
    Int out = read_from_directory.contents
  }
}

task make_directory {
  command {
    echo 1 > int.txt
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Directory d = object { i: read_int(stdout()) }
  }
}

task read_from_directory {
  input {
    Directory d
  }
  command {
    echo ~{d.i}
  }
  runtime {
    docker: "ubuntu:latest"
  }
  output {
    Int contents = read_int(stdout())
  }
}
