version 1.0

task use_this_name {

  input {
    File f
  }

  command {}

  meta {
    email: "skroob@spaceballs.gov"
    author: "President Skroob"
  }

  output {
    File f2 = "a"
  }
}
