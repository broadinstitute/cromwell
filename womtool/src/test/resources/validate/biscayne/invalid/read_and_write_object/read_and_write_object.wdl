version development

workflow read_and_write_object {
  File foo = write_object( { "a": "b" } )
  File foo2 = write_objects( [ { "a": "b" } ]  )

  Map[String, String] foo_back = read_object(foo)
  Map[String, String] foo2_back = read_objects(foo2)

  output {
    File foo_out = foo
  }
}
