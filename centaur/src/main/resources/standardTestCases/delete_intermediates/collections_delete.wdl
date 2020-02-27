version 1.0

struct delete_and_keep {
    File delete
    File keep
}

struct only_keep {
    File keep
}

task collections {
    input {
    }

    command <<<
        echo "I wanna be array deleted" > delete_array.txt
        echo "I wanna be pair deleted" > delete_pair.txt
        echo "I wanna be map deleted" > delete_map.txt
        echo "I wanna be object deleted" > delete_object.txt
        echo "I wanna be struct deleted" > delete_struct.txt
        echo "I wanna be array of array deleted" > delete_array_of_array.txt
        echo "I think I'll go for an array walk" > keep_array.txt
        echo "I think I'll go for a pair walk" > keep_pair.txt
        echo "I think I'll go for a map walk" > keep_map.txt
        echo "I think I'll go for an object walk" > keep_object.txt
        echo "I think I'll go for a struct walk" > keep_struct.txt
        echo "I think I'll go for an array of array walk" > keep_array_of_array.txt
    >>>

    runtime {
        docker: "ubuntu"
    }

    output {
        Array[File] out_array = ["delete_array.txt", "keep_array.txt"]
        Pair[File, File] out_pair = ("delete_pair.txt", "keep_pair.txt")
        Map[String, File] out_map = { "delete": "delete_map.txt", "keep": "keep_map.txt" }
        File out_object_delete = "delete_object.txt"
        File out_object_keep = "keep_object.txt"
        Object out_object = object { delete: out_object_delete, keep: out_object_keep }
        delete_and_keep out_struct = object { delete: "delete_struct.txt", keep: "keep_struct.txt" }
        Array[Array[File]] out_array_of_array = [["delete_array_of_array.txt"], ["keep_array_of_array.txt"]]
    }
}

workflow collections_delete {
    call collections {
        input:
    }

    output {
        Array[File] keep_array = [collections.out_array[1]]
        Pair[String, File] keep_pair = ("keep", collections.out_pair.right)
        Map[String, File] keep_map = { "keep": collections.out_map["keep"] }
        Object keep_object = object { keep: collections.out_object.keep }
        only_keep keep_struct = object { keep: collections.out_struct.keep }
        Array[Array[File]] keep_array_of_array = [collections.out_array_of_array[1]]
    }
}
