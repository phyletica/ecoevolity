#include "catch.hpp"
#include "ecoevolity/path.hpp"

SCENARIO("dirname provides parent directory of path", "[path]") {

    GIVEN("A full file path") {
        std::string path = "/home/jamie/Desktop/file.txt";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("then a full path to the parent directory is returned") {
                REQUIRE(parent_dir == "/home/jamie/Desktop");
            }
        }
    }

    GIVEN("A full file path with leading space") {
        std::string path = "  /home/jamie/Desktop/file.txt";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("the space is preserved") {
                REQUIRE(parent_dir == "  /home/jamie/Desktop");
            }
        }
    }

    GIVEN("A full file path with trailing space") {
        std::string path = "/home/jamie/Desktop/file.txt  ";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("the path to the parent directory is returned with no space") {
                REQUIRE(parent_dir == "/home/jamie/Desktop");
            }
        }
    }

    GIVEN("A file name only") {
        std::string path = "file.txt";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("then an empty string is returned") {
                REQUIRE(parent_dir == "");
            }
        }
    }

    GIVEN("A file name only with leading and trailing space") {
        std::string path = "  file.txt  ";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("then an empty string is returned") {
                REQUIRE(parent_dir == "");
            }
        }
    }

    GIVEN("A path that ends with a separator") {
        std::string path = "../jamie/Desktop/";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("Then just the separator will be dropped") {
                REQUIRE(parent_dir == "../jamie/Desktop");
            }
        }
    }

    GIVEN("A path that ends with a separator and space") {
        std::string path = "../jamie/Desktop/    ";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("the separator and space will be dropped") {
                REQUIRE(parent_dir == "../jamie/Desktop");
            }
        }
    }

    GIVEN("An empty string") {
        std::string path = "";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("then an empty string is returned") {
                REQUIRE(parent_dir == "");
            }
        }
    }

    GIVEN("A string of only spaces") {
        std::string path = "    ";

        WHEN("dirname is called") {
            std::string parent_dir = path::dirname(path);
            
            THEN("then an empty string is returned") {
                REQUIRE(parent_dir == "");
            }
        }
    }
}

SCENARIO("basename returns the final component of a path", "[path]") {

    GIVEN("a full file path") {
        std::string path = "/home/jamie/Desktop/file.txt";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("the file name returned") {
                REQUIRE(base_name == "file.txt");
            }
        }
    }

    GIVEN("a full file path with leading white space") {
        std::string path = "   /home/jamie/Desktop/file.txt";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("the file name returned") {
                REQUIRE(base_name == "file.txt");
            }
        }
    }

    GIVEN("a full file path with trailing white space") {
        std::string path = "/home/jamie/Desktop/file.txt \n";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("the file name returned with space preserved") {
                REQUIRE(base_name == "file.txt \n");
            }
        }
    }

    GIVEN("A file name only") {
        std::string path = "file.txt";

        WHEN("baseame is called") {
            std::string base_name = path::basename(path);
            
            THEN("the file name is returned") {
                REQUIRE(base_name == "file.txt");
            }
        }
    }

    GIVEN("A file name only with leading and trailing space") {
        std::string path = "  file.txt  ";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("the file is returned with spaces preserved") {
                REQUIRE(base_name == "  file.txt  ");
            }
        }
    }

    GIVEN("A path that ends with a separator") {
        std::string path = "../jamie/Desktop/";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("an empty string is returned") {
                REQUIRE(base_name == "");
            }
        }
    }

    GIVEN("A path that ends with a separator and space") {
        std::string path = "../jamie/Desktop/    ";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("the space is returned") {
                REQUIRE(base_name == "    ");
            }
        }
    }

    GIVEN("An empty string") {
        std::string path = "";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("an empty string is returned") {
                REQUIRE(base_name == "");
            }
        }
    }

    GIVEN("A string of only spaces") {
        std::string path = "    ";

        WHEN("basename is called") {
            std::string base_name = path::basename(path);
            
            THEN("spaces are returned") {
                REQUIRE(base_name == "    ");
            }
        }
    }
}

SCENARIO("isabs tests whether a path is absolute", "[path]") {

    GIVEN("A full file path") {
        std::string path = "/home/jamie/Desktop/file.txt";

        WHEN("isabs is called") {
            bool is_abs = path::isabs(path);
            
            THEN("returns true") {
                REQUIRE(is_abs == true);
            }
        }
    }

    GIVEN("A full file path with leading space") {
        std::string path = " /home/jamie/Desktop/file.txt";

        WHEN("isabs is called") {
            bool is_abs = path::isabs(path);
            
            THEN("returns false") {
                REQUIRE(is_abs == false);
            }
        }
    }

    GIVEN("A path that does not being with a separator") {
        std::string path = "home/jamie/Desktop/file.txt";

        WHEN("isabs is called") {
            bool is_abs = path::isabs(path);
            
            THEN("returns false") {
                REQUIRE(is_abs == false);
            }
        }
    }
}

SCENARIO("exists tests whether a path exists", "[path]") {

    GIVEN("A bogus file path") {
        std::string path = "../noone/would/ever/have/this/path/on/their/computer/bhe4625e6brt4evr4";

        WHEN("exists is called") {
            bool e = path::exists(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit file path") {
        std::string path = "data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("exists is called") {
            bool e = path::exists(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit file path starting with 'this dir'") {
        std::string path = "./data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("exists is called") {
            bool e = path::exists(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit dir") {
        std::string path = "data";

        WHEN("exists is called") {
            bool e = path::exists(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit dir ending in separtor") {
        std::string path = "data/";

        WHEN("exists is called") {
            bool e = path::exists(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }
}

SCENARIO("isdir tests whether a path is an existing directory", "[path]") {

    GIVEN("A bogus file path") {
        std::string path = "../noone/would/ever/have/this/path/on/their/computer/bhe4625e6brt4evr4";

        WHEN("isdir is called") {
            bool e = path::isdir(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit file path") {
        std::string path = "data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("isdir is called") {
            bool e = path::isdir(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit file path starting with 'this dir'") {
        std::string path = "./data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("isdir is called") {
            bool e = path::isdir(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit dir") {
        std::string path = "data";

        WHEN("isdir is called") {
            bool e = path::isdir(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit dir ending in separtor") {
        std::string path = "data/";

        WHEN("isdir is called") {
            bool e = path::isdir(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }
}

SCENARIO("isfile tests whether a path is a regular file", "[path]") {

    GIVEN("A bogus file path") {
        std::string path = "../noone/would/ever/have/this/path/on/their/computer/bhe4625e6brt4evr4";

        WHEN("isfile is called") {
            bool e = path::isfile(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit file path") {
        std::string path = "data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("isfile is called") {
            bool e = path::isfile(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit file path starting with 'this dir'") {
        std::string path = "./data/diploid-standard-data-ntax5-nchar5.nex";

        WHEN("isfile is called") {
            bool e = path::isfile(path);
            
            THEN("returns true") {
                REQUIRE(e == true);
            }
        }
    }

    GIVEN("A legit dir") {
        std::string path = "data";

        WHEN("isfile is called") {
            bool e = path::isfile(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }

    GIVEN("A legit dir ending in separtor") {
        std::string path = "data/";

        WHEN("isfile is called") {
            bool e = path::isfile(path);
            
            THEN("returns false") {
                REQUIRE(e == false);
            }
        }
    }
}
