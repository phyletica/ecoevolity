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

    GIVEN("A file name only") {
        std::string path = "file.txt";

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
}
