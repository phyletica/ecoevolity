#include "catch.hpp"
#include "ecoevolity/simple_nuc_data.hpp"

TEST_CASE("Simple test of NucData", "[NucData]") {

    SECTION("splitting string into pre-allocated vector") {
        std::stringstream in_stream;
        in_stream << "   " << std::endl;
        in_stream << "4 17" << std::endl;
        in_stream << "seq_a     ACGTacgtRyMNN--??" << std::endl;
        in_stream << "seq_b     SsKkacgtRyMNN--Xx" << std::endl;
        in_stream << "seq_c     ---TacgtRyMNN--??" << std::endl;
        in_stream << "seq_d     HhBbDdVvWwUu-?Nnn" << std::endl;

        ecoevolity::NucData d;
        d.init_from_phylip_stream(in_stream);
        REQUIRE(d.get_seq("seq_a") == "ACGTacgtRyMNN--??");
        REQUIRE(d.get_seq("seq_b") == "SsKkacgtRyMNN--Xx");
        REQUIRE(d.get_seq("seq_c") == "---TacgtRyMNN--??");
        REQUIRE(d.get_seq("seq_d") == "HhBbDdVvWwUu-?Nnn");
    }
}
