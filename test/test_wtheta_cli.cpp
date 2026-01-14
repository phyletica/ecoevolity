#include "catch.hpp"

#include "ecoevolity/wtheta.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"
#include "ecoevolity/spreadsheet.hpp"


/**
 * These are sanity tests of 'wtheta's output. The do not test the accuracy of
 * Watterson's theta estimates; that is done by test_wattersons_theta.cpp.
 */

TEST_CASE("Testing wtheta cli with Cytrodactylus RAD data",
        "[wtheta]") {

    SECTION("Testing Cytodactylus RAD data from tutorial") {
        RandomNumberGenerator rng = RandomNumberGenerator(48937459);

        std::string config_path = "data/Cyrtodactylus-tutorial-data.yml";
        std::string out_path = "data/Cyrtodactylus-tutorial-data.out";
        std::string expected_out_path = "data/Cyrtodactylus-tutorial-data-expected-wtheta-output.txt";
        REQUIRE(path::exists(config_path));

        char arg0[] = "wtheta";
        char arg1[] = "--relax-missing-sites";
        char * cfg_path = new char[config_path.size() + 1];
        std::copy(config_path.begin(), config_path.end(), cfg_path);
        cfg_path[config_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            &arg1[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        std::ofstream out_stream(out_path);

        // Store original std::cout buffer
        std::streambuf* orig_cout_buffer = std::cout.rdbuf();
        // Redirect std::cout to file
        std::cout.rdbuf(out_stream.rdbuf());

        ret = wtheta_main<RelativeRootCollectionSettings, ComparisonRelativeRootPopulationTreeCollection, BasePopulationTree>(argc, argv);

        // Restore std::out to console
        std::cout.rdbuf(orig_cout_buffer);

        REQUIRE(ret == 0);

        REQUIRE(path::exists(out_path));

        spreadsheet::Spreadsheet wtheta_table;
        wtheta_table.update(out_path);

        spreadsheet::Spreadsheet expected_wtheta_table;
        expected_wtheta_table.update(expected_out_path);

        std::vector<std::string> pops;
        std::vector<std::string> expected_pops;

        std::vector<unsigned int> max_num_alleles;
        std::vector<unsigned int> expected_max_num_alleles;

        std::vector<double> thetas;
        std::vector<double> expected_thetas;

        pops = wtheta_table.get<std::string>("population");
        expected_pops = expected_wtheta_table.get<std::string>("population");

        REQUIRE(pops.size() == expected_pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(pops.at(i) == expected_pops.at(i));
        }

        max_num_alleles = wtheta_table.get<unsigned int>("max_number_of_alleles");
        expected_max_num_alleles = expected_wtheta_table.get<unsigned int>("max_number_of_alleles");

        REQUIRE(max_num_alleles.size() == expected_max_num_alleles.size());
        REQUIRE(max_num_alleles.size() == pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(max_num_alleles.at(i) == expected_max_num_alleles.at(i));
        }

        thetas = wtheta_table.get<double>("wattersons_theta");
        expected_thetas = expected_wtheta_table.get<double>("wattersons_theta");
        
        REQUIRE(thetas.size() == expected_thetas.size());
        REQUIRE(thetas.size() == pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(thetas.at(i) == Approx(expected_thetas.at(i)).epsilon(1e-9));
        }

        delete[] cfg_path;
    }
}

TEST_CASE("Testing wtheta cli with hemi ecoevolity config",
        "[wtheta]") {

    SECTION("Testing hemi ecoevolity config") {
        RandomNumberGenerator rng = RandomNumberGenerator(641645);

        std::string config_path = "data/hemi129-config.yml";
        std::string out_path = "data/hemi129-config.out";
        std::string expected_out_path = "data/hemi129-config-expected-wtheta-output.txt";
        REQUIRE(path::exists(config_path));

        char arg0[] = "wtheta";
        // char arg1[] = "--relax-missing-sites";
        char * cfg_path = new char[config_path.size() + 1];
        std::copy(config_path.begin(), config_path.end(), cfg_path);
        cfg_path[config_path.size()] = '\0';
        char * argv[] = {
            &arg0[0],
            // &arg1[0],
            cfg_path,
            NULL
        };
        int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
        int ret;

        std::ofstream out_stream(out_path);

        // Store original std::cout buffer
        std::streambuf* orig_cout_buffer = std::cout.rdbuf();
        // Redirect std::cout to file
        std::cout.rdbuf(out_stream.rdbuf());

        ret = wtheta_main<RelativeRootCollectionSettings, ComparisonRelativeRootPopulationTreeCollection, BasePopulationTree>(argc, argv);

        // Restore std::out to console
        std::cout.rdbuf(orig_cout_buffer);

        REQUIRE(ret == 0);

        REQUIRE(path::exists(out_path));

        spreadsheet::Spreadsheet wtheta_table;
        wtheta_table.update(out_path);

        spreadsheet::Spreadsheet expected_wtheta_table;
        expected_wtheta_table.update(expected_out_path);

        std::vector<std::string> pops;
        std::vector<std::string> expected_pops;

        std::vector<unsigned int> max_num_alleles;
        std::vector<unsigned int> expected_max_num_alleles;

        std::vector<double> thetas;
        std::vector<double> expected_thetas;

        pops = wtheta_table.get<std::string>("population");
        expected_pops = expected_wtheta_table.get<std::string>("population");

        REQUIRE(pops.size() == expected_pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(pops.at(i) == expected_pops.at(i));
        }

        max_num_alleles = wtheta_table.get<unsigned int>("max_number_of_alleles");
        expected_max_num_alleles = expected_wtheta_table.get<unsigned int>("max_number_of_alleles");

        REQUIRE(max_num_alleles.size() == expected_max_num_alleles.size());
        REQUIRE(max_num_alleles.size() == pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(max_num_alleles.at(i) == expected_max_num_alleles.at(i));
        }

        thetas = wtheta_table.get<double>("wattersons_theta");
        expected_thetas = expected_wtheta_table.get<double>("wattersons_theta");
        
        REQUIRE(thetas.size() == expected_thetas.size());
        REQUIRE(thetas.size() == pops.size());
        for (unsigned int i = 0; i < pops.size(); ++i) {
            REQUIRE(thetas.at(i) == Approx(expected_thetas.at(i)).epsilon(1e-9));
        }

        delete[] cfg_path;
    }
}
