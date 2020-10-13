#include "gtest/gtest.h"
#include "version.h"

#include <string>
using std::string;
TEST(MHMTest, version) {
    //arrange
    //act
    //assert
    string ver(MHM2_VERSION);
    
    EXPECT_STREQ (ver.substr(0,ver.find_last_of('.')).c_str(),  "2.0.1") << "Version is correct";
}
