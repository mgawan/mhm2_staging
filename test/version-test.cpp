#include "gtest/gtest.h"
#include "version.h"

#include <string>
using std::string;
TEST(MHMTest, version) {
    //arrange
    //act
    //assert
    string ver(MHM2_VERSION);
    
    auto first = ver.find_first_of('.');
    auto second = ver.find_first_of('.', first+1);
    EXPECT_STREQ (ver.substr(0,first).c_str(),  "2") << "Major Version is correct";
    EXPECT_STREQ (ver.substr(first+1,second-first-1).c_str(),  "0") << "Minor Version is correct";
    EXPECT_TRUE( ver.find_first_of('.',second+1) != string::npos) << "Has a PATCH"; // no verification on this value
}
