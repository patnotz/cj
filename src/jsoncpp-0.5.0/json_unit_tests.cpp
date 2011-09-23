#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <json/json.h>

// This test is mostly meant as a concise example of how to parse and access
// Json data and to make sure the build is working. It's not meant to an
// exhaustive test of the JsonCPP package.
TEST(Json, readFromStream)
{
  const std::string input =
  "	{\"menu\": {"
  "	  \"id\": 27,"
  "	  \"value\": \"File\","
  "   \"pi\": 0.31415926e1,"
  "	  \"popup\": {"
  "	    \"menuitem\": ["
  "	      {\"value\": \"New\", \"onclick\": \"CreateNewDoc()\"},"
  "	      {\"value\": \"Open\", \"onclick\": \"OpenDoc()\"},"
  "	      {\"value\": \"Close\", \"onclick\": \"CloseDoc()\"}"
  "	    ]"
  "	  }"
  "	}}";
  std::istringstream iss(input);
  const Json::Value nullValue;
  Json::Value root;
  Json::Reader reader;
  const bool parsingSuccessful = reader.parse(iss, root);
  ASSERT_TRUE(parsingSuccessful);
  ASSERT_EQ("File", root["menu"]["value"].asString());
  ASSERT_EQ(27, root["menu"]["id"].asInt());
  ASSERT_DOUBLE_EQ(3.1415926, root["menu"]["pi"].asDouble());
  Json::Value popup = root["menu"]["popup"];
  ASSERT_EQ(3, popup["menuitem"].size());
  ASSERT_EQ("Open", popup["menuitem"].get(1,nullValue)["value"].asString());
}
