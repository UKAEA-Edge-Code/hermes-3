#include <cmath>

#include <bout/bout_types.hxx>
#include <bout/region.hxx>

#include "gtest/gtest.h"

#include "../../include/vantage.hxx"
#include "../../include/component.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

class VantageTest : public FakeMeshFixture {
public:
  VantageTest()
      : FakeMeshFixture(), options({{"units",
                                     {{"eV", 1.0},
                                      {"meters", 1.0},
                                      {"seconds", 1.0},
                                      {"inv_meters_cubed", 1e19}}}}),
        component("test", options, nullptr) {}
  Options options;
  Vantage component;
};

// Test that the component can be created
TEST_F(VantageTest, CreateComponent) {
  Options options;

  Vantage const component("test", options, nullptr);
}