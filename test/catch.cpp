#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <cmath>
#include "smooth.h"

TEST_CASE("Smooth model can be instantiated and configured", "[Smooth]") {

  SECTION("Smooth can be constructed") {
    Smooth smooth;
    REQUIRE(smooth.Size() == 10000);
    REQUIRE(smooth.Field().size() == smooth.Size());
  }
}

TEST_CASE("Smooth mathematical functions are correct", "[Smooth]") {
  Smooth smooth;
  SECTION("Disk support function is correct") {
    REQUIRE(smooth.Disk(500) == Approx(0));
    REQUIRE(smooth.Disk(21.6) == 0.0);
    REQUIRE(smooth.Disk(21.4) > 0.0);
    REQUIRE(smooth.Disk(21.4) < 1.0);
    REQUIRE(smooth.Disk(20.6) > 0.0);
    REQUIRE(smooth.Disk(20.6) < 1.0);
    REQUIRE(smooth.Disk(20.4) == Approx(1.0));
    REQUIRE(smooth.Disk(19.0) == Approx(1.0));
    REQUIRE(smooth.Disk(21.0) == Approx(0.5));
  }
  SECTION("Ring support function is correct") {
    REQUIRE(smooth.Ring(22) == 1.0);
    REQUIRE(smooth.Ring(21.6) == 1.0);
    REQUIRE(smooth.Ring(21.4) > 0.0);
    REQUIRE(smooth.Ring(21.4) < 1.0);
    REQUIRE(smooth.Ring(20.6) > 0.0);
    REQUIRE(smooth.Ring(20.6) < 1.0);
    REQUIRE(smooth.Ring(20.4) == 0.0);
    REQUIRE(smooth.Ring(21.0) == 0.5);
    REQUIRE(smooth.Ring(64.0) == 0.0);
    REQUIRE(smooth.Ring(63.6) == 0.0);
    REQUIRE(smooth.Ring(63.4) > 0.0);
    REQUIRE(smooth.Ring(63.4) < 1.0);
    REQUIRE(smooth.Ring(62.6) > 0.0);
    REQUIRE(smooth.Ring(62.6) < 1.0);
    REQUIRE(smooth.Ring(62.4) == 1.0);
    REQUIRE(smooth.Ring(63.0) == 0.5);
  }
  SECTION("Sigmoid function is correct") {
    double e = std::exp(1.0);
    REQUIRE(Smooth::Sigmoid(1.0, 1.0, 4.0) == 0.5);
    REQUIRE(std::abs(Smooth::Sigmoid(1.0, 0.0, 4.0) - e / (1 + e)) < 0.0001);
    REQUIRE(Smooth::Sigmoid(10000, 1.0, 4.0) == 1.0);
    REQUIRE(std::abs(Smooth::Sigmoid(0.0, 1.0, 0.1)) < 0.001);
  }
  SECTION("Transition function is correct") {
    REQUIRE(std::abs(smooth.Transition(1.0, 0.3) - 1.0) < 0.1);
    REQUIRE(smooth.Transition(1.0, 1.0) == Approx(0));
    REQUIRE(std::abs(smooth.Transition(0.0, 0.3) - 1.0) < 0.1);
    REQUIRE(std::abs(smooth.Transition(0.0, 0.0)) < 0.1);
  }
  SECTION("Wraparound Distance is correct") {
    REQUIRE(smooth.TorusDifference(95, 5, 100) == 10);
    REQUIRE(smooth.TorusDifference(5, 96, 100) == 9);
    REQUIRE(smooth.TorusDifference(5, 10, 100) == 5);
    REQUIRE(smooth.Radius(10, 10, 13, 14) == 5.0);
  }
}

TEST_CASE("NormalisationsAreCorrect") {
  Smooth smooth(100, 100, 10);
  SECTION("Disk Normalisation is correct") {
    // Should be roughly pi*radius*radius,
    REQUIRE(std::abs(smooth.NormalisationDisk() - 314.15) < 1.0);
  }
  SECTION("Ring Normalisation is correct") {
    // Should be roughly pi*outer*outer-pi*inner*inner, pi*100*(9-1), 2513.27
    REQUIRE(std::abs(smooth.NormalisationRing() - 2513.27) < 2.0);
  }
}

TEST_CASE("FillingsAreUnityWhenSeeded") {
  Smooth smooth;
  smooth.SeedConstant(0);
  SECTION("DiskFillingUnityWithDiskSeed") {
    smooth.AddDisk();
    REQUIRE(std::abs(smooth.FillingDisk(0, 0) - 1.0) < 0.1);
  }

  SECTION("Disk Filling Zero With Ring Seed") {
    smooth.AddRing();
    REQUIRE(std::abs(smooth.FillingDisk(0, 0)) < 0.1);
  }
  SECTION("RingFillingUnityWithRingSeed") {
    smooth.AddRing();
    REQUIRE(std::abs(smooth.FillingRing(0, 0) - 1.0) < 0.1);
  }
}

TEST_CASE("FillingFieldHasRangeofValues") {
  Smooth smooth(300, 300);
  smooth.SeedConstant(0);
  smooth.AddRing();
  double min = 1.0;
  double max = 0.0;
  for(int x = 0; x < 300; x++) {
    double filling = smooth.FillingRing(x, 0);
    if(filling < min)
      min = filling;
    if(filling > max)
      max = filling;
  }
  REQUIRE(min < 0.2);
  REQUIRE(max > 0.4);
}

TEST_CASE("Compute Integrals") {
  Smooth smooth(300, 300);
  smooth.SeedConstant(0);

  // check for different positions in the torus
  for(auto const x : {150, 298, 0})
    for(auto const y : {150, 298, 0}) {
      SECTION("At position (" + std::to_string(x) + ", " + std::to_string(y) + ")") {
        SECTION("Ring only") {
          smooth.AddRing(150, 150);

          auto const result = smooth.Integrals(150, 150);
          // 0.1 accuracy because of smoothing
          CHECK(std::get<0>(result) == Approx(0).epsilon(0.1));
          CHECK(std::get<1>(result) == Approx(1).epsilon(0.1));
        }

        SECTION("Disk only") {
          smooth.AddDisk(150, 150);
          auto const result = smooth.Integrals(150, 150);
          CHECK(std::get<0>(result) == Approx(1).epsilon(0.1));
          CHECK(std::get<1>(result) == Approx(0).epsilon(0.1));
        }

        SECTION("Disk and ring") {
          smooth.AddRing(150, 150);
          smooth.AddDisk(150, 150);
          auto const result = smooth.Integrals(150, 150);
          CHECK(std::get<0>(result) == Approx(1).epsilon(0.1));
          CHECK(std::get<1>(result) == Approx(1).epsilon(0.1));
        }
      }
    }
}

TEST_CASE("Update") {
  // just test playing with a single pixel lit up sufficiently that the
  // transition is non-zero in the ring.
  auto const radius = 5;
  Smooth smooth(100, 100, radius);
  smooth.AddPixel(50, 50, 0.3 * smooth.NormalisationRing());
  CHECK(std::get<0>(smooth.Integrals(50, 50))
        == Approx(0.3 * smooth.NormalisationRing() / smooth.NormalisationDisk()));

  // check the integrals are numbers for which Transition gives non-zero result
  // in the ring
  CHECK(std::get<1>(smooth.Integrals(50, 50)) == Approx(0));
  CHECK(std::get<0>(smooth.Integrals(40, 40)) == Approx(0));
  CHECK(std::get<1>(smooth.Integrals(40, 40)) == Approx(0.3));
  CHECK(std::get<0>(smooth.Integrals(42, 39)) == Approx(0));
  CHECK(std::get<1>(smooth.Integrals(42, 39)) == Approx(0.3));

  // Now call update
  smooth.Update();
  auto const field = smooth.Field();
  // And check death in the disk
  CHECK(field[smooth.Index(50, 50)] == Approx(0));
  CHECK(field[smooth.Index(51, 52)] == Approx(0));
  // And check life in the ring
  CHECK(field[smooth.Index(45, 45)] == Approx(smooth.Transition(0, 0.3)));
  CHECK(field[smooth.Index(42, 39)] == Approx(smooth.Transition(0, 0.3)));
  // And check death outside
  CHECK(field[smooth.Index(15, 15)] == Approx(0));
}
