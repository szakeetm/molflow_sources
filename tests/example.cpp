//
// Created by pascal on 8/8/19.
//

#include "../src_shared/SimulationManager.h"
#include "gtest/gtest.h"
//#define MOLFLOW_PATH ""

namespace {

// Tests factorial of positive numbers.
    TEST(SubProcessInit, Zero) {

        {
            SimulationManager simMan("molflow","MFLW");
            EXPECT_EQ(0, simMan.InitSimUnits());
        }

        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbCores = 0;
            simMan.InitSimUnits();
            EXPECT_EQ(0, simMan.simHandles.size());
        }

        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbCores = 4;
            simMan.InitSimUnits();
            EXPECT_EQ(4, simMan.simHandles.size());
            EXPECT_EQ(0, simMan.ExecuteAndWait(COMMAND_EXIT,PROCESS_KILLED));
        }
    }

    TEST(SubProcessCreateAndKill, CPU) {
        {
            SimulationManager simMan("molflow","MFLW");
            simMan.useCPU = true;
            simMan.nbCores = 4;
            simMan.InitSimUnits();
            EXPECT_EQ(4, simMan.simHandles.size());
            simMan.KillAllSimUnits();
            EXPECT_EQ(0, simMan.simHandles.size());
        }
    }

}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);



    return RUN_ALL_TESTS();
}