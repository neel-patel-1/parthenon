//========================================================================================
// Parthenon performance portable AMR framework
// Copyright(C) 2021-2022 The Parthenon collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
// (C) (or copyright) 2020-2022. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================



#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "basic_types.hpp"
#include "config.hpp"
#include "globals.hpp"
#include "interface/update.hpp"
#include "kokkos_abstraction.hpp"

#include "data_pack.hpp"
#include "parthenon_manager.hpp"

/* Package Code */
namespace data_pack {

Packages_t ProcessPackages(std::unique_ptr<ParameterInput> &pin) {
  Packages_t packages;
  packages.Add(data_pack::Initialize(pin.get()));
  return packages;
}

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("data_pack");
  Params &params = pkg->AllParams();

  int n3 = 10;
  int n2 = 20;
  int n1 = 30;

  auto m = Metadata({Metadata::Tensor, Metadata::Requires});
  ParArrayND<double> arr_3d("my_tensor",n3, n2, n1);

  return pkg;
}

TaskStatus DoReduction(Packages_t &packages){
  Real r_val;
  r_val = 10.0;
  packages.Get("data_pack")->AddParam("reduce_val", r_val);
  return TaskStatus::complete;
}

} // namespace data_pack

/* Driver Code */
namespace DataPack{

using namespace parthenon::driver::prelude;

  DriverStatus DataPackDriver::Execute(){
    auto &r_val = pmesh->packages.Get("data_pack")->Param<Real>("reduce_val");
    std::cout<<"Reduction Complete: "<< r_val <<std::endl;
    return DriverStatus::complete;
  }
 
template <typename T>
TaskCollection DataPackDriver::MakeTaskCollection(T &blocks) {
  using data_pack::DoReduction;
  TaskCollection tc;

  TaskRegion &sync_region = tc.AddRegion(0);
  {
    TaskID none(0);
    auto perform_reduce =
        sync_region[0].AddTask(none, DoReduction, pmesh->packages);
  }

  return tc;
}
}// namespace data_pack


/* Main */
int main(int argc, char *argv[]) {
  using parthenon::ParthenonManager;
  using parthenon::ParthenonStatus;
  ParthenonManager pman;

  // Redefine parthenon defaults
  pman.app_input->ProcessPackages = data_pack::ProcessPackages;

  // call ParthenonInit to initialize MPI and Kokkos, parse the input deck, and set up
  auto manager_status = pman.ParthenonInit(argc, argv);
  if (manager_status == ParthenonStatus::complete) {
    pman.ParthenonFinalize();
    return 0;
  }
  if (manager_status == ParthenonStatus::error) {
    pman.ParthenonFinalize();
    return 1;
  }
  // Now that ParthenonInit has been called and setup succeeded, the code can now
  // make use of MPI and Kokkos

  // This needs to be scoped so that the driver object is destructed before Finalize
  {
    // Initialize the driver
    DataPack::DataPackDriver driver(pman.pinput.get(), pman.app_input.get(),
                                              pman.pmesh.get());

    // This line actually runs the simulation
    auto driver_status = driver.Execute();
  }
  // call MPI_Finalize and Kokkos::finalize if necessary
  pman.ParthenonFinalize();

  // MPI and Kokkos can no longer be used

  return (0);
}
