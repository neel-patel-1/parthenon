//========================================================================================
// Parthenon performance portable AMR framework
// Copyright(C) 2021 The Parthenon collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
// (C) (or copyright) 2020-2021. Triad National Security, LLC. All rights reserved.
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
#ifndef EXAMPLE_PARTICLE_LEAPFROG_PARTICLE_LEAPFROG_HPP_
#define EXAMPLE_PARTICLE_LEAPFROG_PARTICLE_LEAPFROG_HPP_

#include <memory>

#include "Kokkos_Random.hpp"

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

namespace data_pack {
  using namespace parthenon::package::prelude;
  using parthenon::Packages_t;
  using parthenon::ParArrayHost;
  using Pack_t = parthenon::MeshBlockVarPack<Real>;

  std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);

  parthenon::TaskStatus DoReduction(Packages_t &packages);

} // namespace data_pack

namespace DataPack{
using namespace parthenon::driver::prelude;
  class DataPackDriver : public Driver {
  public:
    DataPackDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pm)
        : Driver(pin, app_in, pm) {}
    
    /* Not required, but a collection or list of tasks is expected. */
    template <typename T>
    TaskCollection MakeTaskCollection(T &blocks);

    /* Define the drivers execute function*/
    DriverStatus Execute() override;
  };
} //namespace DataPack

#endif // EXAMPLE_PARTICLE_LEAPFROG_PARTICLE_LEAPFROG_HPP_
