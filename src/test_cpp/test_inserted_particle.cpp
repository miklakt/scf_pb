#include "../ss_scf_common.hpp"
#include "../normalization_condition.hpp"
#include "../profiles.hpp"
#include "../topology.hpp"
#include "../particles.hpp"
#include "../inserted_particle.hpp"

#include <iostream>

int main(int argc, char* argv[]){
    particle::Cylinder particle{4,4};
    BrushProfilePlanar brush{0.0, 1000.0, 0.02, topology::kappa(1000.0)};

    InsertedParticleProfiles iparticle{&particle, &brush};
    std::cout << iparticle.osmotic_free_energy(10);

}