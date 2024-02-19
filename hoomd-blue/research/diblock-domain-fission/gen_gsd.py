import sys
import numpy as np
import gsd
import gsd.hoomd

class System:

    def __init__(self,n,N=256,b0=0.75):
        self.L = np.asarray([10, 10, 10])
        self.N = N
        self.n = n
        self.b0 = b0

def main(argv):
    if len(argv) != 1:
        print("Usage: ./gen_gsd.py")
        return

    n = 200
    system = System(n)

    snapshot = gsd.hoomd.Snapshot()

    snapshot.configuration.box = np.concatenate((system.L,[0,0,0]))

    snapshot.particles.N = system.n * system.N

    snapshot.particles.types = ["A", "B"]
    snapshot.particles.position = np.zeros( (snapshot.particles.N,3) )
    snapshot.particles.velocity = np.random.standard_normal( (snapshot.particles.N,3) )
    snapshot.particles.image = np.zeros( (snapshot.particles.N,3),dtype=np.int )
    snapshot.particles.typeid = np.zeros( snapshot.particles.N, dtype = np.int )
    snapshot.bonds.N = system.n * (system.N-1)
    snapshot.bonds.types = ["backbone"]
    snapshot.bonds.typeid = np.zeros( snapshot.bonds.N )
    snapshot.bonds.group = []

    for poly in range(system.n):
        mono = 0
        snapshot.particles.position[ poly * system.N + mono ] = [0,-system.L[1]/2,0] + np.random.standard_normal(3)*1e-3

        for mono in range(1, system.N):
            snapshot.particles.position[ poly * system.N +mono] = snapshot.particles.position[poly*system.N + mono-1] + np.random.standard_normal( 3 ) * system.b0/3

            snapshot.bonds.group.append( [ poly*system.N+mono-1 , poly*system.N+mono])
            if mono > 16:
                snapshot.particles.typeid[poly * system.N +mono] = 1

    snapshot.particles.image += np.rint(snapshot.particles.position/system.L).astype(np.int64)
    snapshot.particles.position -= np.rint(snapshot.particles.position/system.L) * system.L

    snapshot.particles.validate()

    with gsd.hoomd.open( "start.gsd","wb") as f:
        f.append(snapshot)


if __name__ == "__main__":
    main(sys.argv)
