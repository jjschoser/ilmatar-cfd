#include <iostream>

#include "EquationOfState.H"
#include "Euler.H"
#include "FluxSolver.H"
#include "FileHandler.H"
#include "Mesh.H"
#include "Reconstruction.H"
#include "Solver.H"
#include "TestProblems.H"

int main(int argc, char *argv[])
{
    std::string settingsFilename, initFilename, finalFilename;
    REAL finalTime, outInterval;
    std::array<std::array<BoundaryCondition, GRIDDIM>, 2> bc;
    REAL gamma = 1.4;
    #ifdef USE_RIGID
        std::string sdfFilename;
    #endif

    if(argc >= 2)
    {
        settingsFilename = argv[1];
        std::ifstream file(settingsFilename);
        assert(file.is_open());
        std::string finalTimeLine, loBCLine, hiBCLine, gammaLine, outIntervalLine;
        std::getline(file, initFilename);
        std::getline(file, finalFilename);
        std::getline(file, finalTimeLine);
        std::getline(file, loBCLine);
        std::getline(file, hiBCLine);
        std::getline(file, gammaLine);
        std::getline(file, outIntervalLine);
        std::istringstream finalTimeISS(finalTimeLine);
        std::istringstream loBCISS(loBCLine);
        std::istringstream hiBCISS(hiBCLine);
        std::istringstream gammaISS(gammaLine);
        std::istringstream outIntervalISS(outIntervalLine);
        finalTimeISS >> finalTime;
        int loBC, hiBC;
        for(int d = 0; d < GRIDDIM; ++d)
        {
            loBCISS >> loBC;
            hiBCISS >> hiBC;
            bc[0][d] = static_cast<BoundaryCondition>(loBC);
            bc[1][d] = static_cast<BoundaryCondition>(hiBC);
        }
        gammaISS >> gamma;
        outIntervalISS >> outInterval;
        #ifdef USE_RIGID
            std::getline(file, sdfFilename);
        #endif
        file.close();
    }

    const IdealGas eos(gamma);
    const Euler euler(&eos);
    const HLLCSolver fluxSolver(euler);
    const MUSCLHancock recon(euler);

    if(!settingsFilename.empty())
    {
        const int nGhost = 1 + recon.getStencilSize();
        int startStep;
        REAL startTime;
        Mesh<Euler::NVARS> mesh = Mesh<Euler::NVARS>::createFromFile(addPath(settingsFilename, initFilename), startStep, startTime, nGhost);
        #ifdef USE_RIGID
            if(!sdfFilename.empty())
            {
                assert(mesh.loadSDF(addPath(settingsFilename, sdfFilename)));
            }
        #endif
        const REAL cfl = 0.9;  // TODO: allow user to adjust this in settings file
        solve(euler, finalTime, mesh, bc, &fluxSolver, &recon, addPath(settingsFilename, finalFilename), cfl, outInterval, startStep, startTime);
    }
    else
    {
        std::cout << "Running test problems..." << std::endl;
        #if GRIDDIM == 1
            const std::array<int, GRIDDIM> res = {2048};
        #elif GRIDDIM == 2
            const std::array<int, GRIDDIM> res = {512, 512};
        #else  // GRIDDIM == 3
            const std::array<int, GRIDDIM> res = {128, 128, 128};
        #endif

        runSimpleTest(euler, &fluxSolver, &recon, res);
        #if GRIDDIM == 2
            runKelvinHelmholtzTest(euler, &fluxSolver, &recon, res);
            #ifdef USE_RIGID
                const std::array<int, GRIDDIM> shockReflectionRes = {2 * res[0], res[1]};
                runShockReflectionTest(euler, &fluxSolver, &recon, shockReflectionRes);
            #endif
        #endif

        #if GRIDDIM == 3
            #ifdef USE_RIGID
                const std::array<int, GRIDDIM> hypersonicSphereRes = {GRIDDIM_DECL(res[0], 2 * res[1], 2 * res[2])};
                runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, false);
                runHypersonicSphereTest(euler, &fluxSolver, &recon, hypersonicSphereRes, true);

                const std::array<int, GRIDDIM> wingRes = {GRIDDIM_DECL(2 * res[0], res[1], res[2])};
                runWingTest(euler, &fluxSolver, &recon, wingRes);

                const int spaceShuttleRes1D = std::max(512, 4 * res[0]);
                const std::array<int, GRIDDIM> spaceShuttleRes = {GRIDDIM_DECL(spaceShuttleRes1D, spaceShuttleRes1D, spaceShuttleRes1D)};
                runSpaceShuttleTest(euler, &fluxSolver, &recon, spaceShuttleRes);
            #endif
        #endif
    }

    return 0;
}