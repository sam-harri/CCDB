Changelog
=========

v5.8.0
------
**Library**
- Use GLL quadrature points on Tri Dirichlet BCs with CG (!1895)
- Add FieldConvert module to perform local stability analysis for compressible flows (!1319)
- Remove get() accessor from Array data structure (!1937)
- Fix issue with `StdTetExp::v_LocCollapsedToLocCoord` (!1946)
- Fix issue with `NodalTriExp::v_GetTracePhysVals` (!1951)
- Remove deprecated version of v_SetCoeffsToOrientation functions (!1954)
- Fix issue with Dirichlet BCs when using variable P (!1972)
- Patch for implicit-function error in scotch-6.0.4 (!1938)
- Tidy virtual inheritance in NodalTriExp (!1979)
- Fix partial overload virtual function in AssemblyMap, StdRegions, and LocalRegions (!1978)
- Fix partial overload virtual function v_PhysEvaluate (!1980)
- Matrix free ops shape cleanup (!1735) 
- Fix NodalTri processing and static condensation matrix release (!1989)
- Fix third-party Scotch patch (!1998)
- Fixed Modified Arnoldi driver to remove discontinuities from random vectors (!2002)
- Addedd support for backing up filters output when the file name have already existis (!2006)
- Partially revert Geometry::v_ContainsPoint (!2007)

**CI**
- Fix CubeAllElements performance test tolerance (!1943)
- Remove `allow_failure` from compiler warnings and formatting (!1958, !1966)
- remove CI image tag when dockerhub deploy completes (!1960)
- Use recursive strategy for submodule (!1997)
- cleanup CI environment images after packaging (!1991)

**NekMesh**
- Add high-order pyramid and prism support from gmsh (!1956)

**Python**
- Transition bindings to use pybind11 (!1950)

**Documentation**
- Updated the User-guide with additional inofrmation for outflow BC, addressing the issue #103 (!1988)
- Updated the User-guide with additional inofrmation for outflow BC, addressing the issue #103 (!1990)


v5.7.0
-----
**Library**
- Modified MatrixFreeOp library  switch initialisation to use BOOST_PP (!1794)
- Fix memory-leak with LowEnergyBlock preconditioner for time-updated matrices (!1627)
- Fix Fourier expansion integration weights are related test (!1803)
- Introduced the MatrixFree implementation for the PhysInterp1DScaled operator and tidying up MatrixFreeOps (!1812) 
- Separate MeshGraph input/output functions into a new class (!1778)
- Added checkpoint file writing start time in the fieldconvert filter (!1789)
- Fix fieldconvert filter incorrect boundary values (!1789)
- Fix numerical precision issues with filters OutputStartTime (!1789)
- Fix AdaptiveSFD for MPI (!1821)
- Fix interpolation on manifold (!1840)
- Fix IterativeStaticCond when using absolute tolerance (!1850)
- Fix deadlock by scotch with multi-threading support (!1853)
- Fixed L2norm for FilterError (!1871)
- Fix variable p in tetrahedrons (!1881)
- Fix BwdTrans for Pyr with var P (!1886)
- Allow wrapper array around a existing raw pointer (!1848)
- Tweaked some long tests to make them faster (!1918)

**IncNavierStokesSolver**
- Fix initial and boundary conditions in the moving reference frame (!1692, !1820)
- Fix memory-leak for the Mixed_CG_Discontinuous projection when initializing the traceMep (!1806)
- Add synthetic turbulence generation for the incompressible solver (!1664) 
- Fix a uninitialized parameter in VCS (!1880)

**ShallowWaterSolver**
- Implement implicit time-discritization (!1784)

**CompressibleSolver**
- Add synthetic turbulence generator for the compressible solver (!1859)

**NekMesh**
- Added revolve module (!1825)
- Fix Prism Reordering in Process PerAlign (!1899)
- Extend quality measures in ProcessJac and add histogram generation(!1751)
- Reducing run time of some tests in NekMesh(!1922)
- Extend quality measures in ProcessJac and add histogram generation (!1751)
- Reducing run time of some tests in NekMesh (!1922)
- Added a reader for the CGNS input format (!1889)

**FieldConvert**
- Add vortexinducedvelocity module to compute the vortex-induced velocity (!1824)
- Add a module to transform coordinates and vectors for the moving reference frame method (!1830)

**Miscellaneous**
- Use std::stod instead of boost::lexical_cast<NekDouble> (!1819)

**Documentation**
- Add initial documentation for the IncNavierStokesSolver (!1822)
- Updated the supported packages in Userguid (!1904)
- Added a example for RayleighBenardConvection in the user-guide for IncNS (!1919)
- Fix some typos in tutorials (!1929)

**CI and Packaging**
- Debian 10 (BUSTER) is no longer supported (!1902)
- Support is added for Ubuntu Noble Numbat and droped for Bionic Beaver (!1910)
- Removed Fedora 35/36, added Fedora 39/40 (!1909)

v5.6.0
------
**Library**
- Clean-up Set_Rhs_Magnitude function in NekLinSysIter (!1729)
- Consistently use template parameters in VmathArray (!1748)
- Fix issue with CMake and zlib versions >= 1.3.0 (!1744)
- Add 1D demo and test of h-type convergence for a CG projection. (!1738)
- Add 2D projection demo and tests following 1D added in MR !1738. (!1762)
- Tidy up tolerance in NekLinSystIter and NekNonlinSysIter solvers (!1722)
- Enable varcoeffs for Collections (!1701)
- Fix misplaced " in Nektar++Config.cmake (!1742)
- Further tidy-up in linear solver (!1761)
- Use FwdTrans in UnsteadySystem when using Parareal (!1785)
- Automate deployment of README.md to dockerhub (!1786)
- Fix memory leak with Block preconditioner for time-updated matrices (!1737)
- Support for implicit sliding meshes (!1787)
- Fix compilation issue with OpenCASCADE 7.8.0 (!1799)
- Fix MPI communicator for Parallel-in-Time (!1801)
- Fix warnings issued by MSVC compiler in LibUtilities and StdRegions (!1740)
- Fix Parareal convergence monitoring (!1802)
- Avoid copying input/output to/from padded aligned arrays in MatrixFree operators(!1783)

**CompressibleFlowSolver**
- Complete second Frechet derivative implementation (!1761)
- Add conditional updating of elemental Mass and Laplacian matrices for LinearADR matrices (!1766)
- Added routine to order expansion in an optimal manner for MatrixFree/Collection ops (!1770)
- Fix PFASST I/O and pre-initialize coarse preconditioner for Parareal (!1749)
- Remove collection offset arrays since no longer required (!1771)
- Fix summary output (!1779)
- Update Docker images to use bookworm (!1775)
- Update `clang-tidy` and `clang-format` to v16 (!1777)
- Add coverage metric capturing (!1776)

**ShallowWaterSolver**
- Refractoring to reduce code duplication (!1782)

**NekPy**
- Add binding to NekPy to check of geometry elements are valid (!1755)
- Update NekPy to more modern packaging (!1747)
- Add wrapper for selected SolverUtils classes, particularly Filter (!1379)
- Add VTK support for OutputVtk module (!1379)
- Add bindings for SolverUtils::EquationSystem and UnsteadySystem (!1752)
- Return native Python types when getting session variables and parameters (!1380)

**NekMesh**
- Fix optiKind flags in VarOpti for freenodes that are on more than a single curve / surface (!1597)
- Fix VarOpti Surface Node Sliding on the CAD in 2D  (!1569)
- Add feature for r-adaption on user-defined CAD curves (!1349)
- Add feature for r-adaption on user-defined CAD curves (!1349)
- Add unit testing infrustructure and initial example (!1753)
- Added a custom cmake cache file to load defaults for building only NekMesh without the solvers (!1641)
- Extend peralign module to all types of meshes (!1702)

**IncNavierStokesSolver**
- Matrix-Free LinearADR operator for VCSImplicit and others (!1627)
- Make substepping normal velocity evaluation more efficient (!1795)

**FieldConvert**
- Add range function as an option in Xml Input, align python usage and start depracation of -r option (!1791) 

v5.5.0
------
**Library**
- Fix Nektar++Config.cmake to use MPI_CXX (!1224)
- Redesign of Parareal and PFASST driver (!1613)
- Update default global system solver parameters for paralell-in-time (!1649)
- Add member function in MPI communicator to check if time-parallel is enable (!1647)
- Fixed FilterError for homogeneous expansions (!1640)
- Fix ForcingAbsorption for homogeneous expansions (!1650)
- Update AssemblyMap to reduce verbosity when using parallel-in-time (!1651)
- Tidy-up of Collection library (!1622)
- Tidy-up of I/O in BasicUtils (!1623)
- Fix local implementation of CG with Null preconditioner (!1658)
- Update to use C++17 nested namespaces (!1556, !1670)
- Fix a minor bug in ProcessWallNormalData (!1663)
- Fix Explist::v_GetNormals and GetElmtNormalLength (!1625)
- Replace `boost::random` with `std::random` (!1673)
- Add a new feature of Lagrangian points tracking in parallel (!1666)
- Add Robin BC to xxt full solver (!1679)
- Update Parareal file output and tidy-up Parareal/PFASST interpolation (!1678)
- Move from `boost::filesystem` to `std::filesystem` (!1674)
- Use `[[deprecated]]` attribute (!1682)
- Fix QuadExp::v_ComputeTraceNormals(!1685)
- Replace `boost::thread` with `std::thread` (!1687)
- Replace `boost::math::tgamma` by `std::tgamma`(!1686)
- Remove `#include <boost/math/special_functions/fpclassify.hpp>` header (!1686)
- Replace `boost::regex` with `std::regex` (!1676)
- Update SDC scheme for implicit PFASST (!1659)
- Pre-allocate memory for GMRES (!1668)
- Update hdf5 read to only read in selected elements to reduce memory overhead (!1351)
- Remove arbitrary factor in  `GlobalLinSysIterative.cpp` (!1694)
- Some tidy-up in LinearAlgebra (!1699)
- Some further tidy-up in LinearAlgebra (!1700)
- Updated hdf5 to 1.12.3 (!1696)
- Remove redundant tolerance limiter in GMRES (!1707)
- Remove unused tolerance parameter in NekSys class and subclasses (1708)
- Avoid repeatly operator assignment in NekNonlinSysNewton class (!1709)
- Add an exact solution for GetLocCoords of straight-edge quad elements (!1704)
- Add sliding mesh capability (!1605)
- Remove MaxIterations parameter from AssemblyMap (!1710)
- Consistently use relative tolerance for GMRES (!1706)
- Fix use of absolute tolerance for iterative solvers (!1711)
- Move solver info loading form NekSys (and its derived classes) to the callers (!1716)
- Tidy physderiv and helmholtz in MatrixFreeOps (!1653)
- Remove time estimation for Parareal and PFASST drivers (!1691)
- Tidy and clarify implementation of Dirichlet boundary condition in GlobalLinSys (!1724)
- Remove unused file NekLinAlgAlgorithms.hpp (!1728)

**ADRSolver**
- Add support for spatially-constant, but variable direction, diffusion to
  ADRSolver. (!1669)
- Inline Vmath library (!1667)
- Refactor UnsteadyReactionDiffusion as a subclass of UnsteadyDiffusion (!1713)
- Tidy-up ADRSolver (!1715)
- Refactor UnsteadyAdvectionDiffusion as a subclass of UnsteadyAdvection (!1720)
- Update UnsteadyReactionDiffusion solver for explicit time-stepping (!1731)
- Refactor UnsteadyViscousBurgers as a subclass of UnsteadyInviscidBurgers solver (!1718)

**CardiacEPSolver**
- Fix cell model history point filter output after base class change (!1342)
- Add const qualifier to SetUniversalUniqueMap (!1644)
- Add safety check for FinTime parameter for parallel-in-time (!1652)

**IncNavierStokesSolver**
- Save BndElmtExpansion and avoid re-building (!1648)
- Add Simo-advection and a switch for Simo-/Dong-advection to VCSImplicit (!1630)
- Add a new feature of elasticaly mounted object using the moving reference frame (!1495)
- Some tidy-up (!1693)
- Fix memory leak in VCSImplicit due to matrix updating (!1688)
- Rename Simo-/Dong-advection to Extrapolated/Updated respectively for VCSImplicit (!1717)
- Update test for SuccessiveRHS parameter (!1727)
- Fix segmentation error in IncNavierStokesSolver when AeroForces filter is used (!1726)
- Added Examples in solvers/IncNavierStokesSolver/Examples, matching with the user-guide (!1723)
- Extend support for IO timer levels to IncNavierStokesSolver (!1732)
- Fixed issue with substepping when using mixed BCs (!1639)

**CompressibleFlowSolver**
- Add three subsonic inflow boundary conditions: EnforceEntropyPresure, EnforceEntropyVelocity, and EnforceEntropyTotalEnthalpy (!1398)
- Update implicit solver for implicit PFASST (!1659)
- Some tidy-up in the compressible flow solver (!1698)
- Fix ESDIRK scheme for compressible flow solver (!1690)
- Some further tidy-up in compressible flow solver (!1700)

**Documentation**
- Update instructions for VS2022 (!1714)
- Update ADRSolver user guide (!1730)

**ShallowWaterSolver**
- Reduce memory footprint of the Peregrine solver(!1680)
- Fix Peregrine solver for Parareal output(!1689)

**FieldConvert**
- Fix typo in user-guide and fix but for parallel-in-time FieldConvert (!1645)
- Fixed FieldConvert -m addfld (!1500)
- Add a new FieldConvert module to zero a plane in wavespace (!1632)
- Fix tecplot output for line and plane points (!1497)

**CI**
- Disable macOS Intel runner (!1655)
- Upgrade Win10 runners (!1656)
- Add `clang-tidy` support to CI for static analysis (!1556)
- Remove some uncompiled files and add quality step to check (!1672)
- Add `clang-15` tester (!1671)
- Update performance tests to use minimum execution time (!1733)
- Fix missing dependencies for clang-tidy and uncompiled-files jobs (!1736)

**NekMesh**
- Replace deprecated boost::filesystem functions (!1654)
- Remove deprecated AddTraceIntegral interface from ExpList.h (!1646)

**Packaging**
- Remove unnecessary Boost dependencies after C++17 migration (!1734)

**Miscellaneous**
- Add a particle tracker utilitiy using equation system infrastructure (!1310)
- Remove deprecated fs::copy_directory function call (!1662)
- Added a sample python script for wallnormaldata module (!1663)
- Replace boost::lexical_cast with std::to_string for integer variable (!1677)
- Use `[[maybe_unused]]` attribute in solvers (!1683)
- Use `[[maybe_unused]]` attribute in library (!1684)
- Fix sse2 SIMD bug (!1559)
- Improved the Incompressible NS section of User-guide (!1719)

v5.4.0
------
**Library**
- Add unit tests for Polylib (!1533)
- Add the absolution tolerance for the iterative linear system solver (!1549)
- Redesign of Parareal driver (!1509) 
- Add PFASST driver (!1501)
- Add local p-refinement functionality (!1508)
- Tidy-up SDC schemes (!1550)
- Add min/max and fmin/fmax function to the interpretor registry (!1552)
- Fix a _m128d to _m128i conversion bug in sse2.hpp (!1551)
- Register TimeIntegrationMethod to SessionReader (!1539)
- Remove unecessary for loop in NekLinSysIterGMRES (!1554)
- Generalize the use of size_t in TimeIntegration (!1555)
- Fix convergence error in Fixed-point Jacobi algorithm (!1561)
- Tidy-up NekLinSys (!1562)
- Tidy-up TimeIntegration (!1505)
- Fix Adams-Bashforth time integration scheme (!1558)
- Remove unused Kernel directory (!1560)
- Remove unused files in BasicUtils (!1564)
- Fix convergence error in Fixed-point Jacobi algorithm (!1561)
- Update parallel-in-time code base in fieldconvert (!1565)
- Tidy-up Basis and Points classes (!1568)
- Remove unused function GetKey in NekFactory (!1567)
- Add new functionality and code improvements for the local p-refinement (!1571)
- Fix v_Exscan compilation bug (!1566)
- Tidy-up of comm class (!1572)
- Tidy-up Nektar tolerances (!1570)
- Tidy-up UnsteadySystem (!1579)
- Fix StaticCond approach + GMRES when restarting (!1583)
- Tidy-up NekManager (!1587)
- Slightly reduce memory allocation in GMRES (!1586)
- Tidy-up Advection class (!1580)
- Tidy-up diffusion class and subclasses (!1581)
- Tidy-up NekNonlinSys and fix some issues (!1563)
- Avoid pre-allocating memory in GMRES (!1591)
- Added a local space version of GMRES and ConjugateGradient (!1575)
- Added a jacobi preconditioner building on diagonal impelemtation (!1575)
- Add a priori convergence and null-input check to GMRES (!1590)
- Some tidy-up in StdRegions (!1595)
- Tidy-up/removed unused Lapack function definition (!1573)
- Some Tidy-up in MatrixFreeOps (!1599)
- Accelerate interpolation for regular and straight-edge elements (!1283)
- Fix an indexing error in MatrixFreeOps (!1602)
- Update Session File for Parallel-in-Time (!1516)
- Partially fix hdf5 partition for Parareal (!1611)
- Update Collection for Parareal-in-Time (!1607)
- New version of CreateCollection using basisKey (!1603)
- Remove unused Domain.cpp and Domain.h file (!1609)
- Remove useless ReadExpressions and SubstituteExpressions function in sessionReader (!1608)
- Corrected workspace size in triangle BwdTrans (!1610)
- Reactivate Reactivate Movement_fixed_3D_stacked_cylinders_curved_hdf5_par test except on ARM MacOS (!1536)
- Updated the PETSc version to v3.19.3 for THIRDPARTY_USE_PETSC, enabled local PETSc version usage (!1618)

- Rename communicator in LinearAlgebra (!1612)
- Add IProductWRTDerivBase operator for 3DH1D problems (!1483)
- Full support of mixed-order elements in DG (!1606)
- Refactoring collections for allowing variable output size inside the collections' operators, introduced PhysInterp1DScaled in Collections (!1620)

- Use default keyword for destructor in Driver (!1624)
- Add additional test for SDC time-integration(!1621)
- Fix to compiler flags for MSVC (!1604)
- Correct bug in scotch initialisation for substructuring (!1634)
- Add ability to build up Movement objects programmatically (!1600)
- Write out movement data to XML files (!1600)

**CompressibleFlowSolver**
- Fix AUSM3 Riemann solver and add tests to the CI (!1537)
- Set initial condition check points files to zero (0) (!1535) 
- Prevent file overwrite with restart for the Compressible flow solver (!1540)
- Register Projection and GJPStabilisation in SessionReader (!1542)
- Redesign of AUSM Riemann solver class (!1577)
- Redesign of the compressible flow solver. Removal of RinglebFlow and IsentropicVortex subclasses (!1584)
- Redesign of PreconCfs class (!1578)
- Fix virtual functions in ContField class (!1616)
- Tidy-up Compressible flow solver print status (!1615)
- Update of for parallel-in-time (!1589)
- Fix some memory bugs for implicit compressible flow solver with LDGNS (!1617)

**IncNavierStokesSolver**
- Add an option to mask variables for the linear stability problem (!1280)

**ShallowWaterSolver**
- Fix NonlinearPeregrine solver due to a change of API (!1637)

**FieldConvert**
- Add option to use .csv files in pointdatatofld module (!1545)
- Add a new module to output power spectral in a given area (!1271)
- Add a new module to do field averaging (!1271)

**IncNavierStokesSolver**
- Register SolverType in SessionReader (!1541)
- Prevent file overwrite with restart for the IsentropicVortex solver (!1543)
- Fix GJP stabilisation for curved 3D elements (!1593)
- Enable SVV and GJP stabilisation for unresolved scales in VCSImplicit (!1592)
- Fix tolerance in KovaFlow_m10_VCSImplicit_SVV (!1619)

**Miscellaneous**
- Fix compilation against TetGen 1.6 (!1547)
- Tidy-up Parareal initial condition output (!1534)
- Remove double entry in documentation and slightly tidy-up code (!1530)
- Fixed deadlock problem for the LowEnergyBlock preconditioner which appeared for specific combinations of meshes and boundary conditions (!1553)
- Add a check to avoid unnecessary copy in DoOdoProjecton function (!1582)
- Update user guide and developer guide (!1598)
- Some various tidy-up (!1585)
- Added support for performance tests (!1614)
- Integrated performance tests into the CI system and enabled a test (!1629)
- Added additional performance tests and updated test documentation (!1631)
- Fix VmathTimer compilation error (!1635)
- Fix SIMD compilation error (!1636)
- Update ShallowWaterSolver performance test (!1643)

**NekMesh**
- Changed CMake to set NEKTAR_USE_THREAD_SAFETY to ON when MeshGen is set to ON (!1546)
- Fix VarOpti Tolerances (!1569)
- Fixed a bug where MeshElement/Tetrahedron did not assign edge IDs in the constructor (!1596)
- Reduce VarOpti memory usage by removing the storage for GetFunctional in NodeOpti (!1633)

**Documentation**
- Fix missing Tikz external package requires for bookworm (!1638)

**CI**
- Add testing and packaging for Debian bookworm (!1638)
- Added a NekMesh docker image and added the image to the CI (!1641)
**NekPy**
- Add bindings for Movement-related classes (!1600)
- Add bindings for various helpful methods in the SpatialDomains
  library (!1600)

v5.3.0
------
**Library**
- Fixed avx512 back-end for SimdLib (!1333)
- Remove unnecessary IterPerExp methods (!1366)
- Added float to scalar and avx2 back-end, disable avx512, sse2, sve (!1255)
- Updated the library to use m_phys and m_coeff as function arguments (!1412)
- Added float and restored SVE back-end for SimdLib (!1373)
- Fix VmathSIMD by adding optional mapping with # of lanes (!1388)
- Added float and restore avx512 back-end for SimdLib (!1387)
- Fix namespace pollution which causes boost 1.74+ errors (!1389)
- Fix missing copy assignment operator warnings in clang 13+ (!1391)
- Added checkpoint file writing start time in the checkpoint filter (!1401)
- Fix boost 1.77 compatibility errors (!1420)
- Replaced depricated "sprintf" with "std::to_string" (!1406)
- Add compatiblity patch to solve conflict between flex 2.6.3 and scotch 6.0.4 (!1410)
- Add Parareal Driver module (!1317)
- Maintenance for C++-17 compatibility: removed std::unaray_function base class due to removal from the std (!1419)
- Fixed the comment of function Vvtvvtp in VmathArray (!1408)
- Add a FieldConvert utility to compute the divergence of the velocity (!1413)
- Added new filter to calculate variables integral on composite mesh (!1409)
- Overload PhysEvaluate to give first derivatives using barycentric
  interpolation (!1323)
- Non-conformal interface support (!1323)
- Fix a I/O issue related to the IO_InfoSteps parameter (!1422)
- Fix a I/O issue related to the IO_CheckSteps parameter (!1423)
- Fix boost 1.77 compatibility errors (!1402)
- Replaced depricated "sprintf" with "std::to_string" (!1406)
- Add compatiblity patch to solve conflict between flex 2.6.3 and scotch 6.0.4 (!1410)
- Templating FieldUtils::Interpolator class (!1420)
- Fix virtual function overrides in StdRegions and LocalRegions classes (!1435)
- Disable -Werror by default (!1443)
- Add missing override keyword to virtual functions in FieldUtils (!1452)
- Add override keyword to virtual functions in GlobalMapping and MultiRegions (!1450)
- Add fmod and modulus operator to interpreter (!1089)
- Add command line option and environment variable to disable backup field files (!1154)
- Add override keyword to virtual functions in SpatialDomains (!1448)
- Add missing override keyword to virtual functions in Collections (!1453)
- Update tutorial submodule (!1511)
- Add missing override keyword to virtual functions in SolverUtils (!1451)
- Add missing override keyword to virtual functions in LibUtilities (!1459)
- Enable ARM macOS runner, fixes for SCOTCH allocation and PETSc detection on macOS (!1462)
- Add FieldConvert module and filter to project velocity into body-fitted coordinate system (!1467)
- Add Curl operator to ExpList (!1475)
- Changed LinearAdvectionDiffusionReactionSolve arguments to be equivalent to HelmSolve (!1475)
- Added detection for advection matrices in GlobalLinSysIterative with silent switch from CG to GMRES (!1475)
- Added non-symmetric laplacian matrices via varcoeffs (!1475)
- Added matrix-free operator for LinearAdvection (!1475)
- Fix uninitialized coordinates in the Bodyforcing (!1472)
- Fix body-fitted velocity filter and also record the max/min for density,pressure, and temperature field (!1490)
- Fix typos in Vmath and VDmath (!1480)
- Fix minor typo and removed unused functions in LibUtilities/TimeIntegration (!1476) 
- Fix RK5 time integration scheme (!1482)
- Fix fld file import for SingleMode expansion (!1487)
- Fix ESDIRK time integration scheme (!1484)
- Fix IMXGear time-integration scheme for consistent second-order accuracy (!1489)
- Fix ESDIRK time integration scheme (!1484)
- Fix TimeIntegrationDemo.cpp and add ESDIRK tst files to the CI (!1485)
- Add DIRKOrder1, BDFImplicitOrder3, BDFImplicitOrder4, RungeKutta1, and RungeKutta3 schemes to the register (!1485)
- Use DIRK (instead of IMEXdirk) schemes for the start-up phase of high-order BDF and AM schemes (!1485).
- Fix IMEXdirk_1_2_2 and IMEXdirk_2_3_3 time-integration schemes (!1499)
- Add extrapolation time-integration scheme (!1488)
- Fix CNAB/MCNAB time-integration schemes (!1493)
- Slightly tidy-up time integration algorithms (!1496)
- Reduced memory usage in the FilterHistoryPoint (!1458)
- Remove redundant functor typedef (!1498)
- Add missing m_ prefix to member variables in FFTW (!1504)
- Make some virtual functions protected (!1506)
- Remove trailing CONDITIONS tag in xml files (!1510) 
- Disable problematic Movement_fixed_3D_stacked_cylinders_curved_hdf5_par test (!1507)
- Fix I/O issue related to Hdf5 that was unable to open file and fixed similar issue in other IO classes in BasicUtils (!1512)
- Remove unused function SetUpXmlDoc (!1513)
- Add new interpolation function to FieldUtils (!1514)
- Generalize the use of the space communicator (!1518)
- Add parallel-in-time feature to FieldConvert (!1520)
- Add Spectral Deferred Correction (SDC) time integration schemes (!1481)
- Redesign of Spectral Deferred Correction (SDC) algorithm (!1523)
- Modify SessionReader to read restart/exact solution files parallel-in-time (!1521)
- Fix Polylib_test.cpp (!1524)
- Update to Parareal file output (!1517)
- Add convergence criteria to Parareal driver (!1457)
- Add time metadata to tecplot output (!1525)
- Fix segmentation fault when no time integration method specified for unsteady problem (!1526)
- Set adjacent elements for m_bndcondExpansions for both CG and DG (!1491)
- Fix inconsisten treatment of 1D and 2D/3D expansions in DisContField::v_GetBoundaryToElmtMap (!1491)
- Tidy-up parallel-in-time processing in FieldConvert (!1529)

**Python**
- Add wrappers for Interpreter and Equation classes (!1329)

**CompressibleFlowSolver**
- Added Laplacian (NonSmooth) AV to the explicit Navier Stokes solver (!1372)
- Added Physical AV to the implicit Navier Stokes solver (!1372)
- Fixed Segmentation Fault when using C0 Smoother with Shock Capturing (!1394)
- The Incomplete IP method was made the default method for the IP method (!1377).
- Add additional parameters for the Isentropic Vortex equation system (!1323)
- Improve performance of the perconditioner and diffusion operator (!1393)
- Re-add the SFD test with an updated restart file (!1399)
- Improve performance of the block diagonal operator of the preconditioner (!1404)
- ExtractSurface2DCSF utility is updated to use the boost program option (!1407)
- Fix a Wuninitialized-const-reference warning (!1449)
- New implementation of the Stagnation Inflow Boundary Condition (!1478)
- Remove m_root in PreconCfs to avoid possible future conflict with parallel-in-time driver (!1515)
- Update to Parareal file output (!1517)

**CardiacEPSolver**
- Fix a shadowed loop counter variable in the benchmark filter (!1436)
- Update functions in derived classes to be consistent with the base class and add override keyword to virtual functions (!1439)
- Add dummy projection to CardiacEPSolver (!1527)

**IncNavierStokesSolver**
- Replaced depricated "sprintf" with "std::to_string" (!1406)
- Extended Reynolds Stresses filter to passive scalars (!1430)
- Added Implicit Solver (!1475)
- Fixed Taylor-Hood expansion for VCSWeakPressure (!1444)
- Fix filename in LinearisedAdvection (!1479)
- Added scalar advection terms to AdjointSolver (!1466)
- Remove member variables as funtion parameters in LinearisedAdvection solver (!1522)

**VortexWaveInteractionSolver**
- Replaced depricated "sprintf" with "std::to_string" (!1406)

**DummySolver**
- Fix CWIPI test to use DirectFull for projection of received data (!1502)

**NekMesh**
- Replace VTK pointers with VTK smart-pointers to avoid memory leaking, when
exporting in .vtu format (!1386)
- Preserve CAD face labels and save in to session file as a "NAME=" tag on the composites (!1396)
- Fix a header include which caused compilation errors on OCC versions newer than v7.4 (!1395)
- Add option to refine curves in the same manner as the line refinement functionality (!1298)
- Add refined curves and refined lines now prompt the octree to subdivide until the desired refined delta is reached (!1298)
- Fix a segmentation fault with WriteOctree due to missing 'order' parameter (!1418)
- Multi domain input/output for Nekpp and HDF5 file formats (!1323)
- Fix CADSurfOCE curvature bug where negative curvature values could be returned causing incorrect mesh spacing (!1442)
- Fix ProjectCAD bug with findAndProject where the projection was missing and variable was passed without reference (!1442)
- Fix 3d_bl_wing test case for STEP files where the wrong surfaces were selected for the BL (!1442)
- Fix error when setting BL progression to 1.0 due a division by 0 (!1455)
- Changed the BOOLPARAMETERS tag in InputMCF to allow disabling the high order
  surface optimisation with "DisableSurfaceOptimiser" (surface optimisation is
  still enabled by default) (!1455)
- Fix 3d_bl_wing test case for STEP files - updated to use an improved CAD definition for the NACA aerofoil (!1486)

**FieldConvert**
- Add vars and dirs options in the gradient module to specify fields and partial derivative directions (!1415)
- Fix range option so that it also works with hdf5 (!1414)
- Fix halfmodetofourier module with triangles (!1492)
- Fix the output field names of WSS module of FieldConvert, revert !1352 (!1528)

**Miscellaneous**
- Updated gitignore to be friendly with CLion IDE (!1405)
- Correct header section of .cpp, .hpp, and .h files (!1426)
- Linux format .cpp, .hpp, and .h files (!1432)
- Fix wsign compare warning (!1437)
- Fix some Woverloaded-virtual warning (!1439)
- Add missing override keyword to virtual functions in solvers (!1440)
- Fix some Wunused-variable (!1438)
- Fix unused parameter warnings in virtual functions (!1441)  
- Fix a Wreorder warning (!1445)
- Fix some Wimplicit-fallthrough warnings (!1446)
- Switch to using pkg-config for finding PETSc (!1454)
- Use Nektar::LibUtilities::Timer for better accuracy (!1468)
- Make some virtual functions protected (!1469)
- Extend clang-format checks to solvers, utilities, tests and templates (!1434)
- Fix documentation for exponential scheme (!1519)

**CI**
- Enable packaging for Fedora 35, removed Fedora 33/34 from package builds. (!1424)
- Add header checking for \*.cpp, \*.hpp and \*.h files to the CI (!1431)
- Enable packaging for Fedora 36. (!1429)
- Fix XML files indentation (!1428)
- Update solvers CMakeList.txt to fix some warnings detection issue (!1447)
- Remove -fpermissive from NektarCommon.cmake (!1460)
- Remove old distribution versions, added Fedora 35/36 testing to CI (!1461)
- Kill orphan Tester-g processes on Windows and remove source tree after build
(!1471)
- Fixed path issue and warning in the nektar-workbook image (!1470)

v5.2.0
------
**Library**
- Add Arm SVE backend to SIMD library (!1282)
- Added support for manifold  MatrixFree operators (2D in 3D space) (!1304)
- Put in place automatic selection of explicit operations using an opt file (!1304)
- Fixed the moving reference frame rotation (Solver Utils) (!1305)
- Revised FilterAeroForces to accout for the moving reference frame (!1305)
- Add MaxMinFields filter to record the max/min at each quadrature point and output the max/min fields. (!1256)
- Simplify the logic in the MPI pairwise trace exchange (!1307)
- Fix imaginary mode in HalfModeToFourier module (!1247)
- Added a dummy output module OutputStdOut for NekMesh utilities that don't require an output file (!1318)
- Fix compiler errors on ARCHER2 using PrgEnv-cray (!1315)
- Fix cmake SIMD enable/disable options based on architecture (!1320)
- Restrucutred the communicators to reduce direct dependence on session file communicator (!1337)
- Fixed SIMD mask test (!1324)
- Fix memory leak in Timer.cpp (!1330)
- Fix cmake CWIPI option to remove Fortran check (!1331)
- Fix excessive verbose output in GetBndElmtExpansions method (!1341)
- Timer class was updated with safety checks to avoid wrong measurements (!1347)
- Fix to adjust for warnings/errors from Monterey updated compiler (!1355)
- Update `nektar` and `nektar-env` packages to Debian Bullseye (!1356)
- Reformat code with clang-format (!1359)
- Remove unnecessary IterPerExp methods (!1366)
- Fix erronous call to FwdTrans from MR 1366 (!1374)
- Fixed avx512 back-end for SimdLib (!1333)
- Added float to scalar and avx2 back-end, disable avx512, sse2, sve (!1255)
- Change MPI initialisation to allow MPI_Init call outside Nektar++ (!1376)
- Fixed incorrect summary output for diffusion/reaction terms (!1383)

**FieldConvert**
- Add calculation of CFL number for the incompressilbe flow (!1332)
- Added conditional to select the eNearestNeighbour method for 3D interpolation (!1335)
- Fixed the output field names of WSS module of FieldConvert (!1352)
- Add VTU output using VTK library (high-order & multi-block options) (!1343)

**IncNavierStokesSolver**
- Added Boundary conditions for moving reference frame (!1305)
- Added the virtual functions overwriting the FluidInterface for moving reference frame (!1305)
- Add Gradient Jump Penalty (GJP) Stabilisation into the solver (!1290)
- Equation types are registered to the session reader (!1344)
- Added Block-Preconditioner for Full Matrix solve (!1350)
- Update to Parareal file output (!1517)

**ADRSolver:**
- Add Gradient Jump Penalty (GJP) Stabilisation into the Unsteady Advection and Unsteady Advection Diffusion solvers (!1290)

**PulseWaveSolver**
- Parallelised solver (!1337)
	
**NekMesh**
- Allow for one or more blank lines between sections in Tecplot ascii (.dat) files (!1322)
- Small bug-fix for Python API for unused configuration options (!1348)
- Fix bug in ProcessVarOpti/ElUtil for segfault on non-tri or tet meshes (!1381)

**CompressibleFlowSolver**
- Added physical AV, dilatation sensor, Ducros's and smoothing (!1180)
- Added timers around important functions using the Timer class. Timers are available by specifying IO_Timer_Level > -1 (!1347)
- Fixed bug in the calculation of the discontinuity penalty factor for the DiffusionIP implementation (!1368)

**Documentation**
- Fix images not being displayed in HTML documentation and tutorials (!1370)

**CI**
- Remove unused build options (!1360)
- Enable NEKTAR_USE_VTK across full builds and in docker image (!1358)
- Add XML linting and checking in CI pipeline (!1433)

**Packaging**
- Fix various issues with debian unstable and centos8 packaging (!1362)
- Fix missing texlive package dependency for centos packaging (!1382)

v5.1.1
------
**Library**
- Fix a boost headers incompatibility with boost-1.77 (!1297) 
- Add RungeKutta4 as an alternate name for ClassicalRungeKutta4 for time integration method (!1294)

**Python**
- Fix initialisation warning when using HDF5 (!1299)
- Fix issue with implementation of Diffusion IP (!1303)
- Split Helmholtz MatrixFree operator to improve compile times (!1292)
- Fix Boost deprecated header warnings (!1302)
- Add command lines to set starting time and starting checkpoint number of a time-dependent simulation (!1309)
- Fix an index referencing error in the Collections PhysDeriv method for Hex (!1314)

**Python**
- Updates to workbook, fix bugs in StdExpansion and SessionReader with MPI communication being recreated. (!1296)

**BuildSystem**
- Updated third party Lapack version 3.7.1 (!1312)

**CompressibleFlowSolver**
- Fix non-dimensional Sutherland law (!1253)

v5.1.0
------
**Library**
- Restructure library to use local coefficient storage down to the GlobalLinSys
  level. Removed GlobalCoeffs functionality (!963, !1145)
- Corrected the use of communicator in AssemblyMapDG and AssemblyCommDG which
  was not using GetRowComm() (!1144)
- Add interior penalty method to DG framework (!1101)
- Add an error filter for the time-evolution of the L2 and Linf errors (!1147)
- Fix successiveRHS method (!1176)
- Add cachedId in GetExpIndex and use in Fieldconvert (!1167)
- Fix bug in PreconditionerLowEnergy (!1161)
- Fix intel c compiler error in AeroFilters (!1198)
- Fix compilation errors when CWIPI interface enabled (!1207)
- Fix distance in ContainsPoint and GetLocCoords (!1200)
- Fix compiler warning of maybe-uninitialized elType in InputStar (!1217)
- Extend vectoisation to include all elements and initialise collections on first call (!1162)
- Add vectorisation of most element on basix operations (!1158)
- Add constant coefficients to matrix-free Helmholtz operator (!1284)
- Limit MPI methods based on core count (!1208)
- Split out IProduct.cpp and IProductWRTDerivBase.cpp in order to avoid long time compilations (!1228)
- Refactored time integration code using factory pattern (!1034, !1103)
- Fix WriteStream with empty Array/vector (!1233)
- Add interpolation at arbitrary point in 3DH1 (!1233)
  level. Removed GlobalCoeffs functionality (!963, !1159)
- Add interior penalty method to DG framework (!1101)
- Add an error filter for the time-evolution of the L2 and Linf errors (!1147)
- Enable global systems to be generated when using different values of variable
  coefficients (!1159)

**FieldConvert**
- Refactored time integration code using factory pattern (!1034)
- Fix to preprocessor logic for boost with Visual Studio >= 2015 (!1115)
- Fix type consistency and real comparison in SharedArray.hpp, replaced
  num_elements with size() (!1127, !1137, !1141)
- Use base MPI functions instead of the GS library in the trace exchange
  for parallel DG simulations (!1112)
- Replace PhysIntegral with Integral (!1246)
- Change the way periodic boundary conditions in parallel is setup to reduce excessive memory usage (!1235) (!1289)
- Add exponential and fractional-in-time integration schemes (!1106, !1111, !1210)
- Add nonlinear and linear system solvers (!1196)
- Add ESDIRK3 and ESDIRK4 time integration schemes (!1196)
- Add a filter to calculate mean value of solution fields (!1211)
- Fix the time dependent absorption forcing (!1254)
- Enable very high order (>100) quadrature use (!1262)
- Add rotation and improve performance of MovingReferenceFrame forcing (!1185)
- Fix BODYFORCE defined by a file (!1215, !1264)
- Add multi-level partitioning strategy for HDF5 geometry (!1209)
- Fix the URL of ccmio library (!1288)

**FieldConvert**:
- Add phifile module to compute shape functions for the SPM solver (!1065)
- Fix mean and innerProduct modules in 3DH1D cases (!1157)
- Add Python interface (!1081)
- Fix wss module with nparts option and reading of parallel xml files when the root partition is missing(!1197)
- Fix a segment error in the gradient module when the number of fields is smaller than space dimension(!1216)
- Add output of wall normal data from a single point (!1237)
- Add QCriterion for 2D flow (!1243)
- Fix to interppointsdatatofld to allow for mpi processing of large files (!1191)
- Fix the logic of C0Projection:helmsmoothing (!1220)
- Fix extract module for boundaries with periodic boundary conditions (!1277)

**IncNavierStokesSolver**:
- Add MaxMinFields filter to record the max/min at each quadrature point and output the max/min fields. (!1256)
- Fix imaginary mode in HalfModeToFourier module (!1247)

**CardiacEPSolver**
- Added additional parameter sets to Fenton-Karma model (!1119)
- Fix electrogram calculation in 1D/2D domains (!1285)

**IncNavierStokesSolver**
- Add Smoothed Profile Method (SPM) for the formulation of immersed boundaries
  (!1065)
- Add new filter AeroForcesSPM to compute aerodynamic forces in immersed
  boundaries (!1065)
- Add mask function and more baseflow parameters for the linear stability problem (!1201)
- Fix dudt in high-order pressure boundary condition (!1190)
- Add flow rate forcing with a scalar (!1026)

**CompressibleFlowSolver**
- Added the selective frequency damping support for the implicit solver (!!1267)
- Added vectorisation of the Interior Penalty method (!!223)
- Added a simplified implicit solver with naive preconditioner (!!1196)
- Add BRJ preconditioner to the implicit solver (!!1212)
- Fix implicit solver for Euler system (!!1252)
- Updated WallAdiabatic/WallViscous BC to accept time-dependent perturbations on the ghost state (!1248)

**PulseWaveSolver**
- Added viscoelasticity (!1138)
- Added empirical and power laws (!1138)
- Code tidying (!1138)

**Documentation**:
- Updated Windows source build instructions in user guide (!1152)

**Tester**
- Added test metric to check if warnings appear in output and error stream (!1225)

**NekMesh**
- Improved boundary layer splitting and output to CADfix (!938)
- Improve .geo reader and support 3D geometries with voids (!1031)
- Added r-adaptation code (!1109)
- Added Python bindings, change NekMeshUtils to NekMesh (!1149)
- Added pyramid element for the Star-CCM mesh (!1229)
- Added option to use absolute tolerance in peralign (!1225)

**BuildSystem**
- Toggle build type (!1135)
- Updated minimum required CMake version to 3.5.1 (!1152)
- Updated third party Boost version 1.71 (!1152)
- Updated third party OCE version to 0.18.3 (!1234)

v5.0.3
------
**CompressibleFlowSolver**
- Fix repeated output of u,v,w for Euler system

**FieldConvert**
- Fix the Filters output files numbering (!1251, !1261)
- Fix the Filters output files numbering (!1251)
- Fix 2D surfDistance calculation (!1263)

**NekMesh**
- Fix VTK Output for 3D meshes and support XML format (!1258)

**Documentation**
- Fix documentation to note restrictions on use of coupled solver (!1268)
**Library**
- Add robustness to the read expansions (!1239)

v5.0.2
------
**Library**
- Fix bug in StdHexExp FillMode (!1192)

**Documentation**
- Updated Documentation to include HDF5 Mesh Output (!1230)
- Removed Ubuntu Trusty (14.04) from CI and added Focal (20.04) (!1238)

**CI**
- Add Debian Bullseye to CI system (!1181)

**BuildSystem**
- Updated third party zlib version to 1.2.9 to resolve OCE source build issue (!1227)
- Adding SolverUtils as a core library that is built by default (!1240)

v5.0.1
------
**Library**
- Fix incorrect coordinate dimension used in history point filter (!1118)
- Fix compile errors with GCC 9.x (!1108)
- Correct the Energy/Enstropy integral for the 3DH1 flow (!1132)
- Added IsRealEqual method to compare real numbers with relative tolerance.
  Started using it in SharedArray and in NekMesh to fix peralign-extrude tool
  chain (!1134)
- Fix performance of GetExp(coord) by using octree lookup (!1165)
- Fix Collection unit tests (!1160)
- Fix periodic boundary conditions with HDF5 input file (!1163)
- Fix DESTDIR issues for MacPorts (!1179)
- Fix Bodyforcing and history point filter bounds issue (!1184)

**IncNavierStokesSolver**
- Change the baseflow time in the Adjoint advection (!1133)

**FieldConvert**
- Fix OutputTecplot skipping final plane in 3DH1D (!1016)
- Fix Interppoints in 3DH1D (!1140)

**NekMesh**
- Fix compile errors when using intel cc (!1114)

**Documentation**
- Fix error in compilation of developer guide (!1136)

**CI**
- Added checked conversion from double to int in SessionReader (!1113)
- Switched to Gitlab CI (!1120, !1120, !1128, !1129, !1131, !1141)
- Updated bullseye build to remove UCX components (!1203)

v5.0.0
------
**Library**
- Added in sum factorisation version for pyramid expansions and orthogonal
  expansion in pyramids (!750)
- Added detection of 'abort' file to cleanly terminate simulation early (!772)
- Significant overhaul of CMake infrastructure (!770, !804)
- Fix ThridpartyCCM options (!802)
- Fix Windows CRLF tokens in GEO reader and improve comment handling (!805)
- Use chrono in Timer (!807)
- Fix caching of FUNCTION tags that read from file and provide the same
  functionality in FUNCTIONs defined for forcings (!759)
- Transition to C++11 (!795, !847)
- Add patch to tinyxml to fix size_t vs int bug (!820, !1006)
- Add ARPACK thirdparty build capabilities (!828)
- Added native support for csv files in addititon to pts (!760, !835, !906)
- Utilize LAPACK_DIR env variable to find the native blas/lapack install (!827)
- Extend AeroForces filter to compressible flows (!815)
- Remove StdExpansion use from MultiRegion (use Expansions instead). (!831)
- Move steady state check and CFL output from solvers to SolverUtils (!832)
- Remove DG advection implementation from EquationSystem (!832)
- Simplify RawType typedefs (!840)
- Remove unused files from BasicUtils (!841)
- Remove checks for old boost versions which are no longer supported (!841)
- Refactor ParseUtils to be more consistent (!843, !896, !908)
- Added support for using the distance to a specific region (e.g. outlet) in the
  function definitions for the Absorption Forcing (!769)
- Improve performance of DisContField2D::v_ExtractTracePhys (!824)
- Fix small bug in Jacobian Energy (!857)
- fix variable name overriding in file functions (!870)
- Adds CFI CAD engine back-end (!864)
- Adds CFI Mesh IO support (!864)
- Cleanup of CAD system data structures (!864)
- Fix mac OSX on buildbots (!876)
- Fix error from (!826) (!876)
- Fix minor bug in ARPACK thirdparty build cmake (!874)
- Added in sum factorisation version for pyramid expnasions and orthogonal
  expansion in pyramids (!750)
- Adjust boost third-party compilation to account for different toolset
  choices (!886)
- Switch MeshGraph to use factory pattern and add HDF5 geometry support (!900,
  !904, !941)
- Restructure the low energy preconditioner to handle pyramidic and variable
  p expansions (!920)
- Remove requirement for modmetis, switch to SCOTCH by default (!899)
- Switch MeshGraph to use factory pattern and add HDF5 geometry support
  (!900, !904)
- Fix bug in MeshPartition.cpp which caused incorrect array access when
  WeightPartitions was used in parallel (!923)
- Removed instance count from beginning of Array storage to improve memory
  alignment (!921)
- Fix naming issue of duplicate Unit tests (!924)
- Fix warnings about missing virtual destructors in abstract classes (!932)
- Fix ability to have periodic boundary conditions that are aligned by a
  rotation rather than just a translation (!933)
- Added a coupling interface to exchange data between solvers at run time
  and a DummySolver to test the implementations (!853, !931, !950, !973, !1017)
- Fix compilation issue with newer Boost versions and clang (!940)
- If only `NEKTAR_BUILD_LIBRARY` is enabled, only libraries up to and including
  `MultiRegions` will be built by default (!945)
- Dont add doxygen documentation to the all target (!834)
- Fix missing metadata import from Hdf5 files (!971)
- Fix missing flags for periodic BC in DiffusionLDG (!985)
- Add the moving reference frame as a forcing (!987)
- Added rtree for element bounding box lookup to accelerate interpolation (!996,
  !1066)
- Fix integration weights on prisms and pyramids if not using the default
  integration rule (!998)
- Fix missing ContainsPoint in Pyramid expansion (!1000)
- Added path prefixes to find packaged Scotch (!979, !1008)
- Add HDF5 geometry format (!977)
- Combine and generalise demo code in StdRegions and LocalRegions (!993)
- Fix for error output to allow for custom error streams (!944)
- Fixed bug in ReOrientQuadFacePhysMap (!1003)
- Add NekPy Python interface (!962, !990, !989, !1004, !1014, !1061, !1070)
- Fix edge case for ThirdPartyScotch and FindScoth (!1009)
- Fix to populate m_elmtToExpId map if not already set up in GetExpIndex (!1019)
- Added flag to skip periodic BCs while filling Dirichlet BCs in
  ContField3D.cpp (!1018)
- Fix bounding box for interpolation (!1033)
- Added IMEXOrder4, RK5 and AB4 time integration schemes (!1037)
- Fix TriExp.cpp orientation bug (!1048)
- Fix XML attributes in conditions.cpp to be unordered (!1015)
- Fix issue with HDF5 mesh input in serial (!1049)
- Add estimate of filters CPU time (!1044)
- Update CompressibleFlowSolver/Examples/Test_IsentropicVortex1.xml example (!1045)
- Add error if HDG used with periodic BCs (!1071)
- Fix issues related to leading factors, arithmetic order and associativity of
  exponential operator in expression evaluator (!1066)
- Remove use of `using namespace std` in header files (!1066)
- Add error messages for use of ARPACK in serial (!1079)
- Generalise ContainsPoint routine (!1078)
- Homogenized fallthrough to fix issues with gcc 7.4.0 (!1084)

**NekMesh**:
- Add feature to read basic 2D geo files as CAD (!731)
- Add periodic boundary condition meshing in 2D (!733)
- Adjust boundary layer thickness in corners in 2D (!739)
- Add non-O BL meshing in 2D (!757)
- Add ability to compile CCIO library but tar file is not yet openly
  available whist we seek permission from Simens (!799)
- Fix issue with reading CCM files due to definition of default arrays
  rather than a vector (!797)
- Fix inverted triangles and small memory issue in surface meshing (!798)
- Update for the CAD system, more advance self-healing and analysis (!822)
- Additional curve types in GEO reader: BSpline, Circle, Ellipse (!800)
- Fix default command line argument value (!823)
- Add projection meshing module which can curve linear meshes with CAD (!826)
- XML meshes now write with provenance information, including information about
  their source, for debugging purposes (!872)
- Force 3-node loops to avoid degenerate 1-triangle faces (!875)
- Smooth BL normals in 2D when normals intersect or cause invalid macro BL
  elements (!877)
- Revert triangle code to ThirdParty library (!883)
- Fix coinciding nodes issue with very fine meshes (!883)
- Skip CFI groups of bodies and non-numbered nodes (!891)
- Add ability to space out 2D BL nodes to better fit local target Delta (!890)
- Fix automatic peralign call in 2D periodic meshing (!888)
- Fix BL splitting call from MCF (!910)
- Support CFI combined lines (!917)
- Order nodes in Gmsh output (!912)
- Fix manifold face curvature nodes (!913)
- Fix writing 1D surfaces (!930)
- Fix surface string parsin in BL splitting (!937)
- Enable use of distributed packages for triangle and TetGen (!953)
- Fix issue with MLSC after Scotch conversion (!943)
- Add support for Gmsh 4.0 mesh file format (!964)
- Fix issue with extracting 1D curved surface from 2D file (!984)
- Fix surface extraction, added regression test (!994)
- Fix 2D meshing running out of memory due to missing else (!1012)
- Add support for .msh v4.1 file input (!1054)
- Added penalty term to LDG and LDGNS, slight generalization of LDG (!1080)

**FieldConvert**:
- Add input module for Semtex field files (!777)
- Fixed interppoints module (!760)
- Fix OutputTecplot in 2DH1D (!818)
- Move StreamFunction utility to a FieldConvert module (!809)
- Allow using expansion from session file with new `--useSessionExpansion`
  command line option (!842)
- Extend wss module to compressible flows (!810)
- Allow explicitly setting bool options of FieldConvert modules as false (!811)
- Enable output to multiple files (!844)
- Allow using xml file without expansion tag in FieldConvert (!849)
- Add Lambda 2 vortex detection criteria (!882)
- Add module for modifying/adding fields from expressions (!889, !903)
- Add module for evaluating the mean of variables on the domain (!894)
- Add module for counting the total number of DOF (!948)
- Fixed wss module for compressible flows (!958)
- Made Sutherland's law non-dimensional (!972)
- Add module for removing fields from .fld files (!978)
- Fixed nparts option in FieldConvert and automated Info.xml generation (!995)
- Added if statement to fix case of 1D/2D manifold interpolation in 1D/2D space,
  added check on dimensions for interpolation, fixed seg interp (!999)
- Fixed scaling for compressed xml, fixed error printout for mesh only (!1040)
- Add field conversion from Halfmode to SingleMode (!1032)
- Fix double precision output in .dat format (!1059)
- Add phase sampling feature in FilterFieldConvert (!1068)

**IncNavierStokesSolver**
- Replace steady-state check based on difference of norms by check based on
  norm of the difference, to be consistent with the compressible solver (!832)
- Updated SVV to allow for the DGKernel extension (!851)
- Pre-calculate Time invariant portion of Womersley Solution (!814)
- Fix for independent setting of SVV in Homogeneous direction (!936)
- Write flow field based on CFL threshold (!1025)
- Fix unsteady Stokes solver (!1074)

**CompressibleFlowSolver**
- Add 3D regression tests (!567)
- Introduce forcing for quasi-1D Euler simulations (!771)
- Allow performing axi-symmetric Euler and NS simulations (!771, !866)
- Add ability to use an exponential filtering for stabilization with
  seg, quad and hex elements (!771, !862)
- Fix compressible solver with NUMMODES=1 (!868)
- Introduce equations of state to account for real gas effects (!880)
- Made Sutherland's law non-dimensional (!972)
- Modified pressure outlet BCs to allow for the reference static pressure to be
  set from the VALUE fields (!981)
- hp scaling for Laplacian AV (!1013)
- Removed smooth AV (!1072)

**AcousticSolver:**
- Added two new boundary conditions to the APE system: RiemannInvariantBC
  and WhiteNoise (!782)
- Store base flow fields in a discontinuous projection (!918)
- Enabled 1D cases (!918)
- The APE system now uses u_i, c^2 and rho as base flow fields (!918)
- Added the Linearized Euler Equations (LEE) (!918)

**ADRSolver:**
- Fix forcing from file for Poisson solver (!1029)

**APESolver:**
- APESolver was replaced with AcousticSolver (!918)

**PulseWaveSolver**
- Added two new boundary conditions: AInflow and UInflow

**CardiacEPSolver**
- Converted FentonKarma model to dimensional form and added variants (!1011)

**Documentation**:
- Added an initial developer's guide (!1001)
- Updated user guide to reflect current implementation (!1051)
- Added manpages for key solvers and utilities (!1051)

**Tester**
- Fix build with boost 1.67 (!947)
- Various change to tests to decrease test time (!1053)
- Extend to support MPI tests with multiple executables (!1085)

**Packaging:**
- Add Dockerfiles and gitlab CI configuration for automatic builds (!1021,
  !1092, !1098)

v4.4.2
------
**Library**
- Fix evaluation of points (e.g. HistoryPoints, Interpolation to pts) close to
  the interface of two elements (!836)
- Fix deadlock in Hdf5 with homogeneous expansions (!858)
- Fix a few memory leaks in polylib (!863)
- Fix a crash when Interpolator is called on an empty field (!869)
- Fix petsc compile without MPI (!873)
- Fix calculation of BLPoints (!892)
- Fix deadlock in DiffusionLDG (!885)
- Fix uninitialised coefficients in DirectFull solver (!898)
- Updated PETSc to 3.7.7 (!916)
- Fix typecast to an integer which set Lz < 1 to zero when postprocess hdf5 output (!922)
- Fix program options errors on Windows in debug mode (!986)
- Fix potential clobbered output of ModArnoldi EVs when run in parallel (!983)

**IncNavierStokesSolver**
- Add a test for imaginary shift to be only used with Homogenous and SingleMode on. (!928)

**NekMesh**
- Fix missing periodic boundary meshing and boundary layer mesh adjustment
  configurations in 2D (!859)
- Fix 2D BL splitting where out-of-plane nodes would be created (!887)

**Documentation**:
- Fix sign of the viscous term in the velocity correction scheme equations in
  the user guide (!856)
- Fixed anonymous clone URL (!909)
- Add information on the limitations of Imaginary Shift for stability. (!928)

**FieldConvert**
- Allow passing input name with trailing separator (!879)
- Fix the interpcoord option  of the interppointdatatofld module (!952)

**Utilities**
- Fix VtkToPng to account for deprecated VTK API for VTK version > 8.1 (!925)

v4.4.1
------
**Library**
- Remove m_offset_elmt_id and GetOffsetElmtId which fixed problems in 2D when
  quad elements are listed before tri elements (!758)
- Remove the duplicate output of errorutil (!756)
- Fix BLAS CMake dependencies (!763)
- Fix interpolation issue with Lagrange basis functions (!768)
- Fix issue with average fields not working with different polynomial order
  fields (!776)
- Fix rounding of integer parameters (!774)
- Fix Hdf5 output in FilterFieldConvert (!781)
- Fixed extreme memory consumption of Interpolator when interpolating from pts
  to fld or between different meshes (!783)
- Fix deadlock with HDF5 input (!786)
- Fix missing entriess in LibUtilities::kPointsTypeStr (!792)
- Fix compiler warnings with CommDataType (!793)
- Fix ability to set default implementation in Collections and added an option
  to set eNoCollections in FieldConvert as default (!789)
- Fix performance issue in ProcessIsoContour in relation to memory consumption
  (!821)
- Fix performance issue with ExtractPhysToBndElmt (!796)
- Fix available classes being listed multiple times (!817)
- Fix Intel compiler warnings (!837)
- Fix overwriting and backup of chk/fld files on slow file systes (!741)
- Fix DriverAdaptive with second order IMEX (!850)
- Fixed typo in eIMEXGear part (!854)
- Added regression tests for IMEXOrder1, IMEXOrder2, IMEXOrder3, MCNAB,
  IMEXGear, CNAB, 2nd order IMEX-DIRK, 3rd order IMEX-DIRK (!854)
- Fix bug due to subtractive cancellation in polylib routines (!778)

**FieldConvert:**
- Fix issue with field ordering in the interppointdatatofld module (!754)
- Fix issue with FieldConvert when range flag used (!761)
- Fix issue when using output-points combined with noequispaced (!775)
- Fix equispacedoutput for 3DH1D with triangles (!787)

**NekMesh**:
- Fix memory consumption issue with Gmsh output (!747, !762)
- Rework meshing control so that if possible viewable meshes will be dumped
  when some part of the system fails (!756)
- Add manifold meshing option (!756)
- Fix issue with older rea input files (!765)
- Fix memory leak in variational optimiser, add small optimisations (!785)
- Check the dimensionality of the CAD system before running the 2D generator (!780)
- Fix uninitialised memory bug in Nek5000 input module (!801)

**IncNavierStokesSolver**
- Fix an initialisation issue when using an additional advective field (!779)
- Fix MovingBody boundary condition (!852)

**Utilities**
- Fix vtkToFld missing dependency which prevented compiling with VTK 7.1 (!808)

**Documentation**
- Added missing details on artificial viscosity and dealising to compressible
  flow solver user guide (!846)

**Packaging**
- Added missing package for FieldUtils library (!755)

**ADRSolver:**
- Fix UnsteadyAdvectionDiffusion with DG (!855)

v4.4.0
------
**Library**:
- Add support for variable polynomial order for 3D simulations with continuous
  Galerkin discretisation (!604)
- Bump version of gsmpi to suppress autotuning output unless `--verbose` is
  specified (!652)
- Add support for variable polynomial order with periodic boundary conditions
  (!658)
- Statistics are now printed for lowest level of multi-level static condensation
  (!656)
- Sped up interpolataion from pts files and fixed parallel pts import (!584)
- Increased required boost version to 1.56.0 (!584)
- New FieldUtils library allows support for most `FieldConvert` post-processing
  operations during simulation using a new filter (!589)
- Adjust CMake dependencies to reduce compile time (!671)
- Homogeneous1D dealiasing improvements (!622)
- Add support for HDF5 as an alternative output to XML-based output, including
  refactoring of FieldIO, improvements to MPI interface and added communicators
  to boundary conditions (!615)
- Allow expansions to be loaded directly from field file (!617)
- New options for load balancing (DOF or BOUNDARY) in mesh partitioner (!617)
- Rework nodal utilities to support nodal prismatic elements (!660)
- Update Body/Field forces at each timestep (!665)
- Update nodalutil to include quad and hex elements and introduce SPI nodal
  points (!696)
- Add ability to restart time-averaging and Reynolds stresses from checkpoint
  file (!678)
- Extend ExtractDataToCoeffs to support interpolation between basis types for
  quads and hexahedra (!682)
- Enabled MUMPS support in PETSc if a Fortran compiler was found and added 3D
  support to the Helmholtz smoother used e.g. in FieldConverts C0Projection
  module (!714)
- Fix bug in `Vmath::FillWhiteNoise` which caused `ForcingNoise` to have
  a repeated pattern (!718)
- Fix bug in the calculation of the RHS magnitude in CG solver (!721)
- Fix bug in MPI detection for recent CMake on OS X (!725)
- Fix bug in CMake Homebrew and MacPorts detection for OS X (!729)
- Fix bug in FieldUtils when using half mode expansions (!734)
- Do not read the same fld/pts files again for every variable (!670)
- Fix bug in CMake PETSc detection for Ubuntu 16.04/Debian 9 (!735)
- Fix warnings with Intel compiler (!742)

**ADRSolver:**
- Add a projection equation system for C^0 projections (!675)

**APESolver:**
- Use a continuous basefield projection and revert to constant c formulation (!664)
- Added ability to compute CFL number (!664)
- Output Sourceterm (!664)
- Use the Forcing framework to define source terms (!665)

**IncNavierStokesSolver:**
- Add ability to simulate additional scalar fields (!624)
- Improve performance when using homogeneous dealiasing (!622)
- Fix linearised advection for full 3D cases (!708)
- Added a weak pressure formulation following Guermond & Shen (!713)
- Added a convective like outflow boundary condition from Dong (!713)
- Added the ability to specifiy Womersley boundary conditions for pulsatile flow (!472)

**CardiacEPSolver:**
- Added a Python translator utility to generate cell models from CellML (!723)

**FieldConvert:**
- Allow equi-spaced output for 1D and 2DH1D fields (!613)
- Update quality metric to include scaled Jacobian output (!695)
- Allow multiple XML files to be specified in InterpField module (!705)
- Fix issues with isocontour module (!719)
- Fix issue with interpolator routine (!746)

**NekMesh:**
- Modify curve module to allow for spline input (!628)
- Add STL surface writer module (!668)
- New module for inserting an alternate high-order surface into the working
  mesh (!669)
- Add curve projection routines to CAD system (!697)
- Extensive clean-up of NekMeshUtils/MeshElements and extension of makeorder to
  consider CAD information (!698)
- Improvements to mesh linearisation module (!659)
- Add support for Gmsh high-order output (!679)
- Move CAD classes to factory format (!676)
- Add module to check topology of the mesh along with boundary connectivity
  to detect problems such as hanging nodes (!691)
- Add option to `linearise` module to linearise only prisms (!688)
- Add reader for Nek5000 mesh files (!680)
- Add option to `linearise` to use element quality (!690)
- Add flag to `insertsurface` process for non-conforming geometries (!700)
- Bug fix to get two meshgen regression tests working (!700)
- Remove libANN in deference to boost::geometry (!703)
- Refactor library to use NekMesh modules for CAD generation (!704)
- Add `varopti` process module to optimise meshes (!711)
- Add a mesh extract option to the linearise module to visualise the result
  (!712)
- 2D to 3D mesh extrusion module (!715)
- Add new two-dimensional mesher from NACA code or step file (!720)
- Add basic gmsh cad (.geo) reader to the meshing system (!731)
- Fix inverted boundary layer in 2D (!736)
- More sensible element sizing with boundary layers in 2D (!736)
- Change variable names in mcf file to make more sense (!736)
- Fix issues in varopti module so that in can be compiled without meshgen on
  (!736)
- Replace LAPACK Eigenvalue calculation with handwritten function in
  varopti (!738)
- Improved node-colouring algorithm for better load-balancing
  in varopti (!738)
- Simplified calculation of the energy functional in varopti for improved
  performance (!738)

**FieldConvert:**
- Move all modules to a new library, FieldUtils, to support post-processing
  during simulations (!589)
- Add module to stretch homogeneous direction (!609)
- Add module to add composite ID of elements as a field (!674)
- Add reader for Nek5000 field files (!680)

**Tester:**
- Fix output not displayed on segfault or system error (!745)

v4.3.5
------
**Library:**
- Fix bug in DG with hybrid meshes (!694)
- Fix issue with parallel output (!699)
- Fix performance issue with iterative full solver (!693)
- Enforced precision on history point output (!706)

**Documentation**
- Update build instructions in user guide for Windows (!692)

**Tester**
- Fix bug in tester when no parameters specified for test executable (!701)

v4.3.4
------
**Library:**
- Fix performance issue with `v_ExtractDataToCoeffs` for post-processing of
  large simulations (!672)
- Added additional assertions to ensure homogeneous simulations have an even
  number of planes per process (!666)
- Fix compilation with NEKTAR_USE_MESHGEN option
- Fix IterativeFull solver in parallel (!685)
- Fix error message for missing fld file (!689)

**IncNavierStokesSolver:**
- Fix 2nd order time-integration for VCSMapping (!687)

v4.3.4
------
**Library:**
- Fix performance issue with `v_ExtractDataToCoeffs` for post-processing of large
  simulations (!672)

v4.3.3
------
**Library**:
- Auto-detect a shared filesystem and removed --shared-filesystem option (!654)
- Fix filters when using adaptive driver to avoid output being overwritten after
  each adaptive update (!588)
- Minor fix to suppress Xxt output unless `--verbose` is specified (!642)
- Fix of DirectFull solver in case where only Neumann boundary conditions
  are imposed. (!655)

**FieldConvert**:
- Fix to avoid repeated import of field file (!649)
- Fix issue with C^0 projection (!644)
- Fix verbose output when using --procid (!648)

**NekMesh:**
- Fix namespace issue in Star-CCM+ input header in NekMesh (!661)

**CompressibleFlowSolver**:
- Fix issue with residual output (!647)
- Issues with 1D Euler solver fixed (!565)
- Fix deadlocking issue with boundary conditions (!657)

**Packaging**:
- Fix NekMesh dependencies for DEB package (!650)
- Fix PETSc build on newer linux distributions (!646)

v4.3.2
------
**Library**:
- Add small optimisation for DriverAdaptive (!618)
- Updated FFTW build to use the compiler used for building Nektar++ (!629)
- Fix numbering bug in periodic boundary conditions (!631)
- Print error message for invalid equation also in release version (!634)
- HistoryPoints filter now uses closest plane to requested z-coordinate and
  output is produced in physical space (!621).
- Fix minor performance issue with time integration schemes (!632)
- Fix FilterCheckpoint filter to be consistent with `IO_CheckSteps` (!633)
- Fix CMake configuration for building on Windows 10 with VS 2015 (!641)
- Fix `IO_CheckSteps` to avoid missing first checkpoint (!639)
- Fix bug in iterative solver where only root process would ASSERT when
  exceeding the maximum number of iterations (!636)

**FieldConvert**:
- Fix appearence of duplicate messages when running in parallel (!626)
- Fix issue with efficiency when using large number of 3DH1D planes (!627)
- Add module for combining average fields (!620)
- Fix wall shear stress processing module for parallel execution (!635)

**Packaging**:
- Fixes for DEB package dependencies (!630)

v4.3.1
------
**Library**:
- Add `THIRDPARTY_USE_SSL` option to disable use of SSL on systems where CMake
  is not compiled with SSL support. (!602)
- Fixed a number of documentation issues (!586, !593, !596)
- Fix Homogeneous transform when unshuffling is not used. (!599)
- Fix namespace pollution in library header files. (!601)
- Fix issue with METIS compilation on clang 7.3 (!603)
- Fix issue with heterogeneous quadrilaterals (!607)
- Fix bug in modified Arnoldi algorithm causing convergence to be reported when
  number of vectors is less than `nvec` (!608)
- Fix uninitialised array bug in AssemblyMap (!598)
- Fix issue with LAPACK call in eigenvalue calculation (!610)
- Fix FieldConvert processing of partitions in serial (!612)
- Fix use of multi-level static condensation in parallel with periodic
  boundary conditions (!614)
- Fix NaN detection to work in parallel (!605)
- Add additional constructor to ContField3DHomogeneous1D for FieldConvert
  extract module. (!590)

**NekMesh**:
- Fix incorrect link directory on CCMIO library.

**FieldConvert**:
- Fix to FLD input to update the field definitions always, not just when a range
  is specified. (!611)

**Tester**:
- Remove requirement for executable to be specified in .tst file if it is
  overridden on the command-line (!595)

**Packaging**:
- Fix dependency resolution on generation of DEB packages. (!616)

v4.3.0
------
**Library:**
- Changed default XML format to compress mesh data (!533, !547)
- Various fixes for 3D homogeneous post-processing (!531, !529, !528, !526, !521)
- Fix boundary condition imposition for 3D homogeneous 2D HelmSolve (!545)
- Fix range with variable p option (!522)
- Fix bug with hexahedra of heterogeneous order (!520) and reading files (!522)
- Fix history point output formatting (!518)
- Fix for OS X 10.11 (!512)
- Fix `HexGeom::v_GetDir` to support heterogeneous basis functions (!520)
- Added new `NekMeshUtils` library to support new `NekMesh` executable and
  associated CAD routines. Old CAD wrappers in LibUtilities now moved to
  `NekMeshUtils` (!527)
- Fix initialisation bug in ExpList2DH1D and ExpListHomogeneous2D (!528, !529)
- Fix bug in ExpList1D which may lead to invalid .vtu files (!531)
- Make `GetBoundaryToElmtMap` consistent for 3DH1D (!526)
- Add support for PETSc matrix shell to use Nektar++ operations/preconditioners
  (!537)
- Fix bug with initial conditions of CG simulations using variable P (!543)
- Fix bug in 3DH2D with non-zero Dirichlet boundary conditions (!545)
- Added in a method to convert equispaced interpolated points back to
  coefficients which requires the introduction of a new StdRegions matrix.(!561)
- Empty XML tags which would override non-empty XML tags are now ignored (!581)
- Add contribution guide (!551)
- Add a filter to calculate exponential moving averages (!566)

**APESolver:**
- Fix restarting from checkpoint file (!517)

**IncNavierStokesSolver**
- Fix floquet stability analysis for HalfMode case (!536)
- Add a filter to calculate Reynolds stresses (!566)

**FieldConvert:**
- Extended surface distance module to support hexahedra and quads (!524)
- Small fixes in interpolation routine (!515)
- Add support for surface extraction in 3DH1D case (!521)
- Add support for isocontour extraction for 3DH1D (!525)
- Add process module to calculate high-order mesh quality metric (!527).
- Add module to extract one of the planes of 3DH1D (!542)
- Add module to enable mean mode of 3DH1D to be extracted (!530)
- Fix bug in C^0 projection (!541))
- Add command line option to set number of homogeneous planes (!540)
- Add module to project set of points to a fld file(!561)
- Add support for interpolating to a box of points and fix ability to run
  interppointstofld module in parallel when using a plane or box option (!561)
- Add option to output equi-spaced points in VTU format (!550)
- Add module innerproduct (!568)
- Add command line option of `--part-only` and `--part-only-overlapping` (!569)

**NekMesh:**
- `MeshConvert` is now renamed to `NekMesh` to reflect new mesh generation
  functionality (!527).
- Enable face curvature inside core MeshConvert objects (!511)
- Add linearise processing module to remove all curvature from high order
  elements (!509)

**Documentation:**
- Added git submodule for including Nektar++ tutorials in the source tree (!507)

v4.2.0
------
**Library:**
- Add Runge-Kutta SSP schemes for 2nd/3rd order using keys `RungeKutta2_SSP` and
  `RungeKutta3_SSP`. `ClassicalRungeKutta4` is now called `RungeKutta4`. (!481)
- Add rudimentary support for 3D CAD models using OpenCascade - work in progress
  (!486)
- Allow filters to evaluate expressions in their parameter definitions (!489)
- Fix block preconditioner to work with periodic boundary conditions (!420)
- Dump a backtrace when crash occurs and Nektar++ is compiled in FullDebug mode
  (!495)
- Stop the execution of a time-dependent solver if NaN is detected in the
  solution field (!496)
- Fixes to improve robustness of interpolation routines (!499)
- Allow solvers to use multi-level static condensation with Xxt, most useful
  when running a 3DH1D simulation (!502)

**IncNavierStokesSolver:**
- A range of fixes for the coupled stability solver, which now works in parallel
  (!508)

**MeshConvert:**
- Add module to extract prismatic boundary layer elements from mixed prism-tet
  mesh (!493).

**FieldConvert:**
- Add a processing module to calculate height of an element connected to a
  surface, allowing for calculation of y plus values (!488)
- Fixes for equispaced output (!510)

v4.1.0
------
**Library:**
- Add support for interpolating point data from .pts files (!433)
- Fixes for curvilinear element normals (!443)
- Fix consistency issues between FFT and MVM approaches for homogeneous
  expansions (!444)
- Fix a bug in Tecplot output (!445)
- Fix a bug with PETSc and MPI_Finalize (!456)
- Fix bugs with mesh partitioning (!449, !480)
- Fix a bug with non-symmetric SVV parameters for curvilinear elements (!451)
- Fix detection of Intel MKL 2013/2015 (453)
- Fix linearised stability solver in parallel (!454)
- Add a filter for 1D energy spectra (!457)
- Add an incomplete developer guide containing most information from the wiki
  (!459)
- Change user defined boundary conditions to remove dependency on enumerator
  inside SpatialDomains (!460)
- Add a new collections library for optimised evaluation of operators (!461)
- Change minimum version of boost to 1.52.
- Add initial multithreading support (!463)
- Fix third-party boost compilation on OS X (!467)
- Disable some regression tests on 32-bit systems (!468)
- Fix memory issues inside collections (!473)
- Fix collections autotuning (!476)
- Fix VtkToPng utility (!477)
- Add PulseWaveSolver to packaging (!478)
- Fix bug in iterative static condensation solver (!483)
- Fix zlib install path on OS X (!484)
- Fix documentation HTML styling for user and developer guide (!485)
- Add fixes to support native Nektar++ extension in VisIt visulisation software
  (!490)
- Fix warnings on OS X (!491)

**CardiacEPSolver:**
- Fixes for stimuli (!442, !446), conductivity (!441), cell restarts (!458)
- Add a new filter for outputting cell states at specific points over time (!465)

**Linear elastic solver (new):**
- Add solver for linear elasticity equations (!400)

**IncNavierStokesSolver:**
- Add support for moving bodies (!344, !448)
- Fixes for modal energy filter (!427)
- Fix import of mesh file in the Adaptive SFD driver (!440) and other general
  fixes (!452)
- Documentation for high order pressure and outflow boundary conditions (!447)
- Update examples to use correct forcing terms (!470)
- Fixes for half-mode stability (!471)
- Fix static initialisation problem in extrapolation classes (!492)

**CompressibleFlowSolver:**
- Add support for sponge region (!396)
- Add support for adiabiatic walls (!430)
- Add utility to generate boundary layer from similarity solution (!438)

**ShallowWaterSolver:**
- Added a DG solver for the Boussinesq equations of Peregrine (!431)

**APESolver:**
- Add support for variable speed of sound (!438)

**MeshConvert:**
- Fix Star file input for highly stretched elements (!455)
- Add Star input from binary format (!474)
- Tidy up files to align with FieldConvert (!479)

**FieldConvert:**
- Major re-organisation of modules, most post-processing utilities now available
  within FieldConvert (!475)

v4.0.1
------
**Library:**
- Change hybrid parallelisation to use command line options (!368)
- Add support for multi-variable functions in expression evaluator: new
  functions include rad and ang for polar coordinates (!375)
- Add more documentation (!376, !383)
- Various OS X (!377, !378, !382, !425), compiler warning (!432), documentation
  (!434) Windows 7 (!391, !407), CMake (!392, !415), packaging (!435, !436) and
  Intel compiler (!414, !416) fixes
- Refactor of CG and DG assembly maps (!380)
- Fixes for PETSc running in serial (!381, !420)
- Fixes for running Arnoldi solver in parallel (!384)
- Enable MPI tests on Cray machines such as ARCHER (!386)
- Fix issues with extracting face physical values (!393)
- Fix threshold filter (!395)
- HDG can now use block preconditioner (!397)
- Fix issue with singular vertices in parallel (!398)
- Timing executables now use `Timer` class from LibUtilities (!402)
- Fix manifold history points again (!410)
- Fix time output inside energy filter (!412)
- Fix GetExpIndex function (!417)
- Fixes to external project compilation (!419)
- Fixes from CPC paper review (!422)
- Fixes for scotch partitioner tests (!423)
- Fixes for ACML BLAS libraries (!424)
- Allow prepartitioned meshes to be used (!426)
- Enable variable names to be remapped inside files to different names in XML
  functions (!428)

**APESolver:**
- Fixes for tests (!404)
- Add support for advection classes (!408)

**CardiacEPSolver:**
- Add benchmark (!411)
- Fix cardiac exmplaes (!418)

**CompressibleFlowSolver:**
- Add filter for kinetic energy/enstrophy calculation (!388)

**FieldConvert:**
- Support equi-spaced output for simplex elements to reduce storage (!421)

**IncNavierStokesSolver:**
- Unify advection classes with those in `SolverUtils` (!403, !408)

**MeshConvert:**
- Boundary layer refinement now supports hexahedra (!390)
- Improve support for Gmsh high order elements (!401)
- Many fixes for face-interior curvature (!401)
- Add rudimentary test suite (!401)
- New module for imposing curvature based on a scalar function (!401)

v4.0.0
------
**Library:**
- Update boost to 1.55 (!289)
- Fix parallel history points on manifold (!298)
- Add support for scotch partitioner (!311)
- Fixes for thirdparty builds (!319, !330, !353)
- Fix CMake >= 3.0.0 warnings (!320)
- Add support for PETSc library and tidy up global system classes (!322)
- Fixes for 1D Helmholtz solver (!326)
- Fixes for history points (!327) and solver output (!331)
- Fix issue with mesh IDs that do not start from zero (!354)

**CardiacEPSolver:**
- Simplify support for global conductiity (!295)

**FieldConvert:**
- Fixes for parallel operation and interpolation of points (!351)

**IncNavierStokesSolver:**
- Fixes for sponge layer (!272)
- Fix setting of initial conditions (!298)

v3.4.0
------
**Library:**
- New parallel output format. Parallel files are now stored in directories which
  contain partition information. (!100, !102, !236, !242, !249, !256).
- gzip-compressed XML mesh files are now supported with extension .xml.gz (!116,
  !140, !186).
- HDG solvers now run in parallel and have post-processing utilities (!188,
  !230).
- Partitioning can be done only on root process if shared filesystem is
  present with use of `--shared-filesystem` command line option (!220, !250).
- A variety of preconditioners are now supported, including linear space and
  low-energy preconditioning (!148).
- Many changes to geometric factors storage and interpolation (!99, !197).
- Improvements to identification of invalid elements (!208, !227).
- Removed elemental storage to reduce memory consumption by 30-50% for large
  problems (!240).
- Various performance and design improvements for discontinuous formulation (!134).
- Periodic boundary conditions are supported in 3D for both continuous and
  discontinuous formulations (!139, !150, !152, !155, !196).
- Utilities added to mesh converter to help identify pairs of periodic faces
  (!214).
- Preconditioner support for periodic boundary conditions (!231, !239).
- New radiation boundary condition type (!74).
- Some solvers (compressible flow solver, advection-diffusion-reaction solver)
  now support dealiasing options (!78, !146, !167).
- BLAS and vectorisation performance improvements for static-condensed iterative
  solver (!86, !109).
- New driver to improve steady state convergence and add parallel support (!91,
  !235).
- Updated to METIS v5.1.0 (!97, !142, !189).
- Iterative solvers now use previous timestep (when available) to improve
  convergence speed (!106).
- Added CPU timing for timestep loop (!156).
- Added provenance information (date, time, code version, git revision, etc) to
  field file output (!179).
- Disabled long-running regression tests by default (!183).
- Support for command line arguments without parameters (!187).
- Added support for reading boundary conditions from files, and appropriate
  utilities in MeshConvert to extract surfaces (!226).
- Updated XXt and Gs libraries to latest version (!232).
- Fix singularity check for Poisson equations (!74, !154).
- Fixes for 2D Gauss points (!73, !149, !157).
- Fixes to parallel I/O (!77, !218, !264).
- Fixes for parallel implementation (!93, !107, !121, !169, !217, !245, !246).
- Fixes for normal calculation (!94, !135).
- Improved compilation techniques, particularly when compiler includes MPI
  automatically (!80, !82, !84, !85, !113, !114, !131, !141, !166, !210, !241).
- Updated zlib to v1.2.7 (!115).
- Fix for boost 1.5.3 compilation (!120).
- Most compiler warnings silenced with clang/gcc (!81, !92, !103, !123, !201,
  !243).
- Attempts to improve mesh partitioning/load balancing (!160, !170, !175).
- Fixes for Newton iteration to interpolate inside deformed elements (!216,
  !251).
- Fixed curved tetrahedron and hexahedron issue (!219, !248).
- Fixed reading of field files for tetrahedron (!228).
- Fixed uninitialised variable inside SessionrReader (!233).
- Various improvements to support use of Nektar++ externally (!111, !260, !261).
- Fixed base flow reading (!112).

**CardiacEPSolver:**
- Cardiac electrophysiology solver improvements (!87, !95, !96, !108, !119,
  !165, !173, !174, !199, !222).

**CompressibleFlowSolver:**
- Compressible Navier-Stokes equations are now available for both DG and FR
  discretisations (!110, !125, !128).
- Meshes with spatially varying p in both 2D and 3D are now supported (!158).
- Homogeneous Fourier extension is now supported (!180).
- Various fixes (!90, !98, !147, !172).

**DiffusionSolver (new):**
- Added small solver to demonstrate usage of higher library levels outside of
  EquationSystem (!225).

**IncNavierStokesSolver:**
- Major refactoring of time-integration classes (!181, !184).
- Summary information now generated via callbacks (!182).
- Implemented new generic forcing function classes (!194).
- Current time now written out in field files (!198).
- Major refactoring of incompressible Navier-Stokes solver to improve
  readability and performance (212, !213).
- Spectral vanishing viscosity for stabilisation (!101, !104, !211, !263).
- Added filter to compute aerodynamic forces on surfaces (!168, !203, !204).
- Added filter to compute kinetic energy and enstrophy (!207, !257).

**ShallowWaterSolver:**
- Various improvements/modernisations to shallow water solver (!190).

**Utilities:**
- VTK to PNG converter (!122)
- Added scalar gradient utility (!129, !252).
- Added utility to calculate Q-criterion field (!153).
- Added support to XmlToVtk to write Jacobian field (!223).
- Added utility to calculate wall shear stress (!224).
- Fixed vorticity calculator (!138).

**MeshConvert:**
- Added face-interior quadrature and 2D/3D manifold support to spherigon code
  (!130).
- Fixes for boundary layer refinement and prism-to-tetrahedron splitting (!137,
  !202, !206, !244).

**FieldConvert (new):**
- Added new FieldConvert utility which will eventually encompass most existing
  utilities (!255).
