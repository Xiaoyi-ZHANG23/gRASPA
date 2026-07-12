#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>


static inline void WriteBox_LAMMPS(Atoms* System, Components& SystemComponents, ForceField& FF, Boxsize& Box, std::ofstream& textrestartFile, std::vector<std::string>& AtomNames)
{
  size_t NumberOfAtoms = 0;
  size_t printableAtoms = 0;
  // Count both total atoms and printable atoms
  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    // Original count: all atoms
    NumberOfAtoms += SystemComponents.NumberOfMolecule_for_Component[i] * SystemComponents.Moleculesize[i];
    
    // Count printable atoms (with print flag set to true)
    Atoms Data = System[i];
    size_t componentPrintableAtoms = 0;
    for(size_t j = 0; j < Data.size; j++)
    {
      if(Data.Type[j] < SystemComponents.PseudoAtoms.print.size() && SystemComponents.PseudoAtoms.print[Data.Type[j]])
      {
        componentPrintableAtoms++;
      }
    }
    printableAtoms += componentPrintableAtoms;
    fprintf(SystemComponents.OUTPUT, "Printing: Component: %zu [ %s ], NumMol: %zu, Molsize: %zu, Total atoms: %zu, Printable atoms: %zu\n", i, SystemComponents.MoleculeName[i].c_str(), SystemComponents.NumberOfMolecule_for_Component[i], SystemComponents.Moleculesize[i], SystemComponents.NumberOfMolecule_for_Component[i] * SystemComponents.Moleculesize[i], componentPrintableAtoms);
  }
  fprintf(SystemComponents.OUTPUT, "Total atoms: %zu, Printable atoms: %zu\n", NumberOfAtoms, printableAtoms);
  textrestartFile << "# LAMMPS data file written by WriteLammpsdata function in RASPA(written by Zhao Li)" << '\n';
  textrestartFile << "# Total atoms: " << NumberOfAtoms << ", Printable atoms: " << printableAtoms << '\n';
  textrestartFile << printableAtoms << " atoms" << '\n';
  textrestartFile << 0 << " bonds" << '\n';
  textrestartFile << 0 << " angles" << '\n';
  textrestartFile << 0 << " dihedrals" << '\n';
  textrestartFile << 0 << " impropers" << '\n';
  textrestartFile << 0 << " " << Box.Cell[0] << " xlo xhi" << '\n';
  textrestartFile << 0 << " " << Box.Cell[4] << " ylo yhi" << '\n';
  textrestartFile << 0 << " " << Box.Cell[8] << " zlo zhi" << '\n';
  textrestartFile << Box.Cell[3] << " " << Box.Cell[6] << " " << Box.Cell[7] << " xy xz yz" << '\n' << '\n';
  textrestartFile << FF.size << " atom types" << '\n'; 
  textrestartFile << 0 << " bond types" << '\n';
  textrestartFile << 0 << " angle types" << '\n';
  textrestartFile << 0 << " dihedral types" << '\n';
  textrestartFile << 0 << " improper types" << '\n' << '\n';
  //textrestartFile << std::puts("%zu bonds\n", NumberOfAtoms);
  textrestartFile << "Masses" << '\n' << '\n';
  for(size_t i = 0; i < FF.size; i++)
  { 
    double mass = SystemComponents.PseudoAtoms.mass[i];
    textrestartFile << i+1 << " " << mass << " # " << AtomNames[i] << '\n';
  }
  textrestartFile << '\n' << "Pair Coeffs" << '\n' << '\n'; 
  for(size_t i = 0; i < FF.size; i++)
  { 
    const size_t row = i*FF.size+i;
    textrestartFile << i+1 << " " << FF.epsilon[row]*0.00239006 << " " << FF.sigma[row] << " # " << AtomNames[i] << '\n';
  }
}

static inline void WriteAtoms_LAMMPS(Atoms* System, Components& SystemComponents, Boxsize& Box, std::ofstream& textrestartFile, std::vector<std::string>& AtomNames)
{
  textrestartFile << '\n' << "Atoms" << '\n' << '\n';
  size_t Atomcount=0; size_t Molcount=0;
  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    Atoms Data = System[i];
    fprintf(SystemComponents.OUTPUT, "Component %zu, Molsize: %zu\n", i, Data.Molsize);
    double3 Wrap_shift;
    for(size_t j = 0; j < Data.size; j++)
    {
      //Wrap To Box//
      double3 pos = Data.pos[j];
      if(j % Data.Molsize == 0) 
      {
        WrapInBox(pos, Box.Cell, Box.InverseCell, Box.Cubic);
        Wrap_shift = pos - Data.pos[j];
      }
      else
      {
        pos += Wrap_shift;
      }
      // Check print flag before writing atom
      if(Data.Type[j] >= SystemComponents.PseudoAtoms.print.size() || !SystemComponents.PseudoAtoms.print[Data.Type[j]])
      {
        continue; // Skip this atom if print flag is false
      }

      int oldPrecision = std::cout.precision(); // Store current precision
      int newPrecision = 10;
      std::cout.precision(oldPrecision); // Reset to previous precision

      textrestartFile << Atomcount+1 << " " << Molcount+Data.MolID[j]+1 << " " << Data.Type[j]+1 << " " ;
      textrestartFile << std::setprecision(newPrecision) << Data.charge[j] << " ";
      textrestartFile << std::setprecision(newPrecision) << pos.x << " ";
      textrestartFile << std::setprecision(newPrecision) << pos.y << " ";
      textrestartFile << std::setprecision(newPrecision) << pos.z;
      std::cout.precision(oldPrecision);
      textrestartFile << " # " << SystemComponents.MoleculeName[i] << " " << AtomNames[Data.Type[j]] << '\n';
      Atomcount++;
    }
    Molcount+=SystemComponents.NumberOfMolecule_for_Component[i];
  }
}

static inline void WriteAtoms_Restart(Atoms* System, Components& SystemComponents, std::ofstream& textrestartFile, std::vector<std::string>& AtomNames)
{
  textrestartFile << "Reactions: 0" << "\n";
  //size_t Atomcount=0; 
  size_t prevMol = 0;
  for(size_t i = SystemComponents.NComponents.y; i < SystemComponents.NComponents.x; i++)
  {
    Atoms Data = System[i];
    textrestartFile << '\n' << "Component: " <<  i - SystemComponents.NComponents.y << "   Adsorbate " << SystemComponents.NumberOfMolecule_for_Component[i] << " molecules of " << SystemComponents.MoleculeName[i] << '\n';
    textrestartFile << "------------------------------------------------------------------------" << "\n";
    size_t molsize = SystemComponents.Moleculesize[i];
    //First write positions//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-position:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << std::setprecision (15) << Data.pos[j].x << "  " << std::setprecision (15) << Data.pos[j].y << "  " << std::setprecision (15) << Data.pos[j].z << '\n';
    //Then write velocities (Zhao's note: Not implemented yet)//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-velocity:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << 0.0 << "  " << 0.0 << "  " << 0.0 << '\n';
    //Then write force (Zhao's note: Not implemented yet)//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-force:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << 0.0 << "  " << 0.0 << "  " << 0.0 << '\n';
    //Then write charge//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-charge:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << Data.charge[j] << '\n';
    //Then write scale//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-scaling:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << Data.scale[j] << '\n';
    //Finally write fixed (Zhao's note: Not implemented yet)//
    for(size_t j = 0; j < Data.size; j++)
      textrestartFile << "Adsorbate-atom-fixed:" << " " << prevMol + Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << "0 0 0" << '\n';
    prevMol += SystemComponents.NumberOfMolecule_for_Component[i];
  }
}

static inline void WriteComponent_Restart(Atoms* System, Components& SystemComponents, std::ofstream& textrestartFile, Boxsize& Box)
{
  textrestartFile << "Components: " << SystemComponents.NComponents.x - SystemComponents.NumberOfFrameworks << " (Adsorbates " << SystemComponents.TotalNumberOfMolecules-SystemComponents.NumberOfFrameworks << ", Cations 0)" << "\n";
  textrestartFile << "========================================================================\n"; 
  for(size_t i = SystemComponents.NumberOfFrameworks; i < SystemComponents.NComponents.x; i++)
  {
    textrestartFile << "Components 0 (" << SystemComponents.MoleculeName[i] << ") " << "\n";
    int fracID = -1;
    if(SystemComponents.hasfractionalMolecule[i])
    {
      fracID = SystemComponents.Lambda[i].FractionalMoleculeID;
      textrestartFile << "Fractional-molecule-id component " << i-SystemComponents.NumberOfFrameworks << ": " << fracID << "\n";
      textrestartFile << "Lambda-factors component " << i-SystemComponents.NumberOfFrameworks << ": " << SystemComponents.Lambda[i].WangLandauScalingFactor << "\n";
      textrestartFile << "Number-of-biasing-factors component " << i-SystemComponents.NumberOfFrameworks << ": " << SystemComponents.Lambda[i].binsize << "\n";
      textrestartFile << "Biasing-factors component " << i-SystemComponents.NumberOfFrameworks << ": "; 
      for(size_t j = 0; j < SystemComponents.Lambda[i].binsize; j++)
        textrestartFile << SystemComponents.Lambda[i].biasFactor[j] << " ";
      textrestartFile << "\n";
      textrestartFile << "Maximum-CF-Lambda-change component " << i-SystemComponents.NumberOfFrameworks << ": 0.50000\n"; //Zhao's note: continuous lambda not implemented
      textrestartFile << "Maximum-CBCF-Lambda-change component " << i-SystemComponents.NumberOfFrameworks << ": 0.50000\n"; //Zhao's note: continuous lambda not implemented
    }
    textrestartFile << "\n"; 
    textrestartFile << "Maximum-translation-change component " << i-SystemComponents.NumberOfFrameworks << ": " << SystemComponents.MaxTranslation[i].x << " " << SystemComponents.MaxTranslation[i].y << " " << SystemComponents.MaxTranslation[i].z << "\n";
    textrestartFile << "Maximum-translation-in-plane-change component " << i-SystemComponents.NumberOfFrameworks << ": " << "0.000000,0.000000,0.000000" << "\n";
    textrestartFile << "Maximum-rotation-change component " << i-SystemComponents.NumberOfFrameworks << ": " << SystemComponents.MaxRotation[i].x << " " << SystemComponents.MaxRotation[i].y << " " << SystemComponents.MaxRotation[i].z << "\n";
    textrestartFile << "\n";
  }
}

static inline void WriteCellInfo_Restart(Atoms* System, Components& SystemComponents, std::ofstream& textrestartFile, Boxsize& Box)
{
  //size_t Atomcount=0;
  textrestartFile << "Cell info:\n"; 
  textrestartFile << "========================================================================\n";
  textrestartFile << "number-of-unit-cells: 1 1 1\n"; //Zhao's note: not allowing for multiple unit cells for now//
  //First write Unit Box sizes//
  textrestartFile << "unit-cell-vector-a:" << " " << Box.Cell[0] << " " << Box.Cell[1] << " " << Box.Cell[2] << '\n';
  textrestartFile << "unit-cell-vector-b:" << " " << Box.Cell[3] << " " << Box.Cell[4] << " " << Box.Cell[5] << '\n';
  textrestartFile << "unit-cell-vector-c:" << " " << Box.Cell[6] << " " << Box.Cell[7] << " " << Box.Cell[8] << '\n';
  textrestartFile << "\n";
  //Then write Total Box sizes//
  textrestartFile << "cell-vector-a:" << " " << Box.Cell[0] << " " << Box.Cell[1] << " " << Box.Cell[2] << '\n';
  textrestartFile << "cell-vector-b:" << " " << Box.Cell[3] << " " << Box.Cell[4] << " " << Box.Cell[5] << '\n';
  textrestartFile << "cell-vector-c:" << " " << Box.Cell[6] << " " << Box.Cell[7] << " " << Box.Cell[8] << '\n';
  textrestartFile << "\n";
  //Cell lengths//
  textrestartFile << "cell-lengths:" << " " << Box.Cell[0] << " " << Box.Cell[4] << " " << Box.Cell[8] << " " << "\n";
  textrestartFile << "cell-angles:" << " " << 90.00 << " " << 90.00 << " " << 90.00 << " " << "\n";
  textrestartFile << "\n\n";
  //Maximum changes for MC-moves (Zhao's note: Not implemented)//
  textrestartFile << "Maximum changes for MC-moves:" << "\n"; 
  textrestartFile << "========================================================================" << "\n";
  textrestartFile << "Maximum-volume-change: 0.006250" << "\n";
  textrestartFile << "Maximum-Gibbs-volume-change: 0.025000" << "\n";
  textrestartFile << "Maximum-box-shape-change: 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000" << "\n";
  textrestartFile << "\n\n";
  //Acceptance targets for MC-moves (Zhao's note: Not implemented)//
  textrestartFile << "Acceptance targets for MC-moves:" << "\n";
  textrestartFile << "========================================================================" << "\n";
  textrestartFile << "Target-volume-change: 0.500000" << "\n"; 
  textrestartFile << "Target-box-shape-change: 0.500000" << "\n";
  textrestartFile << "Target-Gibbs-volume-change: 0.500000" << "\n";
  textrestartFile << "\n\n";
  //Write Component Data//
  WriteComponent_Restart(System, SystemComponents, textrestartFile, Box);
}

static inline void create_movie_file(Atoms* System, Components& SystemComponents, Boxsize& HostBox, std::vector<std::string> AtomNames, size_t SystemIndex)
{
  //std::ofstream textrestartFile{};
  std::string dirname="Movies/System_" + std::to_string(SystemIndex) + "/";
  std::string fname  = dirname + "/" + "result_" + std::to_string(SystemComponents.CURRENTCYCLE) + ".data";
  std::filesystem::path cwd = std::filesystem::current_path();

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);

  std::ofstream textrestartFile = std::ofstream(fileName, std::ios::out);

  WriteBox_LAMMPS(System, SystemComponents, SystemComponents.FF, HostBox, textrestartFile, AtomNames);
  WriteAtoms_LAMMPS(System, SystemComponents, HostBox, textrestartFile, AtomNames);
  textrestartFile.close();
}

//============================================================================//
//  3D spatial adsorbate-density grid (RASPA2 ComputeDensityProfile3DVTKGrid)  //
//                                                                            //
//  Accumulate_DensityGrid: bin every adsorbate atom of every adsorbate       //
//  component into an nx*ny*nz histogram over the unit cell (one snapshot).    //
//  WriteDensityGrid: normalize to average number density and write one        //
//  Gaussian .cube (triclinic-correct) + one VTK StructuredPoints (orthorhombic//
//  assumption) file per adsorbate component, under DensityGrid/System_<id>/.  //
//============================================================================//

//Accumulate one snapshot of adsorbate positions into the per-component grid.//
//Requires SystemComponents.HostSystem to already hold the current device positions//
//(copied via Copy_AtomData_from_Device) and HostBox.InverseCell to be host-resident.//
static inline void Accumulate_DensityGrid(Components& SystemComponents, Boxsize& HostBox)
{
  const int nx = SystemComponents.DensityGridPoints.x;
  const int ny = SystemComponents.DensityGridPoints.y;
  const int nz = SystemComponents.DensityGridPoints.z;
  const size_t Nvoxel = (size_t)nx * (size_t)ny * (size_t)nz;

  //Lazily allocate per-component accumulators (adsorbate components only)//
  if(SystemComponents.DensityGrid.size() != (size_t)SystemComponents.NComponents.x)
  {
    SystemComponents.DensityGrid.assign((size_t)SystemComponents.NComponents.x, std::vector<double>());
    for(int comp = SystemComponents.NComponents.y; comp < SystemComponents.NComponents.x; comp++)
      SystemComponents.DensityGrid[comp].assign(Nvoxel, 0.0);
  }

  const double* Inv  = HostBox.InverseCell; //fractional s = InverseCell * pos (see ewald_preparation.h:52)//
  const double* Cell = HostBox.Cell;        //cartesian = Cell * fractional (column-major: Cell[0/4/8] = Lx/Ly/Lz)//
  const bool useCenter = SystemComponents.DensityGridMolecularCenter;
  for(int comp = SystemComponents.NComponents.y; comp < SystemComponents.NComponents.x; comp++)
  {
    Atoms& A = SystemComponents.HostSystem[comp];
    std::vector<double>& Grid = SystemComponents.DensityGrid[comp];
    if(Grid.size() != Nvoxel) continue;
    if(useCenter)
    {
      //Bin ONE point per molecule: the molecular geometric center, computed PBC-safely via
      //minimum image relative to the molecule's first atom. For a symmetric molecule (CO2:
      //O=C=O) this is exactly the central atom (C) = center of mass, for any orientation;
      //for a monatomic species (Xe) it is the atom itself. The field then integrates to the
      //average MOLECULE count (= reported loading), not the atom count.//
      const size_t ms = (A.Molsize > 0) ? A.Molsize : 1;
      for(size_t m0 = 0; m0 + ms <= A.size; m0 += ms)
      {
        const double3 ref = A.pos[m0];
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for(size_t k = 0; k < ms; k++)
        {
          const double dx = A.pos[m0+k].x - ref.x, dy = A.pos[m0+k].y - ref.y, dz = A.pos[m0+k].z - ref.z;
          double fx = Inv[0]*dx + Inv[3]*dy + Inv[6]*dz;
          double fy = Inv[1]*dx + Inv[4]*dy + Inv[7]*dz;
          double fz = Inv[2]*dx + Inv[5]*dy + Inv[8]*dz;
          fx -= std::round(fx); fy -= std::round(fy); fz -= std::round(fz);   //minimum image//
          cx += Cell[0]*fx + Cell[3]*fy + Cell[6]*fz;
          cy += Cell[1]*fx + Cell[4]*fy + Cell[7]*fz;
          cz += Cell[2]*fx + Cell[5]*fy + Cell[8]*fz;
        }
        const double px = ref.x + cx/(double)ms, py = ref.y + cy/(double)ms, pz = ref.z + cz/(double)ms;
        double sx = Inv[0]*px + Inv[3]*py + Inv[6]*pz;
        double sy = Inv[1]*px + Inv[4]*py + Inv[7]*pz;
        double sz = Inv[2]*px + Inv[5]*py + Inv[8]*pz;
        sx -= std::floor(sx); sy -= std::floor(sy); sz -= std::floor(sz);
        int ix = (int)(sx * nx); if(ix >= nx) ix = nx - 1; if(ix < 0) ix = 0;
        int iy = (int)(sy * ny); if(iy >= ny) iy = ny - 1; if(iy < 0) iy = 0;
        int iz = (int)(sz * nz); if(iz >= nz) iz = nz - 1; if(iz < 0) iz = 0;
        Grid[(size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz] += 1.0;
      }
    }
    else
    {
      for(size_t j = 0; j < A.size; j++)
      {
        const double3 p = A.pos[j];
        double sx = Inv[0]*p.x + Inv[3]*p.y + Inv[6]*p.z;
        double sy = Inv[1]*p.x + Inv[4]*p.y + Inv[7]*p.z;
        double sz = Inv[2]*p.x + Inv[5]*p.y + Inv[8]*p.z;
        //wrap fractional coordinate into [0,1)//
        sx -= std::floor(sx); sy -= std::floor(sy); sz -= std::floor(sz);
        int ix = (int)(sx * nx); if(ix >= nx) ix = nx - 1; if(ix < 0) ix = 0;
        int iy = (int)(sy * ny); if(iy >= ny) iy = ny - 1; if(iy < 0) iy = 0;
        int iz = (int)(sz * nz); if(iz >= nz) iz = nz - 1; if(iz < 0) iz = 0;
        Grid[(size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz] += 1.0;
      }
    }
  }
  SystemComponents.DensityGridSamples++;
}

//Write the accumulated grid for every adsorbate component (number density, molecules/Angstrom^3).//
static inline void WriteDensityGrid(Components& SystemComponents, Boxsize& HostBox, size_t SystemIndex)
{
  if(!SystemComponents.ComputeDensityGrid) return;
  if(SystemComponents.DensityGridSamples == 0 || SystemComponents.DensityGrid.empty()) return;

  const int nx = SystemComponents.DensityGridPoints.x;
  const int ny = SystemComponents.DensityGridPoints.y;
  const int nz = SystemComponents.DensityGridPoints.z;
  const size_t Nvoxel  = (size_t)nx * (size_t)ny * (size_t)nz;
  const double Nsamp   = (double)SystemComponents.DensityGridSamples;
  const double voxelVolume = HostBox.Volume / (double)Nvoxel; //Angstrom^3 per voxel//
  const double norm    = 1.0 / (Nsamp * voxelVolume);         //counts -> molecules/Angstrom^3//
  const double ANG_TO_BOHR = 1.8897259885789;                 //Gaussian .cube geometry is in Bohr//

  //Lattice vectors as rows of HostBox.Cell (see matrix_multiply_by_vector in maths.cuh)//
  const double* C = HostBox.Cell;
  const double ax = C[0], ay = C[1], az = C[2]; //a-vector//
  const double bx = C[3], by = C[4], bz = C[5]; //b-vector//
  const double cx = C[6], cy = C[7], cz = C[8]; //c-vector//
  const double alen = std::sqrt(ax*ax + ay*ay + az*az);
  const double blen = std::sqrt(bx*bx + by*by + bz*bz);
  const double clen = std::sqrt(cx*cx + cy*cy + cz*cz);

  std::string dirname = "DensityGrid/System_" + std::to_string(SystemIndex) + "/";
  std::filesystem::path cwd = std::filesystem::current_path();
  std::filesystem::create_directories(cwd / dirname);

  for(int comp = SystemComponents.NComponents.y; comp < SystemComponents.NComponents.x; comp++)
  {
    std::vector<double>& Grid = SystemComponents.DensityGrid[comp];
    if(Grid.size() != Nvoxel) continue;
    std::string mname = (comp < (int)SystemComponents.MoleculeName.size()) ? SystemComponents.MoleculeName[comp] : ("Component_" + std::to_string(comp));
    std::string base  = dirname + "DensityGrid_" + std::to_string(comp) + "_" + mname;

    //---- Gaussian .cube : triclinic-correct via voxel vectors (geometry in Bohr) ----//
    {
      std::ofstream f(cwd / (base + ".cube"), std::ios::out);
      f << "gRASPA 3D adsorbate density grid (RASPA2 ComputeDensityProfile3DVTKGrid parity)\n";
      f << "Component " << comp << " (" << mname << "); scalar = number density [molecules/Angstrom^3]; "
        << SystemComponents.DensityGridSamples << " snapshots\n";
      //Natoms (1 marker atom at origin), origin in Bohr//
      f << std::fixed << std::setprecision(6);
      f << "    1    0.000000    0.000000    0.000000\n";
      f << std::scientific << std::setprecision(6);
      f << nx << " " << (ax/nx)*ANG_TO_BOHR << " " << (ay/nx)*ANG_TO_BOHR << " " << (az/nx)*ANG_TO_BOHR << "\n";
      f << ny << " " << (bx/ny)*ANG_TO_BOHR << " " << (by/ny)*ANG_TO_BOHR << " " << (bz/ny)*ANG_TO_BOHR << "\n";
      f << nz << " " << (cx/nz)*ANG_TO_BOHR << " " << (cy/nz)*ANG_TO_BOHR << " " << (cz/nz)*ANG_TO_BOHR << "\n";
      //marker atom: hydrogen at cell origin//
      f << "    1    1.000000    0.000000    0.000000    0.000000\n";
      //volumetric data: ix outer, iy middle, iz inner (matches our flat index)//
      for(int ix = 0; ix < nx; ix++)
        for(int iy = 0; iy < ny; iy++)
        {
          int col = 0;
          for(int iz = 0; iz < nz; iz++)
          {
            double val = Grid[(size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz] * norm;
            f << " " << std::setw(13) << val;
            if(++col % 6 == 0) { f << "\n"; }
          }
          if(col % 6 != 0) f << "\n";
        }
      f.close();
    }

    //---- VTK StructuredPoints (ASCII): assumes ORTHORHOMBIC cell, geometry in Angstrom ----//
    {
      std::ofstream f(cwd / (base + ".vtk"), std::ios::out);
      f << "# vtk DataFile Version 3.0\n";
      f << "gRASPA 3D adsorbate density grid, component " << comp << " (" << mname
        << "), number density [molecules/Angstrom^3], " << SystemComponents.DensityGridSamples << " snapshots\n";
      f << "ASCII\n";
      f << "DATASET STRUCTURED_POINTS\n";
      f << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
      f << "ORIGIN 0 0 0\n";
      f << std::scientific << std::setprecision(8);
      f << "SPACING " << (alen/nx) << " " << (blen/ny) << " " << (clen/nz) << "\n";
      f << "POINT_DATA " << Nvoxel << "\n";
      f << "SCALARS number_density double 1\n";
      f << "LOOKUP_TABLE default\n";
      //VTK ordering: X fastest, then Y, then Z//
      for(int iz = 0; iz < nz; iz++)
        for(int iy = 0; iy < ny; iy++)
          for(int ix = 0; ix < nx; ix++)
            f << (Grid[(size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz] * norm) << "\n";
      f.close();
    }
    fprintf(SystemComponents.OUTPUT, "Wrote 3D density grid for component %d (%s): %s.cube / .vtk (%d x %d x %d, %zu snapshots)\n",
            comp, mname.c_str(), base.c_str(), nx, ny, nz, SystemComponents.DensityGridSamples);
  }
}

static inline void copyFile(const std::string& sourcePath, const std::string& destinationPath)
{
  std::ifstream source(sourcePath, std::ios::binary);
  std::ofstream destination(destinationPath, std::ios::binary);

  if (!source.is_open() || !destination.is_open()) 
  {
    std::cerr << "Cannot Open " << sourcePath << " or " << destinationPath << std::endl;
    return; // false;
  }

  destination << source.rdbuf();

  source.close();
  destination.close();
  //return true;
}

static inline void create_Restart_file(size_t Cycle, Atoms* System, Components& SystemComponents, ForceField& FF, Boxsize& Box, std::vector<std::string>& AtomNames, size_t SystemIndex)
{
  std::ofstream textrestartFile{};
  std::string dirname="Restart/System_" + std::to_string(SystemIndex) + "/";
  std::string fname  = dirname + "/" + "restartfile";
  std::string backup_fname = dirname + "/" + "previous_restartfile";
  std::filesystem::path cwd = std::filesystem::current_path();
  
  //Keep a previous copy of the current restart file//
  copyFile(fname, backup_fname);

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);
  textrestartFile = std::ofstream(fileName, std::ios::out);
  WriteCellInfo_Restart(System, SystemComponents, textrestartFile, Box);
  WriteAtoms_Restart(System, SystemComponents, textrestartFile, AtomNames);
  textrestartFile.close();
}

static inline void WriteAllData(Atoms* System, Components& SystemComponents, std::ofstream& textrestartFile, std::vector<std::string>& AtomNames, size_t i)
{
  //size_t Atomcount=0;
  Atoms Data = System[i];
  size_t molsize = SystemComponents.Moleculesize[i];
  //First write positions//
  //fprintf(SystemComponents.OUTPUT, "Writing AllData for Component %zu, There are %zu atoms, molsize: %zu\n", i, Data.size, molsize);
  textrestartFile << "x y z charge scale scaleCoul Type" << '\n';
  for(size_t j = 0; j < Data.size; j++)
    textrestartFile << " " << Data.MolID[j] << " " << j - Data.MolID[j]*molsize << " " << Data.pos[j].x << "  " << Data.pos[j].y << "  " << Data.pos[j].z << " " << Data.charge[j] << " " << Data.scale[j] << " " << Data.scaleCoul[j] << " " << Data.Type[j] << '\n';
}

static inline void Write_All_Adsorbate_data(size_t Cycle, Atoms* System, Components& SystemComponents, ForceField& FF, Boxsize& Box, std::vector<std::string>& AtomNames, size_t SystemIndex)
{
  //std::ofstream textrestartFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    std::string dirname="AllData/System_" + std::to_string(SystemIndex) + "/";
    std::string fname  = dirname + "/" + "Component_" + std::to_string(i) + ".data";

    std::filesystem::path directoryName = cwd /dirname;
    std::filesystem::path fileName = cwd /fname;
    std::filesystem::create_directories(directoryName);
    std::ofstream textrestartFile = std::ofstream(fileName, std::ios::out);
    WriteAllData(System, SystemComponents, textrestartFile, AtomNames, i);
    textrestartFile.close();
  }
}

//Zhao's note on 112623, moved the following to fxn_main.h//
/*
static inline void Write_Lambda(size_t Cycle, Components& SystemComponents, size_t SystemIndex)
{
  std::ofstream textrestartFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname="Lambda/System_" + std::to_string(SystemIndex) + "/";
  std::string fname  = dirname + "/" + "Lambda_Histogram.data";

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);
  textrestartFile = std::ofstream(fileName, std::ios::out);
  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    if(SystemComponents.hasfractionalMolecule[i])
    {
      textrestartFile << "Component " << i << ": " << SystemComponents.MoleculeName[i] << '\n';
      textrestartFile << "BIN SIZE : " << SystemComponents.Lambda[i].binsize << '\n';
      textrestartFile << "BIN WIDTH: " << SystemComponents.Lambda[i].delta << '\n';
      textrestartFile << "WL SCALING FACTOR: " << SystemComponents.Lambda[i].WangLandauScalingFactor << '\n';
      textrestartFile << "FRACTIONAL MOLECULE ID: " << SystemComponents.Lambda[i].FractionalMoleculeID << '\n';
      textrestartFile << "CURRENT BIN: " << SystemComponents.Lambda[i].currentBin << '\n';
      textrestartFile << "BINS: ";
      for(size_t j = 0; j < SystemComponents.Lambda[i].binsize; j++)
        textrestartFile << j << " ";
      textrestartFile << "\nHistogram: ";
      for(size_t j = 0; j < SystemComponents.Lambda[i].binsize; j++)
        textrestartFile << SystemComponents.Lambda[i].Histogram[j] << " ";
      textrestartFile << "\nBIAS FACTOR: ";
      for(size_t j = 0; j < SystemComponents.Lambda[i].binsize; j++)
        textrestartFile << SystemComponents.Lambda[i].biasFactor[j] << " ";
    }
  }
  textrestartFile.close();
}

static inline void Write_TMMC(size_t Cycle, Components& SystemComponents, size_t SystemIndex)
{
  std::ofstream textTMMCFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname="TMMC/System_" + std::to_string(SystemIndex) + "/";
  std::string fname  = dirname + "/" + "TMMC_Statistics.data";

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);
  textTMMCFile = std::ofstream(fileName, std::ios::out);
  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    if(SystemComponents.Tmmc[i].DoTMMC)
    {
      textTMMCFile << "Component " << i << ": " << SystemComponents.MoleculeName[i] << " -> Updated " << SystemComponents.Tmmc[i].TMUpdateTimes << " times \n";
      textTMMCFile << "Min Macrostate : " << SystemComponents.Tmmc[i].MinMacrostate << '\n';
      textTMMCFile << "Max Macrostate : " << SystemComponents.Tmmc[i].MaxMacrostate << '\n';
      textTMMCFile << "Wang-Landau Factor : " << SystemComponents.Tmmc[i].WLFactor << '\n';
      textTMMCFile << "N NMol Bin CM[-1] CM[0] CM[1] WLBias ln_g TMBias lnpi Forward_lnpi Reverse_lnpi Histogram" << '\n';
      for(size_t j = 0; j < SystemComponents.Tmmc[i].Histogram.size(); j++)
      {
        size_t N   = j / SystemComponents.Tmmc[i].nbinPerMacrostate;
        size_t bin = j % SystemComponents.Tmmc[i].nbinPerMacrostate;
        textTMMCFile << j << " " << N << " " << bin << " " << SystemComponents.Tmmc[i].CMatrix[j].x << " " << SystemComponents.Tmmc[i].CMatrix[j].y << " " << SystemComponents.Tmmc[i].CMatrix[j].z << " " << SystemComponents.Tmmc[i].WLBias[j] << " " << SystemComponents.Tmmc[i].ln_g[j] << " " << SystemComponents.Tmmc[i].TMBias[j] << " " << SystemComponents.Tmmc[i].lnpi[j] << " " << SystemComponents.Tmmc[i].forward_lnpi[j] << " " << SystemComponents.Tmmc[i].reverse_lnpi[j] << " " << SystemComponents.Tmmc[i].Histogram[j] << '\n';
      }
    }
  }
  textrestartFile.close();
}
*/

//============================================================================//
//  Widom insertion-energy grid                                               //
//                                                                            //
//  Spatial map of Widom-sample insertion energies (VDW+Real+Ewald+Tail) and  //
//  Rosenbluth weights, binned by the selected trial's first-bead position    //
//  (fractional coords over the full simulation cell, like the density grid). //
//  3D generalization of RASPA2's 1D ComputeFreeEnergyProfile, whose stored   //
//  quantity is betaF(q) = -ln<W(q)>. Per voxel we keep count, sum(W),        //
//  sum(E_total*W), sum(E_ewald*W) and min(E_total), so the written maps are: //
//    betaF          = -ln( sum(W)/count )                    [dimensionless] //
//    <E>_W          = sum(E_total*W)/sum(W)                  [Kelvin]        //
//    <E_ewald>_W    = sum(E_ewald*W)/sum(W)                  [Kelvin]        //
//  <E>_W is the local version of gRASPA's printed "AVERAGE WIDOM ADSORPTION  //
//  ENERGY" (which is sum(E*W)/sum(W) over ALL samples) -> global cross-check.//
//  Failed CBMC constructions count toward the global sample total with W=0   //
//  (matching RecordRosen's denominator) but cannot be binned: a failed       //
//  multi-trial CBMC sample has no single defined position.                   //
//============================================================================//

//Record ONE production Widom sample. For successful constructions, pos = selected//
//trial's first bead (Cartesian, device-fetched by the caller) and InvCell = the  //
//current device InverseCell (9 doubles, row-major); both ignored on failure.     //
static inline void Accumulate_WidomEnergyGrid(Components& SystemComponents, bool Success, const double* InvCell,
                                              double3 pos, MoveEnergy& E, double Rosenbluth, size_t comp, size_t SystemIndex)
{
  const int nx = SystemComponents.WidomEnergyGridPoints.x;
  const int ny = SystemComponents.WidomEnergyGridPoints.y;
  const int nz = SystemComponents.WidomEnergyGridPoints.z;
  const size_t Nvoxel = (size_t)nx * (size_t)ny * (size_t)nz;
  const size_t NComp  = (size_t)SystemComponents.NComponents.x;

  //Lazily allocate accumulators (adsorbate components only get voxel arrays)//
  if(SystemComponents.WidomGridCounts.size() != NComp)
  {
    SystemComponents.WidomGridSumW.assign(NComp, std::vector<double>());
    SystemComponents.WidomGridSumEW.assign(NComp, std::vector<double>());
    SystemComponents.WidomGridSumEwaldW.assign(NComp, std::vector<double>());
    SystemComponents.WidomGridMinE.assign(NComp, std::vector<double>());
    SystemComponents.WidomGridCounts.assign(NComp, std::vector<size_t>());
    SystemComponents.WidomGridTotalSamples.assign(NComp, 0);
    SystemComponents.WidomGridFailedSamples.assign(NComp, 0);
    SystemComponents.WidomGridGlobalSumW.assign(NComp, 0.0);
    SystemComponents.WidomGridGlobalSumEW.assign(NComp, 0.0);
    SystemComponents.WidomGridRawFiles.assign(NComp, nullptr);
    for(int c = SystemComponents.NComponents.y; c < SystemComponents.NComponents.x; c++)
    {
      SystemComponents.WidomGridSumW[c].assign(Nvoxel, 0.0);
      SystemComponents.WidomGridSumEW[c].assign(Nvoxel, 0.0);
      SystemComponents.WidomGridSumEwaldW[c].assign(Nvoxel, 0.0);
      SystemComponents.WidomGridMinE[c].assign(Nvoxel, 1.0e30);
      SystemComponents.WidomGridCounts[c].assign(Nvoxel, 0);
    }
  }
  if(comp >= NComp || SystemComponents.WidomGridSumW[comp].size() != Nvoxel) return;

  //Global sums: every attempt counts (failed ones as W=0), mirroring RecordRosen//
  SystemComponents.WidomGridTotalSamples[comp]++;
  if(!Success) { SystemComponents.WidomGridFailedSamples[comp]++; return; }

  const double Etot   = E.total();
  const double Eewald = E.HGEwaldE + E.GGEwaldE;
  SystemComponents.WidomGridGlobalSumW[comp]  += Rosenbluth;
  SystemComponents.WidomGridGlobalSumEW[comp] += Etot * Rosenbluth;

  //Bin by fractional coordinate (same convention as Accumulate_DensityGrid)//
  double sx = InvCell[0]*pos.x + InvCell[3]*pos.y + InvCell[6]*pos.z;
  double sy = InvCell[1]*pos.x + InvCell[4]*pos.y + InvCell[7]*pos.z;
  double sz = InvCell[2]*pos.x + InvCell[5]*pos.y + InvCell[8]*pos.z;
  sx -= std::floor(sx); sy -= std::floor(sy); sz -= std::floor(sz);
  int ix = (int)(sx * nx); if(ix >= nx) ix = nx - 1; if(ix < 0) ix = 0;
  int iy = (int)(sy * ny); if(iy >= ny) iy = ny - 1; if(iy < 0) iy = 0;
  int iz = (int)(sz * nz); if(iz >= nz) iz = nz - 1; if(iz < 0) iz = 0;
  const size_t voxel = (size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz;

  SystemComponents.WidomGridSumW[comp][voxel]      += Rosenbluth;
  SystemComponents.WidomGridSumEW[comp][voxel]     += Etot * Rosenbluth;
  SystemComponents.WidomGridSumEwaldW[comp][voxel] += Eewald * Rosenbluth;
  SystemComponents.WidomGridCounts[comp][voxel]++;
  if(Etot < SystemComponents.WidomGridMinE[comp][voxel]) SystemComponents.WidomGridMinE[comp][voxel] = Etot;

  //Optional raw per-sample dump: "the Widom energy at each sampled point"//
  if(SystemComponents.WidomEnergyGridRawDump)
  {
    if(SystemComponents.WidomGridRawFiles[comp] == nullptr)
    {
      std::string dirname = "WidomEnergyGrid/System_" + std::to_string(SystemIndex) + "/";
      std::filesystem::path cwd = std::filesystem::current_path();
      std::filesystem::create_directories(cwd / dirname);
      std::string mname = (comp < SystemComponents.MoleculeName.size()) ? SystemComponents.MoleculeName[comp] : ("Component_" + std::to_string(comp));
      std::string fname = dirname + "WidomSamples_" + std::to_string(comp) + "_" + mname + ".dat";
      SystemComponents.WidomGridRawFiles[comp] = fopen((cwd / fname).c_str(), "w");
      if(SystemComponents.WidomGridRawFiles[comp] != nullptr)
        fprintf(SystemComponents.WidomGridRawFiles[comp],
                "# x[A] y[A] z[A] RosenbluthW E_total[K] E_vdw[K] E_real[K] E_ewald[K] E_tail[K] (first-bead position; energies incl. Ewald reciprocal)\n");
      else
        fprintf(SystemComponents.OUTPUT, "WARNING: WidomEnergyGrid: could not open raw dump file %s, raw dump disabled for this component\n", fname.c_str());
    }
    const double K = 1.2027242847; //internal energy units -> Kelvin//
    if(SystemComponents.WidomGridRawFiles[comp] != nullptr)
      fprintf(SystemComponents.WidomGridRawFiles[comp], "%g %g %g %g %g %g %g %g %g\n",
              pos.x, pos.y, pos.z, Rosenbluth, Etot * K,
              (E.HGVDW + E.GGVDW) * K, (E.HGReal + E.GGReal) * K, Eewald * K, E.TailE * K);
  }
}

//Write the accumulated Widom energy grid: one ASCII .grid (all per-voxel stats) +      //
//three .cube maps per adsorbate component, plus a global cross-check line to OUTPUT.   //
static inline void WriteWidomEnergyGrid(Components& SystemComponents, Boxsize& HostBox, size_t SystemIndex)
{
  if(!SystemComponents.ComputeWidomEnergyGrid) return;
  if(SystemComponents.WidomGridCounts.empty()) return;

  const int nx = SystemComponents.WidomEnergyGridPoints.x;
  const int ny = SystemComponents.WidomEnergyGridPoints.y;
  const int nz = SystemComponents.WidomEnergyGridPoints.z;
  const size_t Nvoxel = (size_t)nx * (size_t)ny * (size_t)nz;
  const double K = 1.2027242847;             //internal energy units -> Kelvin//
  const double ANG_TO_BOHR = 1.8897259885789;
  const double EMPTY_E_K   = 2.0e7 * K;      //unsampled-voxel marker in energy cubes (EnergyGrid cap convention)//
  const double EMPTY_BETAF = 50.0;           //unsampled-voxel marker in betaF cube (W ~ e^-50: inaccessible)//

  //Lattice vectors as rows of HostBox.Cell//
  const double* C = HostBox.Cell;

  std::string dirname = "WidomEnergyGrid/System_" + std::to_string(SystemIndex) + "/";
  std::filesystem::path cwd = std::filesystem::current_path();
  std::filesystem::create_directories(cwd / dirname);

  for(int comp = SystemComponents.NComponents.y; comp < SystemComponents.NComponents.x; comp++)
  {
    if(SystemComponents.WidomGridCounts[comp].size() != Nvoxel) continue;
    if(SystemComponents.WidomGridTotalSamples[comp] == 0) continue;
    std::string mname = (comp < (int)SystemComponents.MoleculeName.size()) ? SystemComponents.MoleculeName[comp] : ("Component_" + std::to_string(comp));
    std::string base  = dirname + "WidomEnergyGrid_" + std::to_string(comp) + "_" + mname;

    //---- ASCII per-voxel stats: voxel centers in Cartesian ----//
    FILE* fgrid = fopen((cwd / (base + ".grid")).c_str(), "w");
    if(fgrid == nullptr)
    { fprintf(SystemComponents.OUTPUT, "WARNING: WidomEnergyGrid: could not open %s.grid for writing, skipping\n", base.c_str()); }
    else
    {
      fprintf(fgrid, "# Widom insertion-energy grid, component %d (%s), %zu samples (%zu failed, unbinned)\n",
              comp, mname.c_str(), SystemComponents.WidomGridTotalSamples[comp], SystemComponents.WidomGridFailedSamples[comp]);
      fprintf(fgrid, "# x[A] y[A] z[A] count meanW betaF(-ln meanW) E_W[K] E_ewald_W[K] minE[K]   (E_W = sum(E*W)/sum(W); empty voxels: count=0, other columns '?')\n");
      for(int ix = 0; ix < nx; ix++)
        for(int iy = 0; iy < ny; iy++)
          for(int iz = 0; iz < nz; iz++)
          {
            const size_t v = (size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz;
            const double fx = (ix + 0.5) / nx, fy = (iy + 0.5) / ny, fz = (iz + 0.5) / nz;
            const double px = C[0]*fx + C[3]*fy + C[6]*fz;
            const double py = C[1]*fx + C[4]*fy + C[7]*fz;
            const double pz = C[2]*fx + C[5]*fy + C[8]*fz;
            const size_t n = SystemComponents.WidomGridCounts[comp][v];
            if(n == 0)
            { fprintf(fgrid, "%g %g %g 0 ? ? ? ? ?\n", px, py, pz); continue; }
            const double meanW = SystemComponents.WidomGridSumW[comp][v] / (double)n;
            const double sumW  = SystemComponents.WidomGridSumW[comp][v];
            const double E_W   = (sumW > 0.0) ? SystemComponents.WidomGridSumEW[comp][v]      / sumW * K : 0.0;
            const double Eew_W = (sumW > 0.0) ? SystemComponents.WidomGridSumEwaldW[comp][v]  / sumW * K : 0.0;
            const double betaF = (meanW > 0.0) ? -std::log(meanW) : EMPTY_BETAF;
            fprintf(fgrid, "%g %g %g %zu %g %g %g %g %g\n", px, py, pz, n, meanW, betaF, E_W, Eew_W,
                    SystemComponents.WidomGridMinE[comp][v] * K);
          }
      fclose(fgrid);
    }

    //---- three .cube maps (triclinic-correct voxel vectors, geometry in Bohr) ----//
    struct CubeSpec { const char* tag; const char* what; };
    const CubeSpec specs[3] = {
      {"EnergyW", "Rosenbluth-weighted mean Widom insertion energy sum(E*W)/sum(W) [Kelvin]"},
      {"BetaF",   "betaF = -ln( sum(W)/count ) per voxel [dimensionless] (RASPA2 FreeEnergyProfile quantity, 3D)"},
      {"EwaldW",  "Rosenbluth-weighted mean Ewald (reciprocal+exclusion) part of the insertion energy [Kelvin]"},
    };
    for(int s = 0; s < 3; s++)
    {
      std::ofstream f(cwd / (base + "_" + specs[s].tag + ".cube"), std::ios::out);
      f << "gRASPA Widom insertion-energy grid, component " << comp << " (" << mname << ")\n";
      f << specs[s].what << "; " << SystemComponents.WidomGridTotalSamples[comp] << " samples; unsampled voxels = "
        << (s == 1 ? EMPTY_BETAF : EMPTY_E_K) << "\n";
      f << std::fixed << std::setprecision(6);
      f << "    1    0.000000    0.000000    0.000000\n";
      f << std::scientific << std::setprecision(6);
      f << nx << " " << (C[0]/nx)*ANG_TO_BOHR << " " << (C[1]/nx)*ANG_TO_BOHR << " " << (C[2]/nx)*ANG_TO_BOHR << "\n";
      f << ny << " " << (C[3]/ny)*ANG_TO_BOHR << " " << (C[4]/ny)*ANG_TO_BOHR << " " << (C[5]/ny)*ANG_TO_BOHR << "\n";
      f << nz << " " << (C[6]/nz)*ANG_TO_BOHR << " " << (C[7]/nz)*ANG_TO_BOHR << " " << (C[8]/nz)*ANG_TO_BOHR << "\n";
      f << "    1    1.000000    0.000000    0.000000    0.000000\n";
      for(int ix = 0; ix < nx; ix++)
        for(int iy = 0; iy < ny; iy++)
        {
          int col = 0;
          for(int iz = 0; iz < nz; iz++)
          {
            const size_t v = (size_t)ix*(size_t)ny*(size_t)nz + (size_t)iy*(size_t)nz + (size_t)iz;
            const size_t n = SystemComponents.WidomGridCounts[comp][v];
            double val;
            if(n == 0) val = (s == 1) ? EMPTY_BETAF : EMPTY_E_K;
            else
            {
              const double sumW  = SystemComponents.WidomGridSumW[comp][v];
              const double meanW = sumW / (double)n;
              if(s == 0)      val = (sumW > 0.0) ? SystemComponents.WidomGridSumEW[comp][v]     / sumW * K : EMPTY_E_K;
              else if(s == 1) val = (meanW > 0.0) ? -std::log(meanW) : EMPTY_BETAF;
              else            val = (sumW > 0.0) ? SystemComponents.WidomGridSumEwaldW[comp][v] / sumW * K : EMPTY_E_K;
            }
            f << " " << std::setw(13) << val;
            if(++col % 6 == 0) { f << "\n"; }
          }
          if(col % 6 != 0) f << "\n";
        }
      f.close();
    }

    //---- close raw dump + global cross-check line (must match printed Widom statistics) ----//
    if(comp < (int)SystemComponents.WidomGridRawFiles.size() && SystemComponents.WidomGridRawFiles[comp] != nullptr)
    { fclose(SystemComponents.WidomGridRawFiles[comp]); SystemComponents.WidomGridRawFiles[comp] = nullptr; }
    const double globW  = SystemComponents.WidomGridGlobalSumW[comp];
    const size_t nTot   = SystemComponents.WidomGridTotalSamples[comp];
    fprintf(SystemComponents.OUTPUT,
            "Widom energy grid, component %d (%s): %zu samples (%zu failed), global <W> = %.10g, global sum(E*W)/sum(W) = %.5f [K] (cross-check vs printed Widom statistics); wrote %s.grid + _EnergyW/_BetaF/_EwaldW.cube\n",
            comp, mname.c_str(), nTot, SystemComponents.WidomGridFailedSamples[comp],
            (nTot > 0) ? globW / (double)nTot : 0.0,
            (globW > 0.0) ? SystemComponents.WidomGridGlobalSumEW[comp] / globW * K : 0.0,
            base.c_str());
  }
}
