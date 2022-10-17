#include "axpy.h"
#include "print_statistics.h"
#include "mc_translation.h"
#include "mc_insertion_deletion.h"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <filesystem>

double cuSoA(int Cycles, Components& SystemComponents, Boxsize Box, Atoms* d_a, Atoms Mol, Atoms NewMol, ForceField FF, double* y, double* dUdlambda, RandomNumber Random, WidomStruct Widom, Units Constants, double init_energy, bool DualPrecision)
{
  
  double tot = 0.0;

  double running_energy = 0.0;

  size_t WidomCount = 0;

  bool DEBUG = false;
  size_t transCount=0;
  printf("There are %zu Molecules, %zu Frameworks\n",SystemComponents.TotalNumberOfMolecules, SystemComponents.NumberOfFrameworks);

  for(size_t i = 0; i < Cycles; i++)
  {
    //Randomly Select an Adsorbate Molecule and determine its Component: MoleculeID --> Component
    //if((SystemComponents.TotalNumberOfMolecules - SystemComponents.NumberOfFrameworks) == 0)
    //  continue;
    size_t SelectedMolecule = (size_t) (get_random_from_zero_to_one()*(SystemComponents.TotalNumberOfMolecules-SystemComponents.NumberOfFrameworks));
    size_t comp = SystemComponents.NumberOfFrameworks; // When selecting components, skip the component 0 (because it is the framework)
    size_t SelectedMolInComponent = SelectedMolecule; size_t totalsize= 0;
    for(size_t ijk = SystemComponents.NumberOfFrameworks; ijk < SystemComponents.Total_Components; ijk++) //Assuming Framework atoms are the top in the Atoms array
    {
      if(SelectedMolInComponent == 0) break;
      totalsize += SystemComponents.NumberOfMolecule_for_Component[ijk];
      if(SelectedMolInComponent >= totalsize)
      {
        comp++;
        SelectedMolInComponent -= SystemComponents.NumberOfMolecule_for_Component[ijk];
      }
    }

    if(SystemComponents.NumberOfMolecule_for_Component[comp] == 0){ //no molecule in the system for this species
      running_energy += Insertion(Box, SystemComponents, d_a, Mol, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, DualPrecision);
      continue;
    }

    double RANDOMNUMBER = get_random_from_zero_to_one();
    if(RANDOMNUMBER < SystemComponents.Moves[comp].TranslationProb)
    {
      transCount++;
      //PERFORM TRANSLATION MOVE//
      running_energy += Translation_Move(Box, SystemComponents, d_a, Mol, NewMol, FF, y, dUdlambda, Random, SelectedMolInComponent, comp);
      if(DEBUG){printf("After Translation: running energy: %.10f\n", running_energy);}
    }
    else if(RANDOMNUMBER < SystemComponents.Moves[comp].TranslationProb + SystemComponents.Moves[comp].RotationProb) //Rotation
    {
      //PERFORM ROTATION MOVE, a test//
      running_energy += Rotation_Move(Box, SystemComponents, d_a, Mol, NewMol, FF, y, dUdlambda, Random, SelectedMolInComponent, comp);
      //printf("After Translation: running energy: %.10f\n", running_energy);
    }
    else if(RANDOMNUMBER < SystemComponents.Moves[comp].WidomProb + SystemComponents.Moves[comp].TranslationProb+SystemComponents.Moves[comp].RotationProb)
    {
      WidomCount ++;
      //printf("Performing Widom\n");
      size_t SelectedTrial=0; bool SuccessConstruction = false; double energy = 0.0; double StoredR = 0.0;
      double Rosenbluth=Widom_Move_FirstBead(Box, SystemComponents, d_a, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, true, false, false, StoredR, &SelectedTrial, &SuccessConstruction, &energy, false); //first false: Reinsertion? second false: Retrace? third false is for using Dual-Precision. For Widom Insertion, don't use it.//
      if(SystemComponents.Moleculesize[comp] > 1 && Rosenbluth > 1e-150)
      {
        size_t SelectedFirstBeadTrial = SelectedTrial; double temp_energy = energy;
        Rosenbluth*=Widom_Move_Chain(Box, SystemComponents, d_a, Mol, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, true, false, &SelectedTrial, &SuccessConstruction, &energy, SelectedFirstBeadTrial, false); //false is for using Dual-Precision. For Widom Insertion, don't use it.//
      }
      //Assume 5 blocks
      size_t BlockIDX = i/(Cycles/SystemComponents.Moves[comp].NumberOfBlocks); //printf("BlockIDX=%zu\n", BlockIDX);
      Widom.Rosenbluth[BlockIDX]+= Rosenbluth;
      Widom.RosenbluthSquared[BlockIDX]+= Rosenbluth*Rosenbluth;
      Widom.RosenbluthCount[BlockIDX]++;
    }
    else if(RANDOMNUMBER < SystemComponents.Moves[comp].ReinsertionProb + SystemComponents.Moves[comp].WidomProb + SystemComponents.Moves[comp].TranslationProb+SystemComponents.Moves[comp].RotationProb)
    {
      if(DEBUG) printf("Before Reinsertion, energy: %.10f\n", running_energy);
      running_energy += Reinsertion(Box, SystemComponents, d_a, Mol, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, DualPrecision);
    }
    else
    {
      // DO GCMC INSERTION //
      if(get_random_from_zero_to_one() < 0.5){ //0.5){
        running_energy += Insertion(Box, SystemComponents, d_a, Mol, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, DualPrecision);}
      else{
        if(DEBUG){printf("Cycle: %zu, DOING DELETION\n", i);}
        running_energy += Deletion(Box, SystemComponents, d_a, Mol, NewMol, FF, Random, Widom, SelectedMolInComponent, comp, DualPrecision);}
    }
    if(i%500==0 &&(SystemComponents.Moves[comp].TranslationTotal > 0))
    {
      printf("i: %zu\n", i);
      Update_Max_Translation(FF, SystemComponents.Moves[comp]);
    }
    if(DEBUG)
    {
      printf("After %zu MOVE: Sum energies\n", i);
      double* xxx; xxx = (double*) malloc(sizeof(double)*2);
      double* device_xxx = CUDA_copy_allocate_double_array(xxx, 2);
      one_thread_GPU_test<<<1,1>>>(Box, d_a, FF, device_xxx); cudaMemcpy(xxx, device_xxx, sizeof(double), cudaMemcpyDeviceToHost);
      printf("Current Total Energy (1 thread GPU): %.10f, running total: %.10f\n", xxx[0], init_energy+running_energy);
      cudaDeviceSynchronize();
      if(abs(xxx[0] - (init_energy+running_energy)) > 0.1) //means that there is an energy drift
      {
        printf("THere is an energy drift at cycle %zu\n", i);
      }
      cudaFree(device_xxx);
    }
  }
  //print statistics
  for(size_t comp = SystemComponents.NumberOfFrameworks; comp < SystemComponents.Total_Components; comp++)
  {
    Print_Translation_Statistics(SystemComponents.Moves[comp], FF);
    Print_Widom_Statistics(Widom, SystemComponents.Moves[comp], FF, Constants.energy_to_kelvin);
    Print_Swap_Statistics(SystemComponents.Moves[comp]);
    printf("TransCount: %zu\n", transCount);
    printf("total-deltaU: %.10f\n", running_energy);
  }
  return running_energy;
}