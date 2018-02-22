
#define GLOBAL_VARIABLES_ARE_BAD
#include <rifine.hh>
#undef GLOBAL_VARIABLES_ARE_BAD

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
// #include <boost/random/mersenne_twister.hpp>

#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <devel/init.hh>
#include <riflib/RotamerGenerator.hh>
#include <riflib/rosetta_field.hh>
#include <riflib/util.hh>
#include <riflib/rotamer_energy_tables.hh>

// #include <numeric/alignment/QCP_Kernel.hh>
#include <parallel/algorithm>
#include <exception>
#include <stdexcept>

#include <scheme/actor/Atom.hh>
#include <scheme/actor/BackboneActor.hh>
#include <scheme/actor/VoxelActor.hh>
#include <scheme/kinematics/Director.hh>
#include <scheme/kinematics/SceneBase.hh>
#include <scheme/nest/pmap/OriTransMap.hh>
#include <scheme/numeric/rand_xform.hh>
// #include <scheme/objective/ObjectiveFunction.hh>
#include <scheme/objective/voxel/FieldCache.hh>
// #include <scheme/objective/voxel/VoxelArray.hh>
// #include <scheme/objective/hash/XformMap.hh>
// #include <scheme/objective/storage/RotamerScores.hh>
#include <scheme/util/StoragePolicy.hh>
#include <scheme/search/HackPack.hh>
#include <scheme/objective/integration/SceneObjective.hh>

#include <riflib/RifFactory.hh>

// rifdock subroutine util.hh
#include <riflib/rifdock_typedefs_rifine.hh>
#include <riflib/rifdock_subroutines/util_rifine.hh>
#include <riflib/scaffold/ScaffoldDataCache.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <chrono>
#include <random>
#include <random>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>


int main( int argc, char *argv[] )
{
   
    
    
    using namespace core::scoring;
    using namespace devel::scheme;
    typedef numeric::xyzVector<core::Real> Vec;
    typedef numeric::xyzMatrix<core::Real> Mat;
    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;
    using devel::scheme::print_header;
    using devel::scheme::RotamerIndex;
    
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////// static shit ??
    /////////////////////////////////////////////////////////////////////////////////
    typedef ::scheme::util::SimpleArray< 3, float > F3;
    typedef ::scheme::util::SimpleArray< 3, int > I3;
    
    
    
    register_options();
    devel::init(argc, argv);
    print_header("setup global options");
    RifDockOpt opt;
    opt.init_from_cli();
    
    
    
    
    for ( int iscaff = 0; iscaff < opt.scaffold_fnames.size(); ++iscaff )
    {
        std::string scaff_fname = opt.scaffold_fnames.at( iscaff );
        std::vector< std::string > scaffold_sequence_glob0; // Scaffold sequence in name3 space
        utility::vector1< core::Size > scaffold_res; // seqpose of residue to design, default whold scaffold
        
        
            // the docking process happens here!!!
        
            std::string scafftag = utility::file_basename( utility::file::file_basename( scaff_fname ) );
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "//////   begin scaffold " << scafftag << " " << iscaff << " of " << opt.scaffold_fnames.size() << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            std::cout << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
            
            
            core::pose::Pose scaffold, scaffold_centered, scaffold_full_centered;
            float scaff_radius = 0.0;
            float redundancy_filter_rg = 0.0;
            
            std::vector< int > scaffres_g2l, scaffres_l2g;
            std::vector< bool > scaffuseres;
            Eigen::Vector3f scaffold_center;
            EigenXform scaffold_perturb = EigenXform::Identity();
            
            // scaffold initialization, scaffold residues, scaffold_twobody, onebody /////////////////////////
            {
                core::import_pose::pose_from_file( scaffold, scaff_fname );
                scaffold_full_centered = scaffold;
                
                for ( int ir = 1; ir <= scaffold.size(); ++ir ) {
                    scaffold_sequence_glob0.push_back( scaffold.residue(ir).name3() );
                }
                
                // scaffold residues
                std::string scaff_res_fname = "";
                if ( opt.scaffold_res_fnames.size() ) {
                    if ( opt.scaffold_res_fnames.size() == opt.scaffold_fnames.size() ) {
                        scaff_res_fname = opt.scaffold_res_fnames.at( iscaff );
                    } else if ( opt.scaffold_res_fnames.size() == 1 ) {
                        scaff_res_fname = opt.scaffold_res_fnames.front();
                    } else {
                        utility_exit_with_message( "-scaffold_res list not the same length as -scaffolds list" );
                    }
                    if ( opt.scaffold_res_use_best_guess ) {
                        utility_exit_with_message( "should only use -scaffold_res_use_best_guess true iff not specifying scaffold_res" );
                    }
                    scaffold_res = devel::scheme::get_res( scaff_res_fname, scaffold );
                } else if ( opt.scaffold_res_use_best_guess ) {
                    scaffold_res = devel::scheme::get_designable_positions_best_guess( scaffold, opt.dont_use_scaffold_loops );
                    std::cout << "using scaffold residues: ";
                    for ( auto ir : scaffold_res ) std::cout << " " << ir << scaffold.residue(ir).name3();
                    std::cout << std::endl;
                } else {
                    for( int ir = 1; ir <= scaffold.size(); ++ir){
                        if( !scaffold.residue(ir).is_protein() ) continue;
                        if( scaffold.residue(ir).name3() == "PRO" ) continue;
                        if( scaffold.residue(ir).name3() == "GLY" ) continue;
                        if( scaffold.residue(ir).name3() == "CYS" ) continue;
                        scaffold_res.push_back(ir);
                    }
                } // scaffold residues selection
                
                if ( opt.scaff2ala )              ::devel::scheme::pose_to_ala( scaffold );
                else if ( opt.scaff2alaselonly )  ::devel::scheme::pose_to_ala( scaffold, scaffold_res );
                std::cout << "rifdock scaffold_res: " << scaffold_res << std::endl;
                
                             
                
                // scaffold center and move the scaffold to the center. ########## setup the pose..
                scaffold_center = pose_center( scaffold, scaffold_res );
                //scaffold_centered = scaffold;
                // scaffold_full_centered = scaffold;
                
                for (int ir = 1; ir <= scaffold.size(); ++ir ) {
                    Vec tmp( scaffold_center[0], scaffold_center[1], scaffold_center[2] );
                    //for ( int ia = 1; ia <= scaffold.residue_type(ir).natoms(); ++ia) {
                    //    core::id::AtomID aid(ia, ir);
                    //    scaffold_centered.set_xyz( aid, scaffold.xyz(aid) - tmp );
                    //}
                    for( int ia = 1; ia <= scaffold_full_centered.residue_type(ir).natoms(); ++ia ){
                        core::id::AtomID aid(ia,ir);
                        scaffold_full_centered.set_xyz( aid, scaffold_full_centered.xyz(aid) - tmp );
                    }
                }
                
                //add_pdbinfo_if_missing( scaffold );
                //add_pdbinfo_if_missing( scaffold_centered );
                //add_pdbinfo_if_missing( scaffold_full_centered );
                
                
                scaffold_full_centered.dump_pdb(scafftag + ".pdb");
            }
    }
    
    
 
    
    
    std::cout << "scaffold_preparation_DONE" << std::endl;
    
    return 0;
}
























