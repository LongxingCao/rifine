// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/RifGeneratorUserHotspots.hh>
    
    #include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/packing.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>

	#include <ObjexxFCL/format.hh>

	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

	#include <core/id/AtomID.hh>
	#include <core/pose/Pose.hh>
	#include <core/scoring/motif/util.hh>
	#include <core/import_pose/import_pose.hh>

  	#include <numeric/xyzMatrix.hh>
	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>
	#include <scheme/numeric/rand_xform.hh>
	#include <scheme/actor/Atom.hh>
	#include <scheme/actor/BackboneActor.hh>

	#include <Eigen/SVD>
	#include <Eigen/Core>
	#include <Eigen/Geometry>

	#include <utility/io/ozstream.hh>



	#include <scheme/objective/hash/XformMap.hh>
	#include <scheme/objective/storage/RotamerScores.hh>
	#include <scheme/actor/BackboneActor.hh>


namespace devel {
namespace scheme {
namespace rif {

	void
	RifGeneratorUserHotspots::modify_rotamer_spec(
		::scheme::chemical::RotamerIndexSpec& rot_spec
		
	){
		for( int i_hotspot_group = 0; i_hotspot_group < this->opts.hotspot_files.size(); ++i_hotspot_group ){
			
			std::string const & hotspot_file = this->opts.hotspot_files[i_hotspot_group];
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose,hotspot_file);

			
			for( int i_hspot_res = 1; i_hspot_res <= pose.size(); ++i_hspot_res ){
				std::string resn = pose.residue(i_hspot_res).name3();
				int parent_key = -1;
				
				int n_proton_chi = 0;
				if (resn == "CYS" || resn == "SER" || resn == "THR" || resn == "TYR"){
					n_proton_chi = 1;
				}

				std::vector<float> mychi(pose.residue(i_hspot_res).nchi());
				for (int n_chi = 0; n_chi < mychi.size(); n_chi++){
					mychi.at(n_chi) = pose.chi(n_chi+1,i_hspot_res);
				}
				
				//check if the rotamer exist in rot_spec
				bool add_this_rotamer = false;
				for( int irot = 0; irot < rot_spec.size(); ++irot ){

					if( rot_spec.resname(irot) != resn ) continue;
					bool duplicate_rotamer = true;
					for (int n_chi = 0; n_chi < rot_spec.get_rotspec(irot).chi_.size(); ++n_chi){
						//std::cout << "checking " << rot_spec.get_rotspec(irot).resname_ << " " << rot_spec.get_rotspec(irot).chi_.at(n_chi) << " " << resn <<" "<< mychi.at(n_chi) << std::endl;
						duplicate_rotamer &= ::scheme::chemical::impl::angle_is_close( rot_spec.get_rotspec(irot).chi_.at(n_chi), mychi.at(n_chi), 5.0f );
					}

					if (duplicate_rotamer){
						std::cout << "duplicated rotamer, not adding: " << i_hotspot_group << " " << i_hspot_res << " " << resn << std::endl;
						add_this_rotamer = false;
						break;
					} 
					else if (!duplicate_rotamer){
						add_this_rotamer = true;
					}					
				}// end loop over all existing rot_spec rotamers

				if (add_this_rotamer){
					std::cout << "Adding input rotamers: " << i_hspot_res << " " << resn << std::endl;
					rot_spec.add_rotamer(resn,mychi,n_proton_chi,parent_key);
				}
			}//end loop over all hotspot res within one hotspot file
		}// end loop over all hotspot files
		//rot_spec.load();
		// utility::io::ozstream outfile("test_out.text");
		// rot_index_spec.save(outfile);
		// rot_index_spec.fill_rotamer_index(*(params->rot_index_p));
		std::cout << "modify_rotamer_sepc DONE " << std::endl;
	}// end modify_rotamer_spec


	void
	RifGeneratorUserHotspots::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){

		typedef ::scheme::actor::BackboneActor<EigenXform> BBActor;

		typedef ::Eigen::Matrix<float,3,1> Pos;
	
		// some sanity checks
		int const n_hspot_groups = this->opts.hotspot_files.size();
		runtime_assert_msg( n_hspot_groups, "no hotspot group files specified!!" );
		runtime_assert_msg( n_hspot_groups<16, "too many hotspot groups!!" );

		std::cout << "this RIF type doesn't support sat groups!!!" << std::endl;

		std::cout << "RifGeneratorUserHotspots opts:" << std::endl;
		std::cout << "    hotspot_sample_cart_bound:  " << this->opts.hotspot_sample_cart_bound;
		std::cout << "    hotspot_sample_angle_bound: " << this->opts.hotspot_sample_angle_bound;
		std::cout << "    hbond_weight:               " << this->opts.hbond_weight;
		std::cout << "    upweight_multi_hbond:       " << this->opts.upweight_multi_hbond;
		std::cout << "    target_center:              "
			<< this->opts.target_center[0] << " "
			<< this->opts.target_center[1] << " "
			<< this->opts.target_center[2] << std::endl;
		//translation to apply to input hotspot files
		Pos target_vec;
		numeric::xyzVector<double> xyz_tgt_cen( this->opts.target_center[0], this->opts.target_center[1], this->opts.target_center[2] );
		target_vec << opts.target_center[0], opts.target_center[1],opts.target_center[2];


		for( auto s : this->opts.hotspot_files ){
			std::cout << "    hotspot_group:              " << s << std::endl;
		}


		// setup the hacky but fast scorer
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer;
		{
			std::vector< ::scheme::chemical::HBondRay > target_donors, target_acceptors;
			for( auto ir : params->target_res ){
				::devel::scheme::get_donor_rays   ( *params->target, ir, params->hbopt, target_donors );
				::devel::scheme::get_acceptor_rays( *params->target, ir, params->hbopt, target_acceptors );
			}
			std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
			{
				rot_tgt_scorer.rot_index_p_ = params->rot_index_p;
				rot_tgt_scorer.target_field_by_atype_ = params->field_by_atype;
				rot_tgt_scorer.target_donors_ = target_donors;
				rot_tgt_scorer.target_acceptors_ = target_acceptors;
				rot_tgt_scorer.hbond_weight_ = this->opts.hbond_weight;
				rot_tgt_scorer.upweight_multi_hbond_ = this->opts.upweight_multi_hbond;
				rot_tgt_scorer.upweight_iface_ = 1.0;

			}
		}

		


    	int const NSAMP = this->opts.hotspot_nsamples;
    	if (NSAMP > opts.dump_hotspot_samples && opts.dump_hotspot_samples > 0) utility_exit_with_message("too many NSAMP");


    	std::mt19937 rng((unsigned int)time(0) + 934875);
    	float const radius_bound = this->opts.hotspot_sample_cart_bound;
    	float const degrees_bound = this->opts.hotspot_sample_angle_bound;
    	float const radians_bound = degrees_bound * M_PI/180.0;


		std::ofstream hotspot_dump_file;
		std::ostringstream os;
		os << "hotspots.pdb";
		std::string s = os.str();
		if (opts.dump_hotspot_samples>=NSAMP){		
			hotspot_dump_file.open (s, std::fstream::in | std::fstream::out | std::fstream::app);
			core::pose::Pose target_pose(*params->target);
			target_pose.dump_pdb(hotspot_dump_file);
		}

    	// std::ostream & out( std::cout );
    	// std::ofstream out;
    	// out.open("rifgen.txt");

		// loop over files (one file is one hotspot group)

		print_header("Building RIF from resampled user hotspots");
		for( int i_hotspot_group = 0; i_hotspot_group < this->opts.hotspot_files.size(); ++i_hotspot_group ){

			std::string const & hotspot_file = this->opts.hotspot_files[i_hotspot_group];
			std::cout << "Hotspot group: " << i_hotspot_group << " - " << hotspot_file << std::endl;
			
			// read hotspot file into pose
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose,hotspot_file);
      		
      		std::cout << "Processing hotspots... " << std::flush; // No endl here!!!!
			// read in pdb files # i_hotspot_group
			for( int i_hspot_res = 1; i_hspot_res <= pose.size(); ++i_hspot_res ){

				std::cout << i_hspot_res << " " << std::flush; // No endl here!!!!

				int input_nheavy = pose.residue(i_hspot_res).nheavyatoms();
				

				EigenXform Xref = ::scheme::chemical::make_stub<EigenXform>(
		            pose.residue(i_hspot_res).xyz( input_nheavy - 2 ) - xyz_tgt_cen,
		            pose.residue(i_hspot_res).xyz( input_nheavy - 1 ) - xyz_tgt_cen,
		            pose.residue(i_hspot_res).xyz( input_nheavy - 0 ) - xyz_tgt_cen
		        );
				
				//get last atom in hotspot residue and iterate over last 3 
      			core::conformation::Atoms::const_iterator iter = pose.residue(i_hspot_res).heavyAtoms_end();
				core::conformation::Atoms::const_iterator end = iter-3;
      	
				//these are the hotspot atoms				
				Pos hot_atom1; Pos hot_atom2; Pos hot_atom3;
				
				//iterate and save xyz into hotspot 
				while(iter >= end){
        			if (iter == end +2){hot_atom1(0,0) = iter->xyz()[0];hot_atom1(1,0) = iter->xyz()[1];hot_atom1(2,0) = iter->xyz()[2];}
        			if (iter == end +1){hot_atom2(0,0) = iter->xyz()[0];hot_atom2(1,0) = iter->xyz()[1];hot_atom2(2,0) = iter->xyz()[2];}
        			if (iter == end +0){hot_atom3(0,0) = iter->xyz()[0];hot_atom3(1,0) = iter->xyz()[1];hot_atom3(2,0) = iter->xyz()[2];}
        			iter --;
      			}
      			hot_atom1 = hot_atom1 - target_vec; 
				hot_atom2 = hot_atom2 - target_vec; 
				hot_atom3 = hot_atom3 - target_vec;	
				//calculate centroid of hot_spot res and translate with target
				Pos hot_cen = (hot_atom1 + hot_atom2 + hot_atom3)/3;
				

				// for each irot that is the right restype (can be had from rot_intex_p)
				int irot_begin = 0, irot_end = params -> rot_index_p -> size();
				for( int irot = irot_begin; irot < irot_end; ++irot ){

					::Eigen::Matrix<float,3,3> rif_res; // this is the rif residue last three atoms matrix
						
					// assign rif_res position by column
					int hatoms = params -> rot_index_p -> nheavyatoms(irot);
					std::vector<SchemeAtom> const & rotamer_atoms( params->rot_index_p->atoms(irot) );
					EigenXform Xrotamer = ::scheme::chemical::make_stub<EigenXform>(
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 3 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 2 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 1 ).position()
                    );

					if (params -> rot_index_p -> resname(irot) == pose.residue(i_hspot_res).name3())
					{
						
						EigenXform impose; //transform for mapping the Rot to Rif
						::Eigen::Matrix<float,3,3> rif_res; // this is the rif residue last three atoms
						int latoms = params -> rot_index_p -> natoms(irot);
						rif_res << rotamer_atoms[hatoms-1].position()[0],rotamer_atoms[hatoms-2].position()[0],rotamer_atoms[hatoms-3].position()[0],rotamer_atoms[hatoms-1].position()[1],rotamer_atoms[hatoms-2].position()[1],rotamer_atoms[hatoms-3].position()[1],rotamer_atoms[hatoms-1].position()[2],rotamer_atoms[hatoms-2].position()[2],rotamer_atoms[hatoms-3].position()[2];
	     	 			Pos rot_cen = (rif_res.col(0) + rif_res.col(1) + rif_res.col(2))/3;
         				
					
						
	     // 	 			if (!this->opts.single_file_hotspots_insertion){
      //     				//svd superimpose
      //     				//matrix to be SVD decomposed
      //     				::Eigen::Matrix<float,3,3> cov_mtx;										
      //     				cov_mtx = (rif_res.col(0) - rot_cen)*(hot_atom1 - hot_cen).transpose() + (rif_res.col(1) - rot_cen)*(hot_atom2 - hot_cen).transpose() + (rif_res.col(2) - rot_cen)*(hot_atom3 - hot_cen).transpose();
      //       			//Eigen SVD
      //       			::Eigen::JacobiSVD<::Eigen::Matrix<float,3,3>> svd(cov_mtx, Eigen::ComputeFullU | Eigen::ComputeFullV);
						// //Rotation Matrix
						// ::Eigen::Matrix<float,3,3> R_mtx = svd.matrixV() * svd.matrixU().transpose();
						// //Translation Matrix
     	// 				::Eigen::Matrix<float,3,1> T_mtx = - R_mtx*rot_cen + hot_cen;
      //   				//Puting R/T together
      //   				::Eigen::Matrix<float,3,4> Tran_mtx;					
				  // 		Tran_mtx.block<3,3>(0,0) = R_mtx;
      //   				Tran_mtx.block<3,1>(0,3) = T_mtx;
      //   				//Final superimpose EigenXform matrix
						// //EigenXform impose;
						// impose.matrix() = Tran_mtx;
						// }
						
						impose = Xref * Xrotamer.inverse();						
						//Additional matrix definition for manipulation
						//Default Rot starting Xform
						EigenXform x_orig_position = EigenXform::Identity();
						//Xform to move rif to and from starting postion
						EigenXform x_2_orig = EigenXform::Identity();
         				x_2_orig.translation() = -hot_cen;
         				EigenXform x_2_orig_inverse = x_2_orig.inverse();
         				
         				

						//impose = Xref * Xrotamer.inverse();//working
						//impose.matrix() = Tran_mtx;//not working



						//EigenXform O_2_orig = EigenXform::Identity();

						int passes = 1;
						EigenXform O_2_orig = EigenXform::Identity();
         				EigenXform tyr_thing = EigenXform::Identity();	
						if (pose.residue(i_hspot_res).name3() == "TYR") {

							Pos the_axis = (hot_atom1 - hot_atom2).normalized();
							O_2_orig.translation() = -hot_atom1;	
							tyr_thing.rotate( Eigen::AngleAxisf(M_PI, the_axis)); 

							passes = 2;

						} else if (pose.residue(i_hspot_res).name3() == "PHE") {
							
							Pos atom6;
							atom6(0,0) = pose.residue(i_hspot_res).xyz( input_nheavy - 5 )[0];
							atom6(1,0) = pose.residue(i_hspot_res).xyz( input_nheavy - 5 )[1];
							atom6(2,0) = pose.residue(i_hspot_res).xyz( input_nheavy - 5 )[2];
							atom6 = atom6 - target_vec;

							Pos the_axis = (hot_atom1 - atom6).normalized();
							O_2_orig.translation() = -hot_atom1;	
							tyr_thing.rotate( Eigen::AngleAxisf(M_PI, the_axis)); 

							passes = 2;
						}
						

						EigenXform O_2_orig_inverse = O_2_orig.inverse();



						for ( int pass = 0; pass < passes; pass++) {
							//std::cout << Tran_mtx.block<3,3>(0,0)*rif_res+temp << std::endl;
							//for( auto const & x_perturb : sample_position_deltas ){
							int num_of_hotspot_inserted = 0;
							//std::cout << "being parallel block" << std::endl;
							#ifdef USE_OPENMP
							#pragma omp parallel for schedule(dynamic,16)
							#endif

							for(int a = 0; a < NSAMP; ++a){							
								EigenXform x_perturb;
								::scheme::numeric::rand_xform_sphere(rng,x_perturb,radius_bound,radians_bound);

								EigenXform building_x_position = impose * x_orig_position;
								if ( pass == 1 ) {
									building_x_position = O_2_orig_inverse * tyr_thing * O_2_orig * building_x_position;
								}

								EigenXform x_position = x_2_orig_inverse * x_perturb * x_2_orig * building_x_position;
									
								// you can check their "energies" against the target like this, obviously substituting the real rot# and position
								float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target( irot, x_position );
								
								if( positioned_rotamer_score < 5){ // probably want this threshold to be an option or something

									//num_of_hotspot_inserted += 1;
									// EigenXform x_position = x_2_orig_inverse * x_perturb  * x_2_orig * impose * x_orig_position; //test
								    // you can check their "energies" against the target like this, obviously substituting the real rot# and position
									//float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target( irot, x_position );


									//target_pose.dump_pdb(s);
									//myfile.open (s, std::fstream::in | std::fstream::out | std::fstream::app);
									//if( positioned_rotamer_score > 0) {positioned_rotamer_score = -1;}

									//std::cout << positioned_rotamer_score << std::endl;
									//accumulator->insert( x_position, positioned_rotamer_score-100, irot, i_hotspot_group, -1 );
									

										//std::cout << "new :                  " << std::endl;
									auto atom_N = x_position * rotamer_atoms[0].position();
									auto atom_CA = x_position * rotamer_atoms[1].position();
									auto atom_C = x_position * rotamer_atoms[2].position(); 

									BBActor bbact( atom_N, atom_CA, atom_C);
									EigenXform new_x_position = bbact.position();

									accumulator->insert( new_x_position, positioned_rotamer_score-4, irot, 
										this -> opts.single_file_hotspots_insertion ? i_hspot_res : i_hotspot_group,
										 -1 );

								 	if (opts.dump_hotspot_samples>=NSAMP){
								 		hotspot_dump_file <<"MODEL        "<<irot<<a<<"                                                                  \n";
										for( auto a : rotamer_atoms ){
										 	a.set_position( x_position * a.position() );
										 	::scheme::actor::write_pdb(hotspot_dump_file, a, params->rot_index_p->chem_index_ );
										}
										hotspot_dump_file <<"ENDMDL                                                                          \n";							
									} //end dumping hotspot atoms

														

									// myfile.close();
								} // end rotamer insertion score cutoff					
							
							} // end NSAMP
							//std::cout << "this is how many inserted: " << i_hotspot_group << " " << irot << " " << num_of_hotspot_inserted << std::endl;
						
						} // end brian ring flip
					
					} // end loop over rotamers which match hotspot res name
				
				} //  end loop over rotamers library

			} // end loop over hotspot group residue (with in one input pdb)
		
		}// end loop over all hotspot input files

		std::cout << std::endl; // This ends the "progress bar"

		if (opts.dump_hotspot_samples>=NSAMP){
			hotspot_dump_file.close();
		}
		//utility_exit_with_message("done");
		// let the rif builder thing know you're done
		accumulator->checkpoint( std::cout );
	
	} //end RifGeneratorUserHotspots

		//auto rif_ptr = accumulator->rif();
        // std::cout << "testing rifbase key iteration" << std::endl;
        // int count = 0;
        // for( auto key : rif_ptr->key_range() ){
        //     EigenXform bin_center = rif_ptr->get_bin_center(key);
        //     right... the edges of the bins.... this is only *mostly* true
        //     runtime_assert( rif_ptr->get_bin_key(bin_center) == key );
        //     std::cout << "BIN: key " << key << " xform.trans: " << bin_center.translation().transpose() << std::endl;
        //     auto rotscores = rif_ptr->get_rotamers_for_key(key);
        //     runtime_assert( rotscores.size() > 0 );
        //     for( auto rot_score : rotscores ){
        //         // can't wait for cxx17 structured bindings!!!
        //         float rotamer_score = rot_score.first;
        //         int rotamer_number = rot_score.second;
        //         std::string resn = params->rot_index_p->resname(rotamer_number);
        //         std::cout << " rotamer " << rotamer_number << " " << resn << ", score " << rotamer_score << std::endl;
        //     }

        //     //if(++count > 10) utility_exit_with_message("aireost");
        // }

		
		//auto rif_ptr = accumulator->rif();
		//std::cout << "Brian" << std::endl;
}
}
}


