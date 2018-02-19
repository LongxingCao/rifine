// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_output_results_hh
#define INCLUDED_riflib_rifdock_subroutines_output_results_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


namespace devel {
namespace scheme {

void
dump_rif_result(
    RifDockData & rdd,
    RifDockResult const & selected_result, 
    std::string const & pdboutfile, 
    int iresl,
    bool quiet = true,
    std::string const & resfileoutfile = "",
    std::string const & allrifrotsoutfile = ""
    ) {


    using ObjexxFCL::format::F;
    using ObjexxFCL::format::I;

    using namespace devel::scheme;
    using std::cout;
    using std::endl;


    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

    ScaffoldIndex si = selected_result.scene_index.scaffold_index;
    ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow( si );
    std::vector<int> const & scaffres_g2l = *(sdc->scaffres_g2l_p);
    std::vector<int> const & scaffres_l2g = *(sdc->scaffres_l2g_p);
    std::vector<std::vector<float> > const & scaffold_onebody_glob0 = *(sdc->scaffold_onebody_glob0_p);
    uint64_t const scaffold_size = scaffres_g2l.size();


    core::pose::Pose pose_from_rif;

    if ( rdd.opt.output_full_scaffold ) {        sdc->setup_both_full_pose( rdd.target ); pose_from_rif = *(sdc->mpc_both_full_pose.get_pose());
    } else if( rdd.opt.output_scaffold_only ) {                                         pose_from_rif = *(sdc->scaffold_centered_p);
    } else if( rdd.opt.output_full_scaffold_only ) {                                    pose_from_rif = *(sdc->scaffold_full_centered_p);
    } else {                                        sdc->setup_both_pose( rdd.target ); pose_from_rif = *(sdc->mpc_both_pose.get_pose());
    }


    rdd.director->set_scene( selected_result.scene_index, iresl, *rdd.scene_minimal );

    devel::scheme::ScoreRotamerVsTarget<
        VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
    > rot_tgt_scorer;
    rot_tgt_scorer.rot_index_p_ = rdd.rot_index_p;
    rot_tgt_scorer.target_field_by_atype_ = rdd.target_field_by_atype;
    rot_tgt_scorer.target_donors_ = *rdd.target_donors;
    rot_tgt_scorer.target_acceptors_ = *rdd.target_acceptors;
    rot_tgt_scorer.hbond_weight_ = rdd.packopts.hbond_weight;
    rot_tgt_scorer.upweight_iface_ = rdd.packopts.upweight_iface;
    rot_tgt_scorer.upweight_multi_hbond_ = rdd.packopts.upweight_multi_hbond;



    EigenXform xposition1 = rdd.scene_minimal->position(1);
    EigenXform xalignout = EigenXform::Identity();
    if( rdd.opt.align_to_scaffold ){
        xalignout = xposition1.inverse();
    }

    xform_pose( pose_from_rif, eigen2xyz(xalignout)           ,        scaffold_size+1, pose_from_rif.size());
    xform_pose( pose_from_rif, eigen2xyz(xalignout*xposition1),                      1,        scaffold_size );


    std::vector< std::pair< int, std::string > > brians_infolabels;

    std::ostringstream packout, allout;
    std::map< int, std::string > pikaa;
    int chain_no = pose_from_rif.num_chains();   
    int res_num = pose_from_rif.size() + 1;
    const std::string chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for( int i_actor = 0; i_actor < rdd.scene_minimal->template num_actors<BBActor>(1); ++i_actor ){
        BBActor bba = rdd.scene_minimal->template get_actor<BBActor>(1,i_actor);
        int const ires = scaffres_l2g.at( bba.index_ );


        {
            std::vector< std::pair< float, int > > rotscores;
            rdd.rif_ptrs[iresl]->get_rotamers_for_xform( bba.position(), rotscores );
            typedef std::pair<float,int> PairFI;
            BOOST_FOREACH( PairFI const & p, rotscores ){
                int const irot = p.second;
                float const onebody = scaffold_onebody_glob0.at( ires ).at( irot );
                float const sc = p.first + onebody;
                float const rescore = rot_tgt_scorer.score_rotamer_v_target( irot, bba.position(), 10.0, 4 );
                if( sc < 0 || rescore + onebody < 0  || p.first + onebody < 0){
                    if ( ! rdd.opt.rif_rots_as_chains) {
                        allout << "MODEL" << endl;
                    }
                    BOOST_FOREACH( SchemeAtom a, rdd.rot_index_p->rotamers_.at( irot ).atoms_ ){
                        a.set_position( xalignout * bba.position() * a.position() ); // is copy
                        a.nonconst_data().resnum = rdd.opt.rif_rots_as_chains ? res_num : ires;
                        a.nonconst_data().chain = rdd.opt.rif_rots_as_chains ? chains.at( chain_no % 52 ) : 'A';
                        ::scheme::actor::write_pdb( allout, a, nullptr );
                    }
                    if (! rdd.opt.rif_rots_as_chains) {
                        allout << "ENDMDL" << endl;
                    } else {
                        allout << "TER" << endl;
                        res_num++;
                        chain_no++;
                    }
                    char oneletter = rdd.rot_index_p->oneletter(irot);
                    if( std::find( pikaa[ires+1].begin(), pikaa[ires+1].end(), oneletter ) == pikaa[ires+1].end() ){
                        pikaa[ires+1] += oneletter;
                    }


                    // Brian
                    std::pair< int, int > sat1_sat2 = rdd.rif_ptrs.back()->get_sat1_sat2(bba.position(), irot);

                    if ( ! quiet ) {

                        std::cout << "seqpos:" << I(3, ires+1);
                        std::cout << " " << oneletter;
                        std::cout << " score:" << F(7, 2, sc);
                        std::cout << " irot:" << I(3, irot);
                        std::cout << " 1-body:" << F(7, 2, onebody );
                        std::cout << " rif score:" << F(7, 2, p.first);
                        std::cout << " rif rescore:" << F(7, 2, rescore);
                        std::cout << " sats:" << I(3, sat1_sat2.first) << " " << I(3, sat1_sat2.second);
                        std::cout << std::endl;
                    }

                    if (sat1_sat2.first > -1) {
                        std::pair< int, std::string > brian_pair;
                        brian_pair.first = ires + 1;
                        brian_pair.second = "HOT_IN:" + str(sat1_sat2.first);
                        brians_infolabels.push_back(brian_pair);
                    }

                }

            }
        }

    }



    // Actually place the rotamers on the pose
    core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
    std::ostringstream resfile, expdb;
    resfile << "ALLAA" << std::endl;
    resfile << "start" << std::endl;
    expdb << "rif_residues ";


    sanity_check_hackpack( rdd, selected_result.scene_index, selected_result.rotamers_, rdd.scene_pt.front());

    std::vector<int> needs_RIFRES;
    for( int ipr = 0; ipr < selected_result.numrots(); ++ipr ){
        int ires = scaffres_l2g.at( selected_result.rotamers().at(ipr).first );
        int irot =                  selected_result.rotamers().at(ipr).second;

        core::conformation::ResidueOP newrsd = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(rdd.rot_index_p->resname(irot)) );
        pose_from_rif.replace_residue( ires+1, *newrsd, true );
        resfile << ires+1 << " A NATRO" << std::endl;
        expdb << ires+1 << (ipr+1<selected_result.numrots()?",":""); // skip comma on last one
        for( int ichi = 0; ichi < rdd.rot_index_p->nchi(irot); ++ichi ){
            pose_from_rif.set_chi( ichi+1, ires+1, rdd.rot_index_p->chi( irot, ichi ) );
        }
        needs_RIFRES.push_back(ires+1);
    }

    // Add PDBInfo labels if they are applicable
    bool using_rosetta_model = selected_result.pose_ != nullptr;
    core::pose::Pose & pose_to_dump( *(selected_result.pose_ ? selected_result.pose_.get() : &pose_from_rif) );
    if( !using_rosetta_model ){
        if( rdd.opt.pdb_info_pikaa ){
            for( auto p : pikaa ){
                std::sort( p.second.begin(), p.second.end() );
                pose_to_dump.pdb_info()->add_reslabel(p.first, "PIKAA" );
                pose_to_dump.pdb_info()->add_reslabel(p.first, p.second );
            }
        } else {
            std::sort(needs_RIFRES.begin(), needs_RIFRES.end());
            for( int seq_pos : needs_RIFRES ){
                pose_to_dump.pdb_info()->add_reslabel(seq_pos, "RIFRES" );
            }
        }
    }

    for ( auto p : brians_infolabels ) {
        pose_to_dump.pdb_info()->add_reslabel(p.first, p.second);
    }

    rdd.scaffold_provider->modify_pose_for_output(si, pose_to_dump);


    // Dump the main output
    utility::io::ozstream out1( pdboutfile );
    out1 << expdb.str() << std::endl;
    pose_to_dump.dump_pdb(out1);
    if ( rdd.opt.dump_all_rif_rots_into_output ) {
        if ( rdd.opt.rif_rots_as_chains ) out1 << "TER" << endl;
        out1 << allout.str();
    }
    out1.close();

    // Dump a resfile
    if( rdd.opt.dump_resfile ){
        utility::io::ozstream out1res( resfileoutfile );
        out1res << resfile.str();
        out1res.close();
    }

    // Dump the rif rots
    if( rdd.opt.dump_all_rif_rots ){
        utility::io::ozstream out2( allrifrotsoutfile );
        out2 << allout.str();
        out2.close();
    }

} // end crappy pdb io


void
dump_search_point(
    RifDockData & rdd,
    SearchPoint const & search_point, 
    std::string const & pdboutfile, 
    int iresl,
    bool quiet) {


    RifDockResult result;
    result.scene_index = search_point.index;

    dump_rif_result( rdd, result, pdboutfile, iresl, quiet, "", "" );
}




void
output_results(
    std::vector< RifDockResult > & selected_results,
    RifDockData & rdd,
    utility::io::ozstream & dokout,
    int64_t npack) {


    using namespace core::scoring;
        using std::cout;
        using std::endl;
        using namespace devel::scheme;
        typedef numeric::xyzVector<core::Real> Vec;
        typedef numeric::xyzMatrix<core::Real> Mat;
        // typedef numeric::xyzTransform<core::Real> Xform;
        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;
        using devel::scheme::print_header;
        using ::devel::scheme::RotamerIndex;

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;

    typedef _RifDockResult<DirectorBase> RifDockResult;
    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;

    if( rdd.opt.align_to_scaffold ) std::cout << "ALIGN TO SCAFFOLD" << std::endl;
    else                        std::cout << "ALIGN TO TARGET"   << std::endl;
    for( int i_selected_result = 0; i_selected_result < selected_results.size(); ++i_selected_result ){
        RifDockResult const & selected_result = selected_results.at( i_selected_result );

// Brian Injection
        ScaffoldIndex si = selected_results[i_selected_result].scene_index.scaffold_index;
        ScaffoldDataCacheOP sdc = rdd.scaffold_provider->get_data_cache_slow( si );


        std::string const & scafftag = sdc->scafftag;

/////
        std::string pdboutfile = rdd.opt.outdir + "/" + scafftag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        if( rdd.opt.output_tag.size() ){
            pdboutfile = rdd.opt.outdir + "/" + scafftag+"_" + rdd.opt.output_tag + "_" + devel::scheme::str(i_selected_result,9)+".pdb.gz";
        }

        std::string resfileoutfile = rdd.opt.outdir + "/" + scafftag+"_"+devel::scheme::str(i_selected_result,9)+".resfile";
        std::string allrifrotsoutfile = rdd.opt.outdir + "/" + scafftag+"_allrifrots_"+devel::scheme::str(i_selected_result,9)+".pdb.gz";

        std::ostringstream oss;
        oss << "rif score: " << I(4,i_selected_result)
            << " rank "       << I(9,selected_result.isamp)
            << " dist0:    "  << F(7,2,selected_result.dist0)
            << " packscore: " << F(7,3,selected_result.packscore)
            // << " score: "     << F(7,3,selected_result.nopackscore)
            // << " rif: "       << F(7,3,selected_result.rifscore)
            << " steric: "    << F(7,3,selected_result.stericscore)
            << " cluster: "   << I(7,selected_result.cluster_score)
            << " rifrank: "   << I(7,selected_result.prepack_rank) << " " << F(7,5,(float)selected_result.prepack_rank/(float)npack)
            << " " << pdboutfile
            << std::endl;
        std::cout << oss.str();
        dokout << oss.str(); dokout.flush();


        dump_rif_result(rdd, selected_result, pdboutfile, rdd.RESLS.size()-1, false, resfileoutfile, allrifrotsoutfile);



    }



}




}}




#endif